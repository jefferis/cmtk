/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUniformVolumePainter.h>
#include <Base/cmtkFixedVector.h>

#include <IO/cmtkVolumeIO.h>

#include <stdio.h>
#include <memory>

#ifdef CMTK_USE_DCMTK
#  include <dcmtk/dcmdata/dctk.h>
#endif

int Dims[3] = { 256, 256, 256 };

void
SetDims( const char* arg )
{
  sscanf( arg, "%d,%d,%d", &Dims[0], &Dims[1], &Dims[2] );
}

float Delta[3] = { 1.0, 1.0, 1.0 };

void
SetDeltas( const char* arg )
{
  sscanf( arg, "%f,%f,%f", &Delta[0], &Delta[1], &Delta[2] );
}

cmtk::ScalarDataType DataType = cmtk::TYPE_USHORT;
cmtk::Types::DataItem Background = 0;

const char* InputImageName = NULL;
bool InputImageGridOnly = false;

const char* OutputFileName = "phantom.nii";

cmtk::UniformVolumePainter::CoordinateModeEnum CoordinateMode = cmtk::UniformVolumePainter::COORDINATES_INDEXED;

int
doMain( const int argc, const char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Generate phantom image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Generate 3D digital phantom images using a selection of drawing commands" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, mk_phantom_3d "[options] command0 [command1 ...]\n"
		       "\t sphere x,y,z radius value\n"
		       "\t box x0,y0,z0 x1,y1,z1 value\n");

    typedef cmtk::CommandLine::Key Key;
    cl.AddCallback( Key( 'D', "dims" ), SetDims, "Set dimensions in voxels" );
    cl.AddCallback( Key( 'V', "voxel" ), SetDeltas, "Set voxel size in [mm]" );

    cmtk::CommandLine::EnumGroup<cmtk::UniformVolumePainter::CoordinateModeEnum>::SmartPtr modeGroup = cl.AddEnum( "coordinates", &CoordinateMode, "Coordinate specification mode." );
    modeGroup->AddSwitch( Key( "indexed" ), cmtk::UniformVolumePainter::COORDINATES_INDEXED, "Use grid indexes to specify coordinates. For each dimension, the valid value range is [0,Dims-1]." );
    modeGroup->AddSwitch( Key( "absolute" ), cmtk::UniformVolumePainter::COORDINATES_ABSOLUTE, "Use absolute volume coordinates. For each dimension, the valid range is [0,FOV]." );
    modeGroup->AddSwitch( Key( "relative" ), cmtk::UniformVolumePainter::COORDINATES_RELATIVE, "Use relative volume coordinates. For each dimension, the valid range is [0,1]." );
    
    cl.AddSwitch( Key( 'c', "char" ), &DataType, cmtk::TYPE_CHAR, "8 bits, signed" );
    cl.AddSwitch( Key( 'b', "byte" ), &DataType, cmtk::TYPE_BYTE, "8 bits, unsigned" );
    cl.AddSwitch( Key( 's', "short" ), &DataType, cmtk::TYPE_SHORT, "16 bits, signed" );
    cl.AddSwitch( Key( 'u', "ushort" ), &DataType, cmtk::TYPE_USHORT, "16 bits, unsigned" );
    cl.AddSwitch( Key( 'i', "int" ), &DataType, cmtk::TYPE_INT, "32 bits signed" );
    cl.AddSwitch( Key( 'f', "float" ), &DataType, cmtk::TYPE_FLOAT, "32 bits floating point" );
    cl.AddSwitch( Key( 'd', "double" ), &DataType, cmtk::TYPE_DOUBLE, "64 bits floating point\n" );

    cl.AddOption( Key( 'B', "bg" ), &Background, "Set the image background value (use to initialize newly created image)." );

    cl.AddOption( Key( 'I', "import" ), &InputImageName, "Import image" );
    cl.AddOption( Key( "import-grid" ), &InputImageName, "Import image grid only, ignore data", &InputImageGridOnly );    
    cl.AddOption( Key( 'o', "outfile" ), &OutputFileName, "File name for output image" );
    
    cl.Parse( argc, argv );

    cmtk::UniformVolume::SmartPtr volume;
    if ( InputImageName )
      {
      if ( InputImageGridOnly )
	{
	volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( InputImageName ) );
	volume->CreateDataArray( DataType );
	volume->GetData()->Fill( Background );
	}
      else
	volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( InputImageName ) );
      }
    else
      {
      const cmtk::Types::Coordinate size[3] = { Delta[0] * (Dims[0]-1),  Delta[1] * (Dims[1]-1),  Delta[2] * (Dims[2]-1) };
      volume = cmtk::UniformVolume::SmartPtr( new cmtk::UniformVolume( cmtk::UniformVolume::IndexType( Dims ), cmtk::FixedVector<3,cmtk::Types::Coordinate>( size ) ) );
      volume->m_MetaInformation[cmtk::META_SPACE] = volume->m_MetaInformation[cmtk::META_SPACE_ORIGINAL] = cmtk::AnatomicalOrientation::ORIENTATION_STANDARD;
      cmtk::TypedArray::SmartPtr data( cmtk::TypedArray::Create( DataType, volume->GetNumberOfPixels() ) );
      volume->SetData( data );
      data->Fill( Background );
      }
    
    cmtk::UniformVolumePainter painter( volume, CoordinateMode );
    
    const char* nextCmd = cl.GetNextOptional();
    while (  nextCmd )
      {
      if ( ! strcmp( nextCmd, "sphere" ) )
	{
	const char* center = cl.GetNextOptional();
	const char* radius = cl.GetNextOptional();
	const char* value = cl.GetNextOptional();

	float cc[3];
	if ( sscanf( center, "%f,%f,%f", &cc[0], &cc[1], &cc[2] ) != 3 )
	  {
	  cmtk::StdErr << "Parameter 'center' of 'sphere' command must be x,y,z\n";
	  return 1;
	  }
	painter.DrawSphere( cmtk::FixedVector<3,cmtk::Types::Coordinate>( cc ), atof( radius ), atof( value ) );
	}

      if ( ! strcmp( nextCmd, "box" ) )
	{
	const char* fromCorner = cl.GetNextOptional();
	const char* toCorner = cl.GetNextOptional();
	const char* value = cl.GetNextOptional();
	
	float boxFrom[3], boxTo[3];
	if ( sscanf( fromCorner, "%f,%f,%f", &boxFrom[0], &boxFrom[1], &boxFrom[2] ) != 3 )
	  {
	  cmtk::StdErr << "Parameter 'corner0' of 'box' command must be three number x,y,z\n";
	  return 1;
	  }
	if ( sscanf( toCorner, "%f,%f,%f", &boxTo[0], &boxTo[1], &boxTo[2] ) != 3 )
	  {
	  cmtk::StdErr << "Parameter 'corner1' of 'box' command must be three numbers x,y,z\n";
	  return 1;
	  }
	painter.DrawBox( cmtk::FixedVector<3,cmtk::Types::Coordinate>( boxFrom ), cmtk::FixedVector<3,cmtk::Types::Coordinate>( boxTo ), atof( value ) );
	}

#ifdef CMTK_USE_DCMTK
      if ( ! strcmp( nextCmd, "mrs-voxel" ) )
	{
	const char* dicom = cl.GetNextOptional();
	const char* value = cl.GetNextOptional();
	
	std::auto_ptr<DcmFileFormat> fileformat( new DcmFileFormat );
	
	fileformat->transferInit();
	OFCondition status = fileformat->loadFile( dicom );
	fileformat->transferEnd();
	
	if ( !status.good() ) 
	  {
	  cmtk::StdErr << "Error: cannot read DICOM file " << dicom << " (" << status.text() << ")\n";
	  exit( 1 );
	  }
	
	DcmDataset *dataset = fileformat->getDataset();
	if ( ! dataset )
	  {
	  exit( 1 );
	  }
	
	Float64 cntr[3];
	Float64 size[3];

	const DcmTagKey keyVoxelLocation( 0x0043, 0x108c); // private GE tag "Voxel Location"
	if ( dataset->findAndGetFloat64( keyVoxelLocation, cntr[0], 0 ).good() )
	  {
	  cntr[0] *= -1; // convert RL to LR
	  dataset->findAndGetFloat64( keyVoxelLocation, cntr[1], 1 ); cntr[1] *= -1; // convert AP to PA
	  dataset->findAndGetFloat64( keyVoxelLocation, cntr[2], 2 );
	  
	  dataset->findAndGetFloat64( keyVoxelLocation, size[0], 3 );
	  dataset->findAndGetFloat64( keyVoxelLocation, size[1], 4 );
	  dataset->findAndGetFloat64( keyVoxelLocation, size[2], 5 );
	  }
	else
	  {
	  // Some image GE MRS DICOM files do not have "Voxel Location" tag, but store location in different undocumented tags.
	  // Order of x/y is switched, and some values have negative signs.
	  dataset->findAndGetFloat64( DcmTagKey( 0x0019, 0x10b3 ), cntr[0] ); // already * -1
	  dataset->findAndGetFloat64( DcmTagKey( 0x0019, 0x10b2 ), cntr[1] ); cntr[1] *= -1;
	  dataset->findAndGetFloat64( DcmTagKey( 0x0019, 0x10b4 ), cntr[2] ); 

	  dataset->findAndGetFloat64( DcmTagKey( 0x0019, 0x10b0 ), size[0] ); 
	  dataset->findAndGetFloat64( DcmTagKey( 0x0019, 0x10af ), size[1] );
	  dataset->findAndGetFloat64( DcmTagKey( 0x0019, 0x10b1 ), size[2] );
	  }

	cmtk::UniformVolumePainter roiPainter( volume, cmtk::UniformVolumePainter::COORDINATES_INDEXED );

	cmtk::FixedVector<3,cmtk::Types::Coordinate> roiCntr( cntr );
	cmtk::FixedVector<3,cmtk::Types::Coordinate> roiSize( size );
	roiSize *= 0.5;

	roiPainter.DrawBox( volume->PhysicalToIndex( roiCntr - roiSize ), volume->PhysicalToIndex( roiCntr + roiSize ), atof( value ) );
	}
#endif // #ifdef CMTK_USE_DCMTK

      nextCmd = cl.GetNextOptional();
      }
    
    cmtk::VolumeIO::Write( *volume, OutputFileName );
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e;
    return 1;
    }
  
  return 0;
}

#include "cmtkSafeMain"
