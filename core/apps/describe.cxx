/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>

#include <cmtkUniformVolume.h>
#include <cmtkFileFormat.h>
#include <cmtkVolumeIO.h>
#include <cmtkVector3D.h>

#include <stdio.h>
#include <list>

#ifdef CMTK_BUILD_MPI
#    include <mpi.h>
#endif

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace describe
{
#endif

bool Verbose = false;

const char* ReadOrientation = NULL;

bool PrintCOM = false;
bool PrintFOM = false;

bool MachineReadable = false;

std::list<cmtk::Vector3D> ProbeListIndex;

const char*
CallbackProbeIndex( const char* arg )
{
  int x, y, z;
  if ( 3 == sscanf( arg, "%d,%d,%d", &x, &y, &z ) )
    {
    ProbeListIndex.push_back( cmtk::Vector3D( x, y, z ) );
    }
  return NULL;
}

int
main( int argc, char *argv[] )
{
#ifdef CMTK_BUILD_MPI
  MPI::Init( argc, argv );
  const int mpiRank = MPI::COMM_WORLD.Get_rank();
  const int mpiSize = MPI::COMM_WORLD.Get_size();
#endif

  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Describe image file formats and parameters" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool prints a detailed description of the input image(s)" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] imageFile0 [imageFile1 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Be verbose" );

    cl.AddSwitch( Key( 'm', "machine-readable" ), &MachineReadable, true, "Print output in format that is easy to parse automatically." );

    cl.AddSwitch( Key( "read-ras" ), &ReadOrientation, "RAS", "Read image in RAS orientation" );
    
    cl.AddSwitch( Key( "print-com" ), &PrintCOM, true, "Compute and print center of mass." );
    cl.AddSwitch( Key( "print-fom" ), &PrintFOM, true, "Compute and print first-order moments." );

    cl.AddCallback( Key( "probe-index" ), CallbackProbeIndex, "Add pixel index for probing." );

    cl.Parse();

    for ( const char* next = cl.GetNext(); next; next = cl.GetNextOptional() ) 
      {
      cmtk::Volume::SmartPtr volume;
      if ( ReadOrientation )
	volume = cmtk::Volume::SmartPtr( cmtk::VolumeIO::ReadOriented( next, ReadOrientation ) );
      else
 	volume = cmtk::Volume::SmartPtr( cmtk::VolumeIO::Read( next ) );
      
      if ( ! volume ) 
	{
	fputs( "Could not read volume data.\n", stderr );
	continue;
	}
      
      const char* orientOriginal = volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION_ORIGINAL].c_str();
      const cmtk::UniformVolume* uniform = cmtk::UniformVolume::SmartPtr::DynamicCastFrom( volume );
      const cmtk::TypedArray* dataArray = volume->GetData();
      
      if ( MachineReadable )
	{
	fprintf( stdout, "FNAME\t%s\n", next );            
	fprintf( stdout, "FORMAT\t%s\n", volume->m_MetaInformation[CMTK_META_FILEFORMAT_ORIGINAL].c_str() );
	fprintf( stdout, "XDIM\t%d\nYDIM\t%d\nZDIM\t%d\n", volume->GetDims( cmtk::AXIS_X ), volume->GetDims( cmtk::AXIS_Y ), volume->GetDims( cmtk::AXIS_Z ) );
	fprintf( stdout, "ORIENT\t%s\n", orientOriginal ? orientOriginal : "UNKNOWN" );
      
	if ( uniform ) 
	  {
	  fprintf( stdout, "GRID\tUniform\nXPIX\t%f\nYPIX\t%f\nZPIX\t%f\nXFOV\t%f\nYFOV\t%f\nZFOV\t%f\n",
		   uniform->m_Delta[0], uniform->m_Delta[1], uniform->m_Delta[2],
		   uniform->Size[0], uniform->Size[1], uniform->Size[2] );
	  }
	
	fprintf( stdout, "XORIGIN\t%f\nYORIGIN\t%f\nZORIGIN\t%f\n", volume->m_Origin.XYZ[0], volume->m_Origin.XYZ[1], volume->m_Origin.XYZ[2] );
	if ( volume->MetaKeyExists(CMTK_META_SPACE_UNITS_STRING ) )
	  fprintf( stdout, "UNITS\t%s\n", volume->m_MetaInformation[CMTK_META_SPACE_UNITS_STRING].c_str() );
	
	if ( dataArray ) 
	  {
	  cmtk::Types::DataItem min = 0, max = 0;
	  if ( dataArray->GetRange( min, max ) ) 
	    {
	    fprintf( stdout, "DTYPE\t%s\nMINDATA\t%f\nMAXDATA\t%f\n", cmtk::DataTypeName[ dataArray->GetType() ], static_cast<float>( min ), static_cast<float>( max ) );
	    } 
	  else
	    {
	    fprintf( stdout, "DTYPE\tNONE\n" );
	    }
	  }
	}
      else
	{
	fprintf( stdout, "File: %s\n", next );            
	fprintf( stdout, "File format: %s\n", volume->m_MetaInformation[CMTK_META_FILEFORMAT_ORIGINAL].c_str() );
	fprintf( stdout, "%d x %d x %d voxels\n", volume->GetDims( cmtk::AXIS_X ), volume->GetDims( cmtk::AXIS_Y ), volume->GetDims( cmtk::AXIS_Z ) );

	fprintf( stdout, "Original image orientation: %s\n", orientOriginal ? orientOriginal : "UNKNOWN" );
      
	const char* spaceUnits = "";
	if ( volume->MetaKeyExists(CMTK_META_SPACE_UNITS_STRING ) )
	  spaceUnits = volume->m_MetaInformation[CMTK_META_SPACE_UNITS_STRING].c_str();
	
	if ( uniform ) 
	  {
	  fprintf( stdout, "Uniform volume\n%f x %f x %f [%s] voxel size\n%f x %f x %f [%s] volume size\n",
		   uniform->m_Delta[0], uniform->m_Delta[1], uniform->m_Delta[2], spaceUnits,
		   uniform->Size[0], uniform->Size[1], uniform->Size[2], spaceUnits );
	  }

	fprintf( stdout, "Volume origin (%f,%f,%f)\n", volume->m_Origin.XYZ[0], volume->m_Origin.XYZ[1], volume->m_Origin.XYZ[2] );
      
	if ( dataArray ) 
	  {
	  cmtk::Types::DataItem min = 0, max = 0;
	  if ( dataArray->GetRange( min, max ) ) 
	    {
	    cmtk::StdErr.printf( "Data type %s, range [%f .. %f]\n", cmtk::DataTypeName[ dataArray->GetType() ], static_cast<float>( min ), static_cast<float>( max ) );
	    } 
	  else
	    {
	    cmtk::StdErr << "Image does not contain valid data.\n";
	    }
	  }
	}
     
      size_t index = 0;
      std::list<cmtk::Vector3D>::const_iterator it = ProbeListIndex.begin();      
      for ( ; it != ProbeListIndex.end(); ++it, ++index )
	{
	cmtk::Types::DataItem data;
	if ( volume->GetDataAt( data, (int)it->XYZ[0], (int)it->XYZ[1], (int)it->XYZ[2] ) )
	  {
	  fprintf( stdout, "Probe %02d = %f (%e)\n", (int)index, data, data );
	  }
	else
	  {
	  fprintf( stdout, "Probe %02d = NAN\n", (int)index );
	  }
	}
     
      if ( PrintCOM || PrintFOM )
	{
	cmtk::Vector3D com, fom;
	if ( PrintFOM )
	  com = volume->GetCenterOfMass( fom );
	else
	  com = volume->GetCenterOfMass();
       
	if ( PrintCOM )
	  {
	  fprintf( stdout, "COM\t%f %f %f\n", com.XYZ[0], com.XYZ[1], com.XYZ[2] );
	  }
	if ( PrintFOM )
	  {
	  fprintf( stdout, "FOM\t%f %f %f\n", fom.XYZ[0], fom.XYZ[1], fom.XYZ[2] );
	  }
	}
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

#ifdef CMTK_BUILD_MPI
  MPI::Finalize();
#endif

  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace describe
} // namespace apps
} // namespace cmtk
#endif


