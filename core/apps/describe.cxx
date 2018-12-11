/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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
#include <Base/cmtkVector3D.h>
#include <Base/cmtkMetaInformationObject.h>

#include <Base/cmtkXform.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkDeformationField.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>

#include <stdio.h>

const char* ReadOrientation = NULL;

bool MachineReadable = false;

void
DescribeImage( const char* path )
{
  cmtk::UniformVolume::SmartPtr volume;
  try
    {
    if ( ReadOrientation )
      volume = cmtk::VolumeIO::ReadOriented( path, ReadOrientation );
    else
      volume = cmtk::VolumeIO::Read( path );
    }
  catch (...) 
    {
    if ( ReadOrientation )
      volume = cmtk::VolumeIO::ReadGridOriented( path, ReadOrientation );
    else
      volume = cmtk::VolumeIO::ReadGrid( path );
    }
  
  const char* orientOriginal = volume->GetMetaInfo( cmtk::META_IMAGE_ORIENTATION_ORIGINAL ).c_str();
  const cmtk::TypedArray* dataArray = volume->GetData();
  
  if ( MachineReadable )
    {
    fprintf( stdout, "FNAME\t%s\n", path );            
    fprintf( stdout, "FORMAT\t%s\n", volume->GetMetaInfo( cmtk::META_FILEFORMAT_ORIGINAL ).c_str() );
    
    if ( volume->MetaKeyExists( cmtk::META_IMAGE_DESCRIPTION ) )
      fprintf( stdout, "DESCRIP\t\"%s\"\n", volume->GetMetaInfo( cmtk::META_IMAGE_DESCRIPTION ).c_str() );
    
    fprintf( stdout, "XDIM\t%d\nYDIM\t%d\nZDIM\t%d\n", volume->GetDims()[cmtk::AXIS_X], volume->GetDims()[cmtk::AXIS_Y], volume->GetDims()[cmtk::AXIS_Z] );
    fprintf( stdout, "ORIENT\t%s\n", orientOriginal ? orientOriginal : "UNKNOWN" );
    
    fprintf( stdout, "GRID\tUniform\nXPIX\t%f\nYPIX\t%f\nZPIX\t%f\nXFOV\t%f\nYFOV\t%f\nZFOV\t%f\n",
	     volume->m_Delta[0], volume->m_Delta[1], volume->m_Delta[2],
	     volume->m_Size[0], volume->m_Size[1], volume->m_Size[2] );
    
    fprintf( stdout, "XORIGIN\t%f\nYORIGIN\t%f\nZORIGIN\t%f\n", volume->m_Offset[0], volume->m_Offset[1], volume->m_Offset[2] );
    if ( volume->MetaKeyExists(cmtk::META_SPACE_UNITS_STRING ) )
      fprintf( stdout, "UNITS\t%s\n", volume->GetMetaInfo( cmtk::META_SPACE_UNITS_STRING ).c_str() );
    
    const cmtk::AffineXform::MatrixType a2p = volume->GetImageToPhysicalMatrix();
    fprintf( stdout, "I2PMAT0\t%f\t%f\t%f\t%f\nI2PMAT1\t%f\t%f\t%f\t%f\nI2PMAT2\t%f\t%f\t%f\t%f\nI2PMAT3\t%f\t%f\t%f\t%f\n", 
	     a2p[0][0], a2p[1][0], a2p[2][0], a2p[3][0], 
	     a2p[0][1], a2p[1][1], a2p[2][1], a2p[3][1], 
	     a2p[0][2], a2p[1][2], a2p[2][2], a2p[3][2],
	     a2p[0][3], a2p[1][3], a2p[2][3], a2p[3][3] );
    
    if ( dataArray ) 
      {
      fprintf( stdout, "DTYPE\t%s\n", cmtk::DataTypeName[ dataArray->GetType() ] );
      if ( dataArray->GetDataSize() )
	{
	const cmtk::Types::DataItemRange range = dataArray->GetRange();
	fprintf( stdout, "MINDATA\t%f\nMAXDATA\t%f\n", static_cast<float>( range.m_LowerBound ), static_cast<float>( range.m_UpperBound ) );
	}
      }
    }
  else
    {
    fprintf( stdout, "File: %s\n", path );            
    fprintf( stdout, "File format: %s\n", volume->GetMetaInfo( cmtk::META_FILEFORMAT_ORIGINAL ).c_str() );
    
    if ( volume->MetaKeyExists( cmtk::META_IMAGE_DESCRIPTION ) )
      fprintf( stdout, "Description: \"%s\"\n", volume->GetMetaInfo( cmtk::META_IMAGE_DESCRIPTION ).c_str() );
    
    fprintf( stdout, "%d x %d x %d voxels\n", volume->GetDims()[cmtk::AXIS_X], volume->GetDims()[cmtk::AXIS_Y], volume->GetDims()[cmtk::AXIS_Z] );
    
    fprintf( stdout, "Original image orientation: %s\n", orientOriginal ? orientOriginal : "UNKNOWN" );
    
    const char* spaceUnits = "";
    if ( volume->MetaKeyExists(cmtk::META_SPACE_UNITS_STRING ) )
      spaceUnits = volume->GetMetaInfo( cmtk::META_SPACE_UNITS_STRING ).c_str();
    
    fprintf( stdout, "Uniform volume\n%f x %f x %f [%s] voxel size\n%f x %f x %f [%s] volume size\n",
	     volume->m_Delta[0], volume->m_Delta[1], volume->m_Delta[2], spaceUnits,
	     volume->m_Size[0], volume->m_Size[1], volume->m_Size[2], spaceUnits );
    
    fprintf( stdout, "Volume origin (%f,%f,%f)\n", volume->m_Offset[0], volume->m_Offset[1], volume->m_Offset[2] );
    
    const cmtk::AffineXform::MatrixType a2p = volume->GetImageToPhysicalMatrix();
    fprintf( stdout, "\nImage-to-physical matrix:\n  %f %f %f %f \n %f %f %f %f \n %f %f %f %f \n %f %f %f %f \n\n", 
	     a2p[0][0], a2p[1][0], a2p[2][0], a2p[3][0], 
	     a2p[0][1], a2p[1][1], a2p[2][1], a2p[3][1], 
	     a2p[0][2], a2p[1][2], a2p[2][2], a2p[3][2],
	     a2p[0][3], a2p[1][3], a2p[2][3], a2p[3][3] );
    
    if ( dataArray ) 
      {
      cmtk::StdOut.printf( "Data type %s", cmtk::DataTypeName[ dataArray->GetType() ] );
      if ( dataArray->GetDataSize() )
	{
	const cmtk::Types::DataItemRange range = dataArray->GetRange();
	cmtk::StdOut.printf( ", range [%f .. %f]\n", static_cast<float>( range.m_LowerBound ), static_cast<float>( range.m_UpperBound ) );
	}
      else
	{
	cmtk::StdOut << "\n";
	}
      } 
    else
      {
      cmtk::StdOut << "Image does not contain valid data.\n";
      }
    }
}

void
DescribeWarpXform( const cmtk::WarpXform& warpXform )
{
  if ( MachineReadable )
    {
    cmtk::StdOut << "DOMAINX\t" << warpXform.m_Domain[0] << "\n";
    cmtk::StdOut << "DOMAINY\t" << warpXform.m_Domain[1] << "\n";
    cmtk::StdOut << "DOMAINZ\t" << warpXform.m_Domain[2] << "\n";
    cmtk::StdOut << "GRIDX\t" << warpXform.m_Dims[0] << "\n";
    cmtk::StdOut << "GRIDY\t" << warpXform.m_Dims[1] << "\n";
    cmtk::StdOut << "GRIDZ\t" << warpXform.m_Dims[2] << "\n";
    cmtk::StdOut << "DELTAX\t" << warpXform.m_Spacing[0] << "\n";
    cmtk::StdOut << "DELTAY\t" << warpXform.m_Spacing[1] << "\n";
    cmtk::StdOut << "DELTAZ\t" << warpXform.m_Spacing[2] << "\n";
    cmtk::StdOut << "OFFSETX\t" << warpXform.m_Offset[0] << "\n";
    cmtk::StdOut << "OFFSETY\t" << warpXform.m_Offset[1] << "\n";
    cmtk::StdOut << "OFFSETZ\t" << warpXform.m_Offset[2] << "\n";

    cmtk::AffineXform::SmartConstPtr initialAffine = warpXform.GetInitialAffineXform();
    if ( initialAffine )
      {
      cmtk::StdOut << "BULK_XLATEX\t" << initialAffine->RetXlate()[0] << "\n";
      cmtk::StdOut << "BULK_XLATEY\t" << initialAffine->RetXlate()[1] << "\n";
      cmtk::StdOut << "BULK_XLATEZ\t" << initialAffine->RetXlate()[2] << "\n";
      cmtk::StdOut << "BULK_ANGLEX\t" << initialAffine->RetAngles()[0] << "\n";
      cmtk::StdOut << "BULK_ANGLEY\t" << initialAffine->RetAngles()[1] << "\n";
      cmtk::StdOut << "BULK_ANGLEZ\t" << initialAffine->RetAngles()[2] << "\n";
      cmtk::StdOut << "BULK_SCALEX\t" << initialAffine->RetScales()[0] << "\n";
      cmtk::StdOut << "BULK_SCALEY\t" << initialAffine->RetScales()[1] << "\n";
      cmtk::StdOut << "BULK_SCALEZ\t" << initialAffine->RetScales()[2] << "\n";
      cmtk::StdOut << "BULK_SHEARX\t" << initialAffine->RetShears()[0] << "\n";
      cmtk::StdOut << "BULK_SHEARY\t" << initialAffine->RetShears()[1] << "\n";
      cmtk::StdOut << "BULK_SHEARZ\t" << initialAffine->RetShears()[2] << "\n";
      }
    }
  else
    {
    cmtk::StdOut << "Domain: " << warpXform.m_Domain[0] << " x " << warpXform.m_Domain[1] << " x " << warpXform.m_Domain[2] << "\n"; 
    cmtk::StdOut << "Control point grid: " << warpXform.m_Dims[0] << " x " << warpXform.m_Dims[1] << " x " << warpXform.m_Dims[2] << "\n"; 
    cmtk::StdOut << "Control point spacing: " << warpXform.m_Spacing[0] << " x " << warpXform.m_Spacing[1] << " x " << warpXform.m_Spacing[2] << "\n"; 
    cmtk::StdOut << "First control point offset: " << warpXform.m_Offset[0] << " x " << warpXform.m_Offset[1] << " x " << warpXform.m_Offset[2] << "\n"; 

    cmtk::AffineXform::SmartConstPtr initialAffine = warpXform.GetInitialAffineXform();
    if ( initialAffine )
      {
      cmtk::StdOut << "Bulk affine transformation:" << "\n";
      cmtk::StdOut << "\tTranslation: " << initialAffine->RetXlate()[0] << ", " << initialAffine->RetXlate()[1] << ", " << initialAffine->RetXlate()[2] << "\n";
      cmtk::StdOut << "\tRotation angles: " << initialAffine->RetAngles()[0] << ", " << initialAffine->RetAngles()[1] << ", " << initialAffine->RetAngles()[2] << "\n";
      cmtk::StdOut << "\tScale factors: " << initialAffine->RetScales()[0] << ", " << initialAffine->RetScales()[1] << ", " << initialAffine->RetScales()[2] << "\n";
      cmtk::StdOut << "\tShear coefficients: " << initialAffine->RetShears()[0] << ", " << initialAffine->RetShears()[1] << ", " << initialAffine->RetShears()[2] << "\n";
      }
    }
}

void
DescribeXform( const char* path )
{
  cmtk::Xform::SmartConstPtr xform = cmtk::XformIO::Read( path );

  if ( ! xform )
    {
    cmtk::StdErr << "Could not read transformation from " << path << "\n";
    }

  const std::string fPath = xform->GetMetaInfo( cmtk::META_XFORM_FIXED_IMAGE_PATH, "" );
  const std::string mPath = xform->GetMetaInfo( cmtk::META_XFORM_MOVING_IMAGE_PATH, "" );

  if ( MachineReadable )
    {
    cmtk::StdOut << "PATH\t" << path << "\n";
    cmtk::StdOut << "NPARAMS\t" << xform->ParamVectorDim() << "\n";
    if ( fPath != "" )
      cmtk::StdOut << "FXIMAGE\t" << fPath << "\n";
    if ( mPath != "" )
      cmtk::StdOut << "MVIMAGE\t" << mPath << "\n";
    }
  else
    {
    cmtk::StdOut << "File: " << path << "\n";
    cmtk::StdOut << "Number of parameters: " << xform->ParamVectorDim() << "\n";
    if ( fPath != "" )
      cmtk::StdOut << "Fixed image path: " << fPath << "\n";
    if ( mPath != "" )
      cmtk::StdOut << "Moving image path: " << mPath << "\n";
    }

  cmtk::SplineWarpXform::SmartConstPtr splineWarpXform = cmtk::SplineWarpXform::SmartConstPtr::DynamicCastFrom( xform );
  if ( splineWarpXform )
    {    
    if ( MachineReadable )
      {
      cmtk::StdOut << "MODEL\tBSPLINE\n";
      }
    else
      {
      cmtk::StdOut << "Transformation model: B-spline free-form deformation\n";
      }
    DescribeWarpXform( *splineWarpXform );
    }

  cmtk::DeformationField::SmartConstPtr dfield = cmtk::DeformationField::SmartConstPtr::DynamicCastFrom( xform );
  if ( dfield )
    {    
    if ( MachineReadable )
      {
      cmtk::StdOut << "MODEL\tDFIELD\n";
      }
    else
      {
      cmtk::StdOut << "Transformation model: Deformation field\n";
      }
    DescribeWarpXform( *dfield );
    }

  cmtk::AffineXform::SmartConstPtr affine = cmtk::AffineXform::SmartConstPtr::DynamicCastFrom( xform );
  if ( affine )
    {    
    if ( MachineReadable )
      {
      cmtk::StdOut << "MODEL\tAFFINE\n";
      cmtk::StdOut << "XLATEX\t" << affine->RetXlate()[0] << "\n";
      cmtk::StdOut << "XLATEY\t" << affine->RetXlate()[1] << "\n";
      cmtk::StdOut << "XLATEZ\t" << affine->RetXlate()[2] << "\n";
      cmtk::StdOut << "ANGLEX\t" << affine->RetAngles()[0] << "\n";
      cmtk::StdOut << "ANGLEY\t" << affine->RetAngles()[1] << "\n";
      cmtk::StdOut << "ANGLEZ\t" << affine->RetAngles()[2] << "\n";
      cmtk::StdOut << "SCALEX\t" << affine->RetScales()[0] << "\n";
      cmtk::StdOut << "SCALEY\t" << affine->RetScales()[1] << "\n";
      cmtk::StdOut << "SCALEZ\t" << affine->RetScales()[2] << "\n";
      cmtk::StdOut << "SHEARX\t" << affine->RetShears()[0] << "\n";
      cmtk::StdOut << "SHEARY\t" << affine->RetShears()[1] << "\n";
      cmtk::StdOut << "SHEARZ\t" << affine->RetShears()[2] << "\n";
      }
    else
      {
      cmtk::StdOut << "Transformation model: Affine\n";
      cmtk::StdOut << "Translation: " << affine->RetXlate()[0] << ", " << affine->RetXlate()[1] << ", " << affine->RetXlate()[2] << "\n";
      cmtk::StdOut << "Rotation angles: " << affine->RetAngles()[0] << ", " << affine->RetAngles()[1] << ", " << affine->RetAngles()[2] << "\n";
      cmtk::StdOut << "Scale factors: " << affine->RetScales()[0] << ", " << affine->RetScales()[1] << ", " << affine->RetScales()[2] << "\n";
      cmtk::StdOut << "Shear coefficients: " << affine->RetShears()[0] << ", " << affine->RetShears()[1] << ", " << affine->RetShears()[2] << "\n";
      }
    }

  cmtk::StdOut << "\n";
}

int
doMain( int argc, const char *argv[] )
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Describe image and transformation file formats and parameters" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool prints a detailed description of the input files as either image(s) or transformation(s)." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "describe [options] file0 [file1 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'm', "machine-readable" ), &MachineReadable, true, "Print output in format that is easy to parse automatically." );
    cl.AddSwitch( Key( "read-ras" ), &ReadOrientation, "RAS", "Read all images in RAS orientation" );

    cl.Parse( argc, const_cast<const char**>( argv ) );

    for ( const char* next = cl.GetNext(); next; next = cl.GetNextOptional() ) 
      {
      const cmtk::FileFormatID id = cmtk::FileFormat::Identify( next );
      switch ( id )
	{
	case cmtk::FILEFORMAT_NEXIST:
	  cmtk::StdErr << "File does not exist: " << next << "\n";
	  break;
	case cmtk::FILEFORMAT_UNKNOWN:
	  cmtk::StdErr << "Unknown file format: " << next << "\n";
	  break;
	default:
	  if ( cmtk::FileFormat::IsImage( id ) )
	    {
	    DescribeImage( next );
	    }
	  else if ( cmtk::FileFormat::IsXform( id ) )
	    {
	    DescribeXform( next );
	    }
	  else
	    {
	    cmtk::StdErr << "File is neither image nor transformation: " << next << "\n";
	    }
	}
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  return 0;
}

#include "cmtkSafeMain"

