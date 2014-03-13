/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkProgress.h>
#include <System/cmtkThreads.h>

#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkDeformationField.h>
#include <Base/cmtkXformList.h>

#include <IO/cmtkXformIO.h>
#include <IO/cmtkXformListIO.h>
#include <IO/cmtkVolumeIO.h>

#include <stdio.h>

#include <vector>
#include <string>
#include <numeric>

#ifdef CMTK_USE_GCD
#  include <dispatch/dispatch.h>
#endif

bool Mask = false;
bool OutputAbsolute = false;

std::string RefFileName;
std::string OutFileName;

std::string TargetImagePath;
bool TargetImageFSL = false;

std::vector<std::string> InputXformPaths;

std::string Downsample;

cmtk::Types::Coordinate InversionToleranceFactor = 0.1;

int
doMain ( const int argc, const char *argv[] ) 
{
  cmtk::Threads::GetNumberOfThreads();

  cmtk::UniformVolume::SmartPtr volume;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Transformation to Deformation Field" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Convert parametric rigid or nonrigid transformation to deformation field, sampled at pixel locations of a given reference image" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'm', "mask" ), &Mask, true, "Use reference image pixels as a binary mask." );

    cl.AddOption( Key( "target-image-fsl" ), &TargetImagePath, "Path to optional target image to be used in FSL. This is the image that the transformation (and thus the created deformation field) maps to. "
		  "When creating a deformation field to be used with FSL's \"applywarp\" tool, is may be necessary to determine whether the x coordinate of the resulting deformation vectors has to be flipped.", &TargetImageFSL );
    cl.AddOption( Key( "inversion-tolerance-factor" ), &InversionToleranceFactor, "Factor for numerical tolerance of B-spline inversion [multiples of minimum grid pixel size; default=0.1]" );
    cl.AddOption( Key( "downsample" ), &Downsample, "Downsample grid by factors 'x,y,z' or by single factor 'xyz'" );

    cl.AddSwitch( Key( "output-absolute" ), &OutputAbsolute, true, "Make output deformation field with absolute target location vectors, rather than relative offset vectors." );

    cl.AddParameter( &OutFileName, "OutputPath", "Path for the output deformation field." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddParameter( &RefFileName, "ReferenceImage", "Input reference grid path. The dimensions and pixel size of this image determine the geometry of the output." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameterVector( &InputXformPaths, "XformList", "List of concatenated transformations. Insert '--inverse' to use the inverse of the transformation listed next. "
			   "(If the first transformation in the sequence is inverted, then '--inverse' must be preceded by '--', i.e., use '-- --inverse xform.path')." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );  
 
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  if ( Mask )
    volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( RefFileName ) );
  else
    volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( RefFileName ) );
  if ( ! volume ) 
    {
    cmtk::StdErr << "Could not read reference volume " << RefFileName << "\n";
    throw cmtk::ExitException(1);
    }          

  // Do we have a target image given (presumably for FSL use at this point)?
  cmtk::AffineXform::MatrixType targetImageMatrix = cmtk::AffineXform::MatrixType::Identity();
  if ( !TargetImagePath.empty() )
    {
    cmtk::UniformVolume::SmartPtr targetImageGrid( cmtk::VolumeIO::ReadGrid( TargetImagePath ) );
    if ( ! targetImageGrid ) 
      {
      cmtk::StdErr << "Could not read target volume grid from " << TargetImagePath << "\n";
      throw cmtk::ExitException(1);
      }

    // Get initial matrix from index to physical space
    targetImageMatrix = targetImageGrid->GetImageToPhysicalMatrix().GetInverse();

    if ( TargetImageFSL )
      {
      // Warn if not "absolute" output.
      if ( ! OutputAbsolute )
	{
	cmtk::StdErr << "WARNING: generating deformation field for FSL, but \"--output--absolute\" is not active.\n";
	}
      
      // Does the matrix have a positive determinant? Then this is what FSL calls "neurological orientation" and we need to flip x
      if ( targetImageGrid->GetImageToPhysicalMatrix().GetTopLeft3x3().Determinant() > 0 )
	{
	cmtk::AffineXform::MatrixType targetImageMatrixFlip = cmtk::AffineXform::MatrixType::Identity();
	targetImageMatrixFlip[0][0] = -1;
	targetImageMatrixFlip[3][0] = (targetImageGrid->m_Dims[0]-1) * targetImageGrid->m_Delta[0];
	targetImageMatrix = targetImageMatrix * targetImageMatrixFlip;
	}
      }

    targetImageGrid = targetImageGrid->GetReoriented( "RAS" );
    targetImageMatrix = targetImageGrid->GetImageToPhysicalMatrix() * targetImageMatrix;
    }    

  cmtk::XformList xformList = cmtk::XformListIO::MakeFromStringList( InputXformPaths );  
  xformList.SetEpsilon( InversionToleranceFactor * volume->GetMinDelta() );
  
  const cmtk::XformList& xformListRef = xformList; // need this to work around GCD bug

  if ( !Downsample.empty() )
    {
    int factors[3] = { 1, 1, 1 };
    const size_t nFactors = sscanf( Downsample.c_str(), "%6d,%6d,%6d", factors, factors+1, factors+2 );
    if ( nFactors == 1 )
      {
      factors[1] = factors[2] = factors[0];
      }
    else
      {
      if ( nFactors != 3 )
	{
	cmtk::StdErr << "ERROR: downsampling factors must either be three integers, x,y,z, or a single integer\n";
	throw cmtk::ExitException( 1 );
	}
      }
    volume = cmtk::UniformVolume::SmartPtr( volume->GetDownsampledAndAveraged( factors ) );
    }
  
  cmtk::DeformationField::SmartPtr dfield( new cmtk::DeformationField( volume ) );
  
  const cmtk::DataGrid::IndexType& dims = volume->GetDims();
  cmtk::Progress::Begin( 0, dims[cmtk::AXIS_Z], 1, "Deformation field generation" );

#ifdef CMTK_USE_GCD
  dispatch_apply( dims[2], dispatch_get_global_queue(0, 0), ^(size_t z){
#else
#pragma omp parallel for
  for ( int z = 0; z < dims[cmtk::AXIS_Z]; ++z )
#endif
    {
    cmtk::Xform::SpaceVectorType v0, v1;

    cmtk::Progress::SetProgress( z );

    size_t offset = 3 * z * dims[cmtk::AXIS_X] * dims[cmtk::AXIS_Y];
    for ( int y = 0; y < dims[cmtk::AXIS_Y]; ++y )
      {
      for ( int x = 0; x < dims[cmtk::AXIS_X]; ++x, offset+=3 ) 
	{
	v1 = v0 = volume->GetGridLocation( x, y, z );

	bool invalid = true;
	if ( (!Mask) || (volume->GetDataAt( x, y, z ) > 0) )
	  {
	  invalid = !xformListRef.ApplyInPlace( v1 );
	  }

	if ( !invalid )
	  {
	  // do any transformations to account for target image coordinate changes due to reorientation
	  if ( TargetImageFSL )
	    {
	    v1 *= targetImageMatrix;
	    }

	  if ( !OutputAbsolute )
	    v1 -= v0;
	  }
	else
	  {
	  v1 = cmtk::Vector3D( std::numeric_limits<cmtk::Types::Coordinate>::quiet_NaN() );
	  }

	// store final vector.
	if ( TargetImageFSL )
	  {
	  const size_t flipOffset = offset + 3 * ((dims[cmtk::AXIS_X] - 1) - 2 * x);
	  dfield->m_Parameters[flipOffset+0] = v1[0];
	  dfield->m_Parameters[flipOffset+1] = v1[1];
	  dfield->m_Parameters[flipOffset+2] = v1[2];
	  }
	else
	  {
	  dfield->m_Parameters[offset+0] = v1[0];
	  dfield->m_Parameters[offset+1] = v1[1];
	  dfield->m_Parameters[offset+2] = v1[2];
	  }
	}
      }
    }
#ifdef CMTK_USE_GCD
		  });
#endif

  cmtk::Progress::Done();

  cmtk::XformIO::Write( dfield, OutFileName );

  return 0;
}

#include "cmtkSafeMain"
