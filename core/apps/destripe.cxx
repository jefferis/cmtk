/*
//
//  Copyright 2011 SRI International
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
//  $Revision: 3001 $
//
//  $LastChangedDate: 2011-03-16 12:14:22 -0700 (Wed, 16 Mar 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkDebugOutput.h>

#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkDataGrid.h>
#include <Base/cmtkMathUtil.h>
#include <Base/cmtkGaussianKernel.h>

#include <vector>

int SliceAxis = -1;

double KernelFWHM = 1;

const char* InputFilePath = NULL;
const char* OutputFilePath = NULL;

int
doMain ( const int argc, const char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Destripe volume image data." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This program corrects stripe artifacts in acquired image stacks which can result from between-slice intensity scale differences." );
    
    typedef cmtk::CommandLine::Key Key;    

    cl.BeginGroup( "SliceOrient", "Slice Orientation" );
    cmtk::CommandLine::EnumGroup<int>::SmartPtr sliceGroup = cl.AddEnum( "slice-axis", &SliceAxis, "Define slice direction axis: this is the through-slice direction of the acquisition." );
    sliceGroup->AddSwitch( Key( "guess-from-input" ), -1, "Guess from input image" );
    sliceGroup->AddSwitch( Key( 'a', "axial" ), (int)cmtk::AXIS_Z, "Sliced axial images" );
    sliceGroup->AddSwitch( Key( 's', "sagittal" ),(int)cmtk::AXIS_X, "Sliced sagittal images" );
    sliceGroup->AddSwitch( Key( 'c', "coronal" ), (int)cmtk::AXIS_Y, "Sliced coronal images" );
    sliceGroup->AddSwitch( Key( 'x', "slice-x" ), (int)cmtk::AXIS_X, "Sliced along x axis" );
    sliceGroup->AddSwitch( Key( 'y', "slice-y" ), (int)cmtk::AXIS_Y, "Sliced along y axis" );
    sliceGroup->AddSwitch( Key( 'z', "slice-z" ), (int)cmtk::AXIS_Z, "Sliced along z axis" );

    cl.AddOption( Key( "kernel-fwhm" ), &KernelFWHM, "Gaussian kernel full width at half maximum (FWHM)." );
    
    cl.AddParameter( &InputFilePath, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &OutputFilePath, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT | cmtk::CommandLine::PROPS_OPTIONAL );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return 1;
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InputFilePath ) );
  if ( ! volume )
    {
    cmtk::StdErr << "ERROR: could not read volume " << InputFilePath << "\n";
    return 1;
    }

  const cmtk::DataGrid::IndexType volumeDims = volume->GetDims();

  // guess slice orientation - if two dimensions are equal, the thir is usually the slice direction
  if ( SliceAxis == -1 )
    {
    if ( volumeDims[0] == volumeDims[1] )
      SliceAxis = 2;
    else if ( volumeDims[0] == volumeDims[2] )
      SliceAxis = 1;
    else if ( volumeDims[1] == volumeDims[2] )
      SliceAxis = 0;
    else
      {
      cmtk::StdErr << "ERROR: cannot guess slice axis when all three image dimensions are different.\n";
      return 1;
      }
    }

  const int nSlices = volumeDims[SliceAxis];
  std::vector<cmtk::Types::DataItem> sliceProjection( nSlices );
  
  const int idxX = (SliceAxis==1 || SliceAxis==2) ? 0 : 1;
  const int idxY = (SliceAxis==0 || SliceAxis==2) ? 1 : 2;

#pragma omp parallel for
  for ( int slice = 0; slice < nSlices; ++slice )
    {
    cmtk::DataGrid::IndexType idx;
    idx[SliceAxis] = slice;

    cmtk::Types::DataItem sum = 0;
    size_t count = 0;
    cmtk::Types::DataItem value;
    for ( idx[idxX] = 0; idx[idxX] < volumeDims[idxX]; ++idx[idxX] )
      {
      for ( idx[idxY] = 0; idx[idxY] < volumeDims[idxY]; ++idx[idxY] )
	{
	if ( volume->GetDataAt( value, volume->GetOffsetFromIndex( idx ) ) )
	  {
	  sum += value;
	  ++count;
	  }
	}
      }

    if ( count )
      sliceProjection[slice] = sum / count;
    else
      sliceProjection[slice] = 0;      
    }

  const std::vector<cmtk::Types::DataItem> kernel = cmtk::GaussianKernel<cmtk::Types::DataItem>::GetHalfKernel( cmtk::Units::GaussianFWHM( KernelFWHM ) );

  std::vector<cmtk::Types::DataItem> smoothed( nSlices, 0.0 );
  for ( int slice = 0; slice < nSlices; ++ slice )
    {
    cmtk::Types::DataItem kernelSum = kernel[0];
    smoothed[slice] = kernel[0] * sliceProjection[slice];
    for ( int ofs = 1; ofs < static_cast<int>( kernel.size() ); ++ofs )
      {
      if ( slice - ofs  >= 0 )
	{
	kernelSum += kernel[ofs];
	smoothed[slice] += kernel[ofs] * sliceProjection[slice-ofs];
	}
      if ( slice + ofs < nSlices )
	{
	kernelSum += kernel[ofs];
	smoothed[slice] += kernel[ofs] * sliceProjection[slice+ofs];
	}
      }
    smoothed[slice] /= kernelSum;
    }

#pragma omp parallel for
  for ( int slice = 0; slice < nSlices; ++slice )
    {
    cmtk::DataGrid::IndexType idx;
    idx[SliceAxis] = slice;

    const cmtk::Types::DataItem correction = smoothed[slice] / sliceProjection[slice];

    cmtk::Types::DataItem value;
    for ( idx[idxX] = 0; idx[idxX] < volumeDims[idxX]; ++idx[idxX] )
      {
      for ( idx[idxY] = 0; idx[idxY] < volumeDims[idxY]; ++idx[idxY] )
	{
	if ( volume->GetDataAt( value, volume->GetOffsetFromIndex( idx ) ) )
	  {
	  volume->SetDataAt( value * correction, volume->GetOffsetFromIndex( idx ) );
	  }
	}
      }
    }

  for ( int slice = 0; slice < nSlices; ++slice )
    {
    cmtk::StdOut << sliceProjection[slice] << " ";
    }
  cmtk::StdOut << "\n\n";

  for ( int slice = 0; slice < nSlices; ++slice )
    {
    cmtk::StdOut << smoothed[slice] << " ";
    }
  cmtk::StdOut << "\n\n";

  cmtk::VolumeIO::Write( *volume, OutputFilePath );
  
  return 0;
}

#include "cmtkSafeMain"
