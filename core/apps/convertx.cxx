/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2010 SRI International
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
#include <cmtkTimers.h>

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>

#include <cmtkStudyList.h>
#include <cmtkClassStreamStudyList.h>

#include <math.h>
#include <cmtkMathFunctionWrappers.h>

#include <cmtkImageOperation.h>
#include <cmtkImageOperationFlip.h>
#include <cmtkImageOperationApplyMask.h>
#include <cmtkImageOperationErodeDilate.h>
#include <cmtkImageOperationBoundaryMap.h>
#include <cmtkImageOperationDownsample.h>
#include <cmtkImageOperationHistogramEqualization.h>
#include <cmtkImageOperationCropRegion.h>
#include <cmtkImageOperationCropThreshold.h>
#include <cmtkImageOperationScaleToRange.h>
#include <cmtkImageOperationThreshold.h>
#include <cmtkImageOperationMedianFilter.h>
#include <cmtkImageOperationMedialSkeleton.h>
#include <cmtkImageOperationGaussFilter.h>
#include <cmtkImageOperationDistanceMap.h>

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <set>
#include <map>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace convert
{
#endif

bool Verbose = false;

int
main( int argc, char* argv[] )
{
  cmtk::Types::DataItem paddingDataValue = 0.0;
  bool paddingDataFlag = false;

  const char* imagePathIn = NULL;
  const char* imagePathOut = NULL;

  try 
    {
    cmtk::CommandLine cl( argc, argv );  
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Convert between image file formats and data types. Also apply simple, general-purpose image operations in the process. "
		       "An arbitrary number of operations can be specified on the command line, which will be applied exactly in the order given." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] infile outfile" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.BeginGroup( "Input", "Input Image Controls" );
    cl.AddOption( Key( 'N', "set-padding" ), &paddingDataValue, "Set padding data for input image. All pixels in the input image that have this value will be ignored in all operations.", 
		  &paddingDataFlag );
    cl.EndGroup();

    cl.BeginGroup( "Conversion", "Data Type Conversion" );
    cl.AddCallback( Key( "char" ), &cmtk::ImageOperationConvertType::NewChar, "8 bits, signed integer" );
    cl.AddCallback( Key( "byte" ), &cmtk::ImageOperationConvertType::NewByte, "8 bits, unsigned integer" );
    cl.AddCallback( Key( "short" ), &cmtk::ImageOperationConvertType::NewShort, "16 bits, signed integer" );
    cl.AddCallback( Key( "ushort" ), &cmtk::ImageOperationConvertType::NewUShort, "16 bits, unsigned integer" );
    cl.AddCallback( Key( "int" ), &cmtk::ImageOperationConvertType::NewInt, "32 bits signed integer" );
    cl.AddCallback( Key( "float" ), &cmtk::ImageOperationConvertType::NewFloat, "32 bits floating point" );
    cl.AddCallback( Key( "double" ), &cmtk::ImageOperationConvertType::NewDouble, "64 bits floating point\n" );
    cl.EndGroup();
    
    cl.BeginGroup( "Flipping", "Image Flipping" );
    cl.AddCallback( Key( "flip-x" ), &cmtk::ImageOperationFlip::NewX, "Flip (mirror) along x-direction" );
    cl.AddCallback( Key( "flip-y" ), &cmtk::ImageOperationFlip::NewY, "Flip (mirror) along y-direction" );
    cl.AddCallback( Key( "flip-z" ), &cmtk::ImageOperationFlip::NewZ, "Flip (mirror) along z-direction" );
    cl.EndGroup();

    cl.BeginGroup( "MaskAndThreshold", "Image Masking and Thresholding" );
    cl.AddCallback( Key( "mask" ), &cmtk::ImageOperationApplyMask::New, "Binary mask file name: eliminate all image pixels where mask is 0." );
    cl.AddCallback( Key( "mask-inverse" ), &cmtk::ImageOperationApplyMask::NewInverse, "Inverse binary mask file name eliminate all image pixels where mask is NOT 0." );

    cl.AddCallback( Key( "thresh-below" ), &cmtk::ImageOperationThreshold::NewBelow, "Set all values below threshold to threshold value." );
    cl.AddCallback( Key( "thresh-above" ), &cmtk::ImageOperationThreshold::NewAbove, "Set all values above threshold to threshold value." );
    cl.AddCallback( Key( "thresh-below-to-padding" ), &cmtk::ImageOperationThreshold::NewBelowToPadding, "Set all values below threshold to padding value." );
    cl.AddCallback( Key( "thresh-above-to-padding" ), &cmtk::ImageOperationThreshold::NewAboveToPadding, "Set all values above threshold to padding value." );
    cl.EndGroup();

    cl.BeginGroup( "Intensity", "Intensity Transformations" );
    cl.AddCallback( Key( "scale-to-range" ), &cmtk::ImageOperationScaleToRange::New, "Scale image intensities to range 'from:to', e.g., '0:255' before conversion to byte data." );
    cl.AddCallback( Key( "histogram-equalization" ), &cmtk::ImageOperationHistogramEqualization::New, "Apply histogram equalization." );
    cl.AddCallback( Key( "histogram-equalization-nbins" ), &cmtk::ImageOperationHistogramEqualization::NewBins, "Apply histogram equalization with <int> number of bins." );

    cl.EndGroup();
    
    cl.BeginGroup( "Morphological", "Morphological Operations" );
    cl.AddCallback( Key( "erode" ), &cmtk::ImageOperationErodeDilate::NewErode, "Morphological erosion operator" );
    cl.AddCallback( Key( "dilate" ), &cmtk::ImageOperationErodeDilate::NewDilate, "Morphological dilation operator" );
    cl.AddCallback( Key( "boundary-map" ), &cmtk::ImageOperationBoundaryMap::New, "Create boundary map" );
    cl.AddCallback( Key( "multi-boundary-map" ), &cmtk::ImageOperationBoundaryMap::NewMulti, "Create multi-valued boundary map" );

    cl.AddCallback( Key( "distance-map" ), &cmtk::ImageOperationDistanceMap::NewUnsigned, "Compute unsigned Euclidean distance map. Input image is interpreted as binary mask." );
    cl.AddCallback( Key( "signed-distance-map" ), &cmtk::ImageOperationDistanceMap::NewSigned, "Compute signed (inside=negative, outside=positive) Euclidean distance map" );
    cl.AddCallback( Key( "medial-skeleton" ), &cmtk::ImageOperationMedialSkeleton::New, "Compute medial skeleton of binary mask" );
    cl.EndGroup();

    cl.BeginGroup( "Filtering", "Filter Operations" );
    cl.AddCallback( Key( "median-filter" ), &cmtk::ImageOperationMedianFilter::New, "Median filter. This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "gaussian-filter-sigma" ), &cmtk::ImageOperationGaussFilter::NewSigma, 
		    "Filter image with Gaussian kernel. This operation takes a single real-valued parameter, which specifies the kernel coefficient sigma in world units [e.g., mm] as the parameter." );
    cl.AddCallback( Key( "gaussian-filter-fwhm" ), &cmtk::ImageOperationGaussFilter::NewFWHM, 
		    "Filter image with Gaussian kernel. This operation takes a single real-valued parameter, which specifies the kernel full width at half maximum in world units [e.g., mm]." );
    cl.EndGroup();

    cl.BeginGroup( "Grid", "Grid Operations" );
    cl.AddCallback( Key( "downsample" ), &cmtk::ImageOperationDownsample::New, "Downsample image by factors 'Fx,Fy,Fz' or by single factor 'Fxyz'" );
    cl.AddCallback( Key( "crop-by-index" ), &cmtk::ImageOperationCropRegion::New, "Crop image to a region specified by a set of six grid index coordinates given as comma-separated integers x0,y0,z0,x1,y1,z1" );
    cl.AddCallback( Key( "crop-by-threshold" ), &cmtk::ImageOperationCropThreshold::New, "Crop image to region determined via a given threshold. "
		    "The resulting image will contain all pixels larger than the given parameter." );
    cl.AddCallback( Key( "crop-by-threshold-write-region" ), &cmtk::ImageOperationCropThreshold::NewWriteRegion, 
		    "Crop image to region determined via a given threshold and write cropping region to standard output." );
    cl.AddCallback( Key( "crop-by-threshold-write-xform" ), &cmtk::ImageOperationCropThreshold::NewWriteXform, 
		    "Crop image to region determined via a given threshold and write cropping transformation to standard output." );
    cl.EndGroup();

    cl.AddParameter( &imagePathIn, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &imagePathOut, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );

    if ( ! cl.Parse() ) 
      return 1;
    }
  catch ( cmtk::CommandLine::Exception ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( imagePathIn, Verbose ) );
  if ( ! volume ) 
    {
    cmtk::StdErr << "ERROR: could not read image " << imagePathIn << "\n";
    exit( 1 );
    }
  else
    {
    cmtk::TypedArray::SmartPtr volumeData = volume->GetData();
    if ( ! volumeData ) 
      {
      cmtk::StdErr << "ERROR: image seems to contain no data.\n";
      exit( 1 );
      }
    }
  
  if ( paddingDataFlag ) 
    {
    volume->GetData()->SetPaddingValue( paddingDataValue );
    }
  
  volume = cmtk::ImageOperation::ApplyAll( volume );

  if ( Verbose )
    {
    cmtk::StdErr << "Writing to file " << imagePathOut << "\n";
    }
  
  cmtk::VolumeIO::Write( volume, imagePathOut, Verbose );
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace convert
} // namespace apps
} // namespace cmtk
#endif
