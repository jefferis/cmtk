/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
#include <System/cmtkDebugOutput.h>
#include <System/cmtkTimers.h>

#include <Base/cmtkUniformVolume.h>
#include <IO/cmtkVolumeIO.h>

#include <IO/cmtkStudyList.h>
#include <IO/cmtkClassStreamStudyList.h>

#include <math.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkImageOperationConvertType.h>
#include <Base/cmtkImageOperationFlip.h>
#include <IO/cmtkImageOperationApplyMask.h>
#include <Base/cmtkImageOperationErodeDilate.h>
#include <Base/cmtkImageOperationBoundaryMap.h>
#include <Base/cmtkImageOperationConnectedComponents.h>
#include <Base/cmtkImageOperationDownsample.h>
#include <Base/cmtkImageOperationHistogramEqualization.h>
#include <Base/cmtkImageOperationHoughTransformSphere.h>
#include <Base/cmtkImageOperationCropRegion.h>
#include <Base/cmtkImageOperationCropThreshold.h>
#include <Base/cmtkImageOperationScaleToRange.h>
#include <Base/cmtkImageOperationThreshold.h>
#include <Base/cmtkImageOperationPruneHistogram.h>
#include <Base/cmtkImageOperationSetPadding.h>
#include <Base/cmtkImageOperationMapValues.h>
#include <Base/cmtkImageOperationReplace.h>
#include <Base/cmtkImageOperationRegionFilter.h>
#include <Base/cmtkImageOperationMedialSkeleton.h>
#include <Base/cmtkImageOperationGaussFilter.h>
#include <Base/cmtkImageOperationLaplaceFilter.h>
#include <Base/cmtkImageOperationDistanceMap.h>
#include <Base/cmtkImageOperationRevert.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include <map>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

int
doMain( const int argc, const char* argv[] )
{
  const char* imagePathIn = NULL;
  const char* imagePathOut = NULL;

#ifdef CMTK_USE_SQLITE
  const char* updateDB = NULL;
#endif

  try 
    {
    cmtk::CommandLine cl;  
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Convert between image file formats and data types." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool converts between image file formats and pixel data types. It can also apply simple, general-purpose image operations in the process. "
		       "An arbitrary number of operations can be specified on the command line, which will be applied exactly in the order given." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "convertx [options] infile outfile" );

    typedef cmtk::CommandLine::Key Key;
#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &updateDB, "Path to image/transformation database that should be updated with the newly created image." );
    cl.EndGroup();
#endif

    cl.BeginGroup( "Input", "Input Image Controls" );
    cl.AddCallback( Key( "set-padding" ), &cmtk::ImageOperationSetPadding::New, "Set padding value: all pixels in the input image that have this value will be ignored in all subsequent operations." );
    cl.AddCallback( Key( "unset-padding" ), &cmtk::ImageOperationSetPadding::NewUnset, "Unset padding value: for all subsequent operations, all pixels will be treated according to their value." );
    cl.EndGroup();

    cl.BeginGroup( "Conversion", "Data Type Conversion" );
    cl.AddCallback( Key( "char" ), &cmtk::ImageOperationConvertType::NewChar, "8 bits, signed integer" );
    cl.AddCallback( Key( "byte" ), &cmtk::ImageOperationConvertType::NewByte, "8 bits, unsigned integer" );
    cl.AddCallback( Key( "short" ), &cmtk::ImageOperationConvertType::NewShort, "16 bits, signed integer" );
    cl.AddCallback( Key( "ushort" ), &cmtk::ImageOperationConvertType::NewUShort, "16 bits, unsigned integer" );
    cl.AddCallback( Key( "int" ), &cmtk::ImageOperationConvertType::NewInt, "32 bits signed integer" );
    cl.AddCallback( Key( "uint" ), &cmtk::ImageOperationConvertType::NewUInt, "32 bits unsigned integer" );
    cl.AddCallback( Key( "float" ), &cmtk::ImageOperationConvertType::NewFloat, "32 bits floating point" );
    cl.AddCallback( Key( "double" ), &cmtk::ImageOperationConvertType::NewDouble, "64 bits floating point\n" );
    cl.EndGroup();
    
    cl.BeginGroup( "Mappings", "Value Mappings" );
    cl.AddCallback( Key( "map-values" ), &cmtk::ImageOperationMapValues::New, "Apply mapping function to pixel values. Mapping is defined as 'VAL0[,VAL1,...][:NEWVAL]' to map values VAL0, VAL1, etc. to new value NEWVAL. "
		    "If NEWVAL is not given, values are set to padding." );
    cl.AddCallback( Key( "map-values-only" ), &cmtk::ImageOperationMapValues::NewExclusive, "Apply mapping function to pixel values and replace unmapped pixels with padding. Multiple such mapping rules can be concatenated as "
		    "RULE0+RULE1[+...]; all concatenated rules will be applied simultaneously."
		    "Mapping is defined as 'VAL0[,VAL1,...][:NEWVAL]' to map values VAL0, VAL1, etc. to new value NEWVAL. If NEWVAL is not given, values are set to padding. Multiple such mapping rules can be concatenated as "
		    "RULE0+RULE1[+...]; all concatenated rules will be applied simultaneously." );
    cl.AddCallback( Key( "replace-padding" ), &cmtk::ImageOperationReplace::NewReplacePadding, "Replace padded pixel data with given value." );
    cl.AddCallback( Key( "replace-inf-nan" ), &cmtk::ImageOperationReplace::NewReplaceInfNaN, "Replace all infinite and not-a-number pixels with given value." );
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
    cl.AddCallback( Key( "binarize-thresh" ), &cmtk::ImageOperationThreshold::NewBinarize, "Set all values below threshold to 0, all values equal or above to 1." );
    cl.AddCallback( Key( "prune-histogram" ), &cmtk::ImageOperationPruneHistogram::New, "Threshold image by 'intensity histogram pruning', i.e., for given argument n [histogram bins] "
		    "determine thresholds such that the 1/n-th fraction of highest and lowest voxels are thresholded." );
    cl.AddCallback( Key( "prune-histogram-high" ), &cmtk::ImageOperationPruneHistogram::NewHigh, "Like '--prune-histograms', but only remove high intensities." );
    cl.AddCallback( Key( "prune-histogram-low" ), &cmtk::ImageOperationPruneHistogram::NewLow, "Like '--prune-histograms', but only remove low intensities." );
    cl.EndGroup();

    cl.BeginGroup( "Intensity", "Intensity Transformations" );
    cl.AddCallback( Key( "scale-to-range" ), &cmtk::ImageOperationScaleToRange::New, "Scale image intensities to range 'from:to', e.g., '0:255' before conversion to byte data." );
    cl.AddCallback( Key( "histogram-equalization" ), &cmtk::ImageOperationHistogramEqualization::New, "Apply histogram equalization." );
    cl.AddCallback( Key( "histogram-equalization-nbins" ), &cmtk::ImageOperationHistogramEqualization::NewBins, "Apply histogram equalization with <int> number of bins." );

    cl.EndGroup();
    
    cl.BeginGroup( "Morphological", "Morphological Operations" );
    cl.AddCallback( Key( "revert" ), &cmtk::ImageOperationRevert::New, "Revert a binary mask, i.e., exchange foreground and background." );
    cl.AddCallback( Key( "erode" ), &cmtk::ImageOperationErodeDilate::NewErode, "Morphological erosion operator" );
    cl.AddCallback( Key( "dilate" ), &cmtk::ImageOperationErodeDilate::NewDilate, "Morphological dilation operator" );
    cl.AddCallback( Key( "connected-components" ), &cmtk::ImageOperationConnectedComponents::New, "Create connected components map with regions numbered by decreasing component size" );
    cl.AddCallback( Key( "boundary-map" ), &cmtk::ImageOperationBoundaryMap::New, "Create boundary map" );
    cl.AddCallback( Key( "multi-boundary-map" ), &cmtk::ImageOperationBoundaryMap::NewMulti, "Create multi-valued boundary map" );

    cl.AddCallback( Key( "distance-map" ), &cmtk::ImageOperationDistanceMap::NewUnsigned, "Compute unsigned Euclidean distance map. Input image is interpreted as binary mask." );
    cl.AddCallback( Key( "signed-distance-map" ), &cmtk::ImageOperationDistanceMap::NewSigned, "Compute signed (inside=negative, outside=positive) Euclidean distance map" );
    cl.AddCallback( Key( "medial-skeleton" ), &cmtk::ImageOperationMedialSkeleton::New, "Compute n-dimensional medial skeleton of binary mask. Argument must be either 1 or 2." );
    cl.EndGroup();

    cl.BeginGroup( "Filtering", "Filter Operations" );
    cl.AddCallback( Key( "median-filter" ), &cmtk::ImageOperationRegionFilter::NewMedian, "Median filter. This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "mean-filter" ), &cmtk::ImageOperationRegionFilter::NewMean, "Regional mean filter. This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "fast-mean-filter" ), &cmtk::ImageOperationRegionFilter::NewFastMean, "Regional mean filter (fast, linear time implementation). This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "variance-filter" ), &cmtk::ImageOperationRegionFilter::NewVariance, "Regional variance filter. "
		    "This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "fast-variance-filter" ), &cmtk::ImageOperationRegionFilter::NewFastVariance, "Fast (linear-time) regional variance filter. "
		    "This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "third-moment-filter" ), &cmtk::ImageOperationRegionFilter::NewThirdMoment, "Regional third moment filter. "
		    "This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "standard-deviation-filter" ), &cmtk::ImageOperationRegionFilter::NewStandardDeviation, "Regional standard deviation filter. "
		    "This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "smoothness-filter" ), &cmtk::ImageOperationRegionFilter::NewSmoothness, "Regional 'smoothness' filter. "
		    "This operation takes the filter radius in pixels as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.AddCallback( Key( "gaussian-filter-sigma" ), &cmtk::ImageOperationGaussFilter::NewSigma, 
		    "Filter image with Gaussian kernel. This operation takes a single real-valued parameter, which specifies the kernel coefficient sigma in world units [e.g., mm] as the parameter." );
    cl.AddCallback( Key( "gaussian-filter-fwhm" ), &cmtk::ImageOperationGaussFilter::NewFWHM, 
		    "Filter image with Gaussian kernel. This operation takes a single real-valued parameter, which specifies the kernel full width at half maximum in world units [e.g., mm]." );
    cl.AddCallback( Key( "laplace-filter" ), &cmtk::ImageOperationLaplaceFilter::New, "Filter image with edge-enhancing Laplacian kernel." );
    cl.EndGroup();
    
    cl.BeginGroup( "Feature", "Feature Detection Operations" );
    cl.AddCallback( Key( "hough-sphere" ), &cmtk::ImageOperationHoughTransformSphere::New, "Hough transform to detect sphere with a given radius. "
		    "The radius in world units (typically mm) is provided as an argument to this operation." );
    cl.EndGroup();
    
    cl.BeginGroup( "Grid", "Grid Operations" );
    cl.AddCallback( Key( "downsample-select" ), &cmtk::ImageOperationDownsample::NewSelect, "Downsample image by pixel selection using per-axis factors 'Fx,Fy,Fz' or using single factor 'Fxyz'" );
    cl.AddCallback( Key( "downsample-average" ), &cmtk::ImageOperationDownsample::NewAverage, "Downsample image by averaging using per-axis factors 'Fx,Fy,Fz' or using single factor 'Fxyz'" );
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

    if ( ! cl.Parse( argc, argv ) ) 
      return 1;
    }
  catch ( const cmtk::CommandLine::Exception& ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( imagePathIn ) );
  cmtk::TypedArray::SmartPtr volumeData = volume->GetData();

  volume = cmtk::ImageOperation::ApplyAll( volume );

  cmtk::DebugOutput( 1 ) << "Writing to file " << imagePathOut << "\n";  
  cmtk::VolumeIO::Write( *volume, imagePathOut );
  return 0;
}

#include "cmtkSafeMain"
