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

#include <cmtkDataGridMorphologicalOperators.h>
#include <cmtkDataGridFilter.h>
#include <cmtkUniformVolumeFilter.h>

#include <cmtkStudyList.h>
#include <cmtkClassStreamStudyList.h>

#include <math.h>
#include <cmtkMathFunctionWrappers.h>
#include <cmtkTypedArrayFunctionHistogramEqualization.h>

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

cmtk::ScalarDataType DataType = cmtk::TYPE_NONE;

std::set<cmtk::Types::DataItem> SelectValues;

void
Select( const char* arg ) 
{
  SelectValues.insert( atoi( arg ) ); 
}

std::map<cmtk::Types::DataItem,cmtk::Types::DataItem> MapValues;
cmtk::Types::DataItem MapFrom;
bool MapOnly = false;

void
Map( const char* arg ) 
{
  MapValues.insert( std::make_pair<cmtk::Types::DataItem,cmtk::Types::DataItem>( MapFrom, atoi( arg ) ) );
}

const char* InFileName = NULL;
const char* OutFileName = NULL;

const char* CropFileName = NULL;

cmtk::Types::DataItem ReplaceInfNaN = 0;
bool ReplaceInfNaNFlag = false;

bool CropImages = false;
int CropImagesRegionFrom[3] = { 0,0,0 };
int CropImagesRegionTo[3] = { 0,0,0 };

bool AutoCrop = false;
cmtk::Types::DataItem AutoCropThreshold = 0;
const char* CropXformOutFileName = NULL;

const char* Downsample = NULL;

bool AddGaussianNoise = 0;
cmtk::Types::Coordinate AddGaussianNoiseFWHM = 0;

cmtk::Types::DataItem
AddGaussianNoiseFunction( const cmtk::Types::DataItem value )
{
  return value + cmtk::MathUtil::NormalRandom( AddGaussianNoiseFWHM / sqrt( 8 * log( 2.0 ) ) );
}

cmtk::Types::DataItem
InvertFunction( const cmtk::Types::DataItem value )
{
  return 1.0 / value;
}

void
CallbackCropImages( const char* arg )
{
  CropImages = (6 == sscanf( arg, "%d,%d,%d,%d,%d,%d",
			     &CropImagesRegionFrom[0], &CropImagesRegionFrom[1], &CropImagesRegionFrom[2],
			     &CropImagesRegionTo[0], &CropImagesRegionTo[1], &CropImagesRegionTo[2] ) );
}

bool BoundaryMap = false;
bool BoundaryMapMultiValued = false;

std::list<int> ErodeOrDilate;

void
CallbackErode( const char* arg )
{
  ErodeOrDilate.push_back( -atoi( arg ) );
}

void
CallbackDilate( const char* arg )
{
  ErodeOrDilate.push_back( atoi( arg ) );
}

cmtk::Types::DataItem Threshold = 0;
bool HaveThreshold = false;
cmtk::Types::DataItem ThresholdMax = 0;
bool HaveThresholdMax = false;

bool ThresholdToPadding = false;

cmtk::Types::DataItem RescaleMin = 0;
bool HaveRescaleMin = false;
cmtk::Types::DataItem RescaleMax = 0;
bool HaveRescaleMax = false;

cmtk::Types::DataItem RescaleSlope = 1.0;
bool HaveRescaleSlope = false;

cmtk::Types::DataItem BinarizeThreshold = 0;
bool HaveBinarizeThreshold = false;

int MedianFilterRadius = 0;

bool Revert = false;
bool Invert = false;
bool HistogramEqualization = false;
int PruneHistogramBinsHigh = 0;

cmtk::Types::DataItem Gamma = 1.0;
bool ApplyGamma = false;

cmtk::Types::Coordinate GaussFilterSigma = 1.0;
bool ApplyGaussFilter = false;

bool ApplySobelFilter = false;

bool ApplyLog = false;
bool ApplyExp = false;
bool ApplySquare = false;
bool ApplySqrt = false;
bool ApplyAbs = false;

cmtk::Types::DataItem PaddingData = 0.0;
bool HavePaddingData = false;

cmtk::Types::DataItem ReplacePaddingData = 0.0;
bool HaveReplacePaddingData = false;

const char* MaskFileName = NULL;
bool MaskInverse = false;

bool FlipX = false;
bool FlipY = false;
bool FlipZ = false;

int EliminatePaddingVoting = 0;

int
main( int argc, char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );  
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Convert between image formats" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] infile outfile" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'N', "set-padding" ), &PaddingData, "Set Padding data for input image.", &HavePaddingData );
    cl.AddOption( Key( 'R', "replace-padding" ), &ReplacePaddingData, "Replace Padding data for input image.", &HaveReplacePaddingData );
    cl.AddOption( Key( "replace-inf-nan" ), &ReplaceInfNaN, "Replace all infinite and NaN values with given value.", &ReplaceInfNaNFlag );
    
    cl.AddSwitch( Key( 'c', "char" ), &DataType, cmtk::TYPE_CHAR, "8 bits, signed" );
    cl.AddSwitch( Key( 'b', "byte" ), &DataType, cmtk::TYPE_BYTE, "8 bits, unsigned" );
    cl.AddSwitch( Key( 's', "short" ), &DataType, cmtk::TYPE_SHORT, "16 bits, signed" );
    cl.AddSwitch( Key( 'u', "ushort" ), &DataType, cmtk::TYPE_USHORT, "16 bits, unsigned" );
    cl.AddSwitch( Key( 'i', "int" ), &DataType, cmtk::TYPE_INT, "32 bits signed" );
    cl.AddSwitch( Key( 'f', "float" ), &DataType, cmtk::TYPE_FLOAT, "32 bits floating point" );
    cl.AddSwitch( Key( 'd', "double" ), &DataType, cmtk::TYPE_DOUBLE, "64 bits floating point\n" );
    
    cl.AddSwitch( Key( "flip-x" ), &FlipX, true, "Flip (mirror) along x-direction" );
    cl.AddSwitch( Key( "flip-y" ), &FlipY, true, "Flip (mirror) along y-direction" );
    cl.AddSwitch( Key( "flip-z" ), &FlipZ, true, "Flip (mirror) along z-direction" );

    cl.AddSwitch( Key( 'B', "bmap" ), &BoundaryMap, true, "Create boundary map" );
    cl.AddSwitch( Key( "bmap-multi" ), &BoundaryMapMultiValued, true, "Create multi-valued boundary map" );
    cl.AddCallback( Key( 'E', "erode" ), CallbackErode, "Erosion iterations" );
    cl.AddCallback( Key( 'D', "dilate" ), CallbackDilate, "Dilation iterations" );
    cl.AddOption( Key( "eliminate-padding-voting" ), &EliminatePaddingVoting, "Eliminate padding-data pixels by [n] iterations of neighborhood voting." );
    cl.AddOption( Key( 't', "thresh" ), &Threshold, "Threshold value", &HaveThreshold );
    cl.AddOption( Key( "thresh-max" ), &ThresholdMax, "Maximum threshold value", &HaveThresholdMax );
    cl.AddSwitch( Key( "thresh-to-padding" ), &ThresholdToPadding, true, "Threshold to Padding (padding) instead of threshold values." );
    cl.AddOption( Key( 'T', "binarize-thresh" ), &BinarizeThreshold, "Binarize with threshold value", &HaveBinarizeThreshold );
    cl.AddOption( Key( 'm', "median" ), &MedianFilterRadius, "Median filter with given radius" );
    cl.AddOption( Key( 'G', "gauss-filter" ), &GaussFilterSigma, "Apply smoothing with Gaussian kernel [sigma].", &ApplyGaussFilter );
    cl.AddSwitch( Key( "sobel-filter" ), &ApplySobelFilter, true, "Apply Sobel edge detection filter." );
    cl.AddSwitch( Key( 'r', "revert" ), &Revert, true, "Revert value range" );
    cl.AddSwitch( Key( 'H', "hist-eq" ), &HistogramEqualization, true, "Histogram equalization." );
    cl.AddOption( Key( 'P', "prune-histogram" ), &PruneHistogramBinsHigh, "Prune high intensity range for use with <n> bin histograms." );

    cl.AddSwitch( Key( "abs" ), &ApplyAbs, true, "Apply fabs() function." );
    cl.AddSwitch( Key( 'l', "log" ), &ApplyLog, true, "Apply log() function." );
    cl.AddSwitch( Key( 'e', "exp" ), &ApplyExp, true, "Apply exp() function." );
    cl.AddSwitch( Key( "sqrt" ), &ApplySqrt, true, "Apply sqrt() function." );
    cl.AddSwitch( Key( "square" ), &ApplySquare, true, "Square values." );
    cl.AddSwitch( Key( "invert" ), &Invert, true, "Replace each value with (1.0 / value)" );

    cl.AddOption( Key( "add-gaussian-noise" ), &AddGaussianNoiseFWHM, "Add gaussian noise [parameter: FWHM]", &AddGaussianNoise );

    cl.AddOption( Key( 'g', "gamma" ), &Gamma, "Apply gamma correction.", &ApplyGamma );
    cl.AddOption( Key( 'M', "mask" ), &MaskFileName, "Binary mask file name: eliminate all image pixels where mask is 0." );
    cl.AddOption( Key( "mask-inverse" ), &MaskFileName, "Inverse binary mask file name eliminate all image pixels where mask is NOT 0.", &MaskInverse );
    cl.AddCallback( Key( "crop" ), CallbackCropImages, "Crop image: x0,y0,z0,x1,y1,z2" );
    cl.AddOption( Key( 'C', "cropfile" ), &CropFileName, "Crop with grid indexes from file." );
    cl.AddOption( Key( 'A', "auto-crop" ), &AutoCropThreshold, "Automatic cropping with [threshold value]", &AutoCrop );
    cl.AddOption( Key( "crop-xform-out" ), &CropXformOutFileName, "Output filename for shift transformation from original to cropped image" );

    cl.AddOption( Key( "downsample" ), &Downsample, "Downsample image by factors 'x,y,z' or by single factor 'xyz'" );

    cl.AddCallback( Key( 'S', "select" ), Select, "Select specific values (can be repeated)." );

    cl.AddOption( Key( "rescale-min" ), &RescaleMin, "Rescale minimum value", &HaveRescaleMin );
    cl.AddOption( Key( "rescale-max" ), &RescaleMax, "Rescale maximum value", &HaveRescaleMax );
    cl.AddOption( Key( "rescale-slope" ), &RescaleSlope, "Rescale by slope (factor)", &HaveRescaleSlope );

    cl.AddOption( Key( "map-from" ), &MapFrom, "Map from value (must be followed by --map-to and can be repeated)." );
    cl.AddCallback( Key( "map-to" ), Map, "Map to value (alternates with --map-from and can be repeated)." );
    cl.AddSwitch( Key( "map-only" ), &MapOnly, true,"Output only successfully mapped values; set all others to zero" );
    
    if ( ! cl.Parse() ) return 1;
    
    InFileName = cl.GetNext();
    OutFileName = cl.GetNext();
    }
  catch ( cmtk::CommandLine::Exception ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InFileName, Verbose ) );
  if ( ! volume ) 
    {
    cmtk::StdErr << "ERROR: could not read image " << InFileName << "\n";
    exit( 1 );
    }
  
  cmtk::TypedArray::SmartPtr volumeData = volume->GetData();
  if ( ! volumeData ) 
    {
    cmtk::StdErr << "Image seems to contain no data.\n";
    exit( 1 );
    }
  
  if ( HavePaddingData ) 
    {
    volumeData->SetPaddingValue( PaddingData );
    }
  
  if ( HaveReplacePaddingData ) 
    {
    volumeData->ReplacePaddingData( ReplacePaddingData );
    }

  // note: we need to either do this before cropping, or crop the
  // mask file as well.
  if ( MaskFileName ) 
    {
    cmtk::UniformVolume::SmartPtr maskVolume( cmtk::VolumeIO::ReadOriented( MaskFileName, Verbose ) );
    if ( maskVolume && maskVolume->GetData() ) 
      {
      if ( Verbose )
	cmtk::StdErr << "Applying mask from file.\n";
      
      const std::string maskOrientation = maskVolume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
      const std::string inputImageOrientation = volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
      if ( maskOrientation != inputImageOrientation )
	{
	if ( Verbose )
	  {
	  cmtk::StdErr << "INFO: reorienting mask from '" << maskOrientation << "' to '" << inputImageOrientation << "' to match input image.\n";
	  }
	maskVolume = cmtk::UniformVolume::SmartPtr( maskVolume->GetReoriented( inputImageOrientation.c_str() ) );
	}
      
      const cmtk::TypedArray* maskData = maskVolume->GetData();
      for ( size_t i = 0; i < volumeData->GetDataSize(); ++i )
	if ( maskData->IsPaddingOrZeroAt( i ) != MaskInverse ) 
	  volumeData->SetPaddingAt( i );
      }
    else
      {
      cmtk::StdErr << "ERROR: could not read mask from file " << MaskFileName << "\nProgram will terminate now, just to be safe.\n";
      exit( 1 );
      }
    }

  if ( CropFileName ) 
    {
    std::ifstream cropStream( CropFileName );
    int cropFrom[3], cropTo[3];
    cropStream >> cropFrom[0] >> cropFrom[1] >> cropFrom[2] >> cropTo[0] >> cropTo[1] >> cropTo[2];
    
    if ( Verbose )
      cmtk::StdErr.printf( "Cropping between [%d,%d,%d] and [%d,%d,%d].\n", cropFrom[0], cropFrom[1], cropFrom[2], cropTo[0], cropTo[1], cropTo[2] );
    
    volume->SetCropRegion( cropFrom, cropTo );
    volume = cmtk::UniformVolume::SmartPtr( volume->GetCroppedVolume() );
    volumeData = volume->GetData();
    }
  else
    {
    if ( CropImages )
      {
      volume->SetCropRegion( CropImagesRegionFrom, CropImagesRegionTo );
      volume = cmtk::UniformVolume::SmartPtr( volume->GetCroppedVolume() );
      volumeData = volume->GetData();
      }
    }

  if ( AutoCrop )
    {
    volume->AutoCrop( AutoCropThreshold, true /*recrop*/ );

    if ( Verbose )
      {
      int cropFrom[3], cropTo[3];
      volume->GetCropRegion( cropFrom, cropTo );
      printf( "AutoCrop %d,%d,%d,%d,%d,%d\n", cropFrom[0], cropFrom[1], cropFrom[2], cropTo[0], cropTo[1], cropTo[2] );
      }

    if ( CropXformOutFileName )
      {
      cmtk::Types::Coordinate cropFrom[3], cropTo[3];
      volume->GetCropRegion( cropFrom, cropTo );

      cmtk::StudyList slist;
      cmtk::Study::SmartPtr refstudy = slist.AddStudy( OutFileName );
      cmtk::Study::SmartPtr fltstudy = slist.AddStudy( InFileName );
      
      const cmtk::Types::Coordinate parameters[15] = { cropFrom[0], cropFrom[1], cropFrom[2], 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 };
      
      cmtk::AffineXform::SmartPtr affineXform( new cmtk::AffineXform( parameters ) );
      slist.AddXform( refstudy, fltstudy, affineXform );
      cmtk::ClassStreamStudyList::Write( CropXformOutFileName, &slist );
      }

    volume = cmtk::UniformVolume::SmartPtr( volume->GetCroppedVolume() );
    volumeData = volume->GetData();
    }

if ( Downsample )
  {
  int factors[3] = { 1, 1, 1 };
  const size_t nFactors = sscanf( Downsample, "%d,%d,%d", factors, factors+1, factors+2 );
  if ( nFactors == 1 )
    {
    factors[1] = factors[2] = factors[0];
    }
  else
    {
    if ( nFactors != 3 )
      {
      cmtk::StdErr << "ERROR: downsampling factors must either be three integers, x,y,z, or a single integer\n";
      exit( 1 );
      }
    }
  volume = cmtk::UniformVolume::SmartPtr( volume->GetDownsampled( factors ) );
  volumeData = volume->GetData();
  }

  if ( ReplaceInfNaNFlag )
    {
    cmtk::Types::DataItem value = 0;
    for ( size_t i = 0; i < volumeData->GetDataSize(); ++i )
      {
      if ( volumeData->Get( value, i ) )
	{
	if ( !finite( value ) )
	  {
	  volumeData->Set( ReplaceInfNaN, i );
	  }
	}
      }
    }
  
  if ( EliminatePaddingVoting )
    {
    cmtk::DataGridMorphologicalOperators ops( volume );
    ops.EliminatePaddingVoting();
    }
  
  if ( HaveThreshold || HaveThresholdMax ) 
    {
    if ( Verbose )
      {
      if ( HaveThreshold )
	cmtk::StdErr.printf( "Thresholding at min level %f.\n", Threshold );
      if ( HaveThresholdMax )
	cmtk::StdErr.printf( "Thresholding at max level %f.\n", ThresholdMax );
      }
    
    cmtk::Types::DataItem min, max;
    volumeData->GetRange( min, max );
    
    if ( ! HaveThreshold ) 
      Threshold = min;
    if ( ! HaveThresholdMax ) 
      ThresholdMax = max;
    
    if ( ThresholdToPadding )
      volumeData->ThresholdToPadding( Threshold, ThresholdMax );
    else
      volumeData->Threshold( Threshold, ThresholdMax );
    //    volume->FillCropBackground( 100 );
    }
  
  if ( HaveRescaleSlope ) 
    {
    if ( Verbose )
      {
      if ( HaveRescaleSlope )
	cmtk::StdErr.printf( "Rescaling by factor %f.\n", RescaleSlope );
      }
    
    volumeData->Rescale( RescaleSlope, 0.0 );
    }
  
  if ( HaveRescaleMin || HaveRescaleMax ) 
    {
    if ( Verbose )
      {
      if ( HaveRescaleMin )
	cmtk::StdErr.printf( "Rescaling to min level %f.\n", RescaleMin );
      if ( HaveRescaleMax )
	cmtk::StdErr.printf( "Rescaling to max level %f.\n", RescaleMax );
      }

    cmtk::Types::DataItem min, max;
    volumeData->GetRange( min, max );

    if ( ! HaveRescaleMin ) 
      RescaleMin = min;
    if ( ! HaveRescaleMax ) 
      RescaleMax = max;

    const cmtk::Types::DataItem factor = (RescaleMax - RescaleMin) / (max-min);
    const cmtk::Types::DataItem offset = RescaleMin - min * factor;
    volumeData->Rescale( factor, offset );
    }
  
  if ( HaveBinarizeThreshold ) 
    {
    if ( Verbose )
      cmtk::StdErr.printf( "Binarizing with threshold level %f.\n", BinarizeThreshold );
    volumeData->Binarize( BinarizeThreshold );
    }
  
  if ( FlipX )
    {
    volume->ApplyMirrorPlane( cmtk::AXIS_X );
    }
  if ( FlipY )
    {
    volume->ApplyMirrorPlane( cmtk::AXIS_Y );
    }
  if ( FlipZ )
    {
    volume->ApplyMirrorPlane( cmtk::AXIS_Z );
    }
  
  if ( MedianFilterRadius ) 
    {
    if ( Verbose )
      cmtk::StdErr.printf( "Median filter with radius %d.\n", MedianFilterRadius );

    volume->SetData( cmtk::DataGridFilter( volume ).GetDataMedianFiltered( MedianFilterRadius ) );
    }
  
  if ( ApplyGaussFilter ) 
    {
    if ( Verbose )
      cmtk::StdErr.printf( "Gaussian filter with sigma = %f [mm].\n", GaussFilterSigma );

    volume->SetData( cmtk::UniformVolumeFilter( volume ).GetDataGaussFiltered( GaussFilterSigma ) );
    }
  
  if ( ApplySobelFilter ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Sobel edge detection filter.\n";

    volume->SetData( cmtk::DataGridFilter( volume ).GetDataSobelFiltered() );
    }
  
  if ( ApplyAbs ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Applying fabs() function.\n";
    volumeData->ApplyFunctionDouble( cmtk::Wrappers::Abs );
    }
  
  if ( ApplyLog ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Applying log() function.\n";
    volumeData->ApplyFunctionDouble( cmtk::Wrappers::Log );
    }
  
  if ( ApplyExp ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Applying exp() function.\n";
    volumeData->ApplyFunctionDouble( cmtk::Wrappers::Exp );
    }
  
  if ( ApplySqrt ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Applying sqrt() function.\n";
    volumeData->ApplyFunctionDouble( cmtk::Wrappers::Sqrt );
    }
  
  if ( ApplySquare ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Applying square function.\n";
    volumeData->ApplyFunctionDouble( cmtk::MathUtil::Square<double> );
    }

  if ( Invert ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Applying 1/x function.\n";
    volumeData->ApplyFunction( InvertFunction );
    }

  if ( PruneHistogramBinsHigh )
    {
    volumeData->PruneHistogram( true /*pruneHi*/, false /*pruneLo*/, PruneHistogramBinsHigh );
    }
  
  if ( HistogramEqualization ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Histogram equalization.\n";
    volumeData->ApplyFunctionObject( cmtk::TypedArrayFunctionHistogramEqualization( *volumeData ) );
    }
  
  if ( Revert ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Inverting value range.\n";
    cmtk::Types::DataItem min, max;
    volume->GetData()->GetRange( min, max );
    volume->GetData()->Rescale( -1, max + min );
    }
  
  if ( ApplyGamma ) 
    {
    if ( Verbose )
      cmtk::StdErr.printf( "Applying gamma correction with gamma = %f.\n", Gamma );
    volume->GetData()->GammaCorrection( Gamma );
    }
  
  if ( DataType == cmtk::TYPE_NONE )
    DataType = volumeData->GetType();

  if ( !SelectValues.empty() ) 
    {
    if ( Verbose )
      cmtk::StdErr.printf( "Selecting %d given values.\n", SelectValues.size() );
    
    size_t dataSize = volumeData->GetDataSize();
    cmtk::TypedArray::SmartPtr selectData( cmtk::TypedArray::Create( cmtk::ScalarDataType( DataType ), dataSize ) );
    
    cmtk::Types::DataItem item;
    for ( size_t i = 0; i < dataSize; ++i ) 
      {
      if ( volumeData->Get( item, i ) )
	if ( SelectValues.find( item ) != SelectValues.end() )
	  selectData->Set( item, i );
	else
	  selectData->Set( 0, i );
      else
	selectData->Set( 0, i );
      }
    
    volume->SetData( selectData );
    volumeData = selectData;
    }
  
  if ( MapValues.size() ) 
    {
    if ( Verbose )
      cmtk::StdErr << "Mapping values according to substitution table.\n";
    
    size_t dataSize = volumeData->GetDataSize();
    cmtk::TypedArray::SmartPtr mapData( cmtk::TypedArray::Create( cmtk::ScalarDataType( DataType ), dataSize ) );
    
    cmtk::Types::DataItem item;
    for ( size_t i = 0; i < dataSize; ++i ) {
    if ( volumeData->Get( item, i ) ) 
      {
      std::map<cmtk::Types::DataItem,cmtk::Types::DataItem>::const_iterator it = MapValues.find( item );
      if ( it != MapValues.end() )
      mapData->Set( it->second, i );
      else
	if ( MapOnly ) 
	  mapData->Set( 0, i );
	else
	  mapData->Set( item, i );
      } 
    else
      mapData->Set( 0, i );
    }
    
    volume->SetData( mapData );
    volumeData = mapData;
    }
  
  if ( BoundaryMap || BoundaryMapMultiValued )
    {
    cmtk::DataGridMorphologicalOperators ops( volume );
    volume->SetData( ops.GetBoundaryMap( BoundaryMapMultiValued ) );
    }

  switch ( DataType ) 
    {
    case cmtk::TYPE_CHAR:
    case cmtk::TYPE_BYTE:
    case cmtk::TYPE_SHORT:
    case cmtk::TYPE_USHORT:
    case cmtk::TYPE_INT:
    case cmtk::TYPE_FLOAT:
    case cmtk::TYPE_DOUBLE:
      if ( DataType != volumeData->GetType() ) 
	{
	if ( Verbose )
	  cmtk::StdErr << "Converting to new data type.\n";
	
	volume->SetData( cmtk::TypedArray::SmartPtr( volumeData->Convert( cmtk::ScalarDataType( DataType ) ) ) );
	}
      break;
    default:
      break;
    }
  
  if ( AddGaussianNoise )
    {
    srandom( static_cast<unsigned int>( cmtk::Timers::GetTimeProcess() ) );
    volumeData->ApplyFunction( AddGaussianNoiseFunction );
    }
  
  std::list<int>::const_iterator it = ErodeOrDilate.begin();
  for ( ; it != ErodeOrDilate.end(); ++it )
    {
    if ( *it < 0 ) 
      {
      if ( Verbose )
	{
	cmtk::StdErr << "Eroding by " << -(*it) << " pixels.\n";
	}
      cmtk::DataGridMorphologicalOperators ops( volume );
      volume->SetData( ops.GetEroded( -(*it) ) );
      }
    else
      {
      if ( Verbose )
	{
	cmtk::StdErr << "Dilating by " << (*it) << " pixels.\n";
	}
      cmtk::DataGridMorphologicalOperators ops( volume );
      volume->SetData( ops.GetDilated( (*it) ) );
      }
    }

  if ( Verbose )
    {
    cmtk::StdErr << "Writing to file " << OutFileName << "\n";
    }
  
  cmtk::VolumeIO::Write( volume, OutFileName, Verbose );
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace convert
} // namespace apps
} // namespace cmtk
#endif
