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
#include <cmtkVolumeIO.h>

#include <cmtkBestDirectionOptimizer.h>
#include <cmtkEntropyMinimizationIntensityCorrectionFunctional.h>

#include <math.h>
#include <vector>
#include <algorithm>

bool Verbose = false;

const char* ImportBiasFieldAdd = NULL;
const char* ImportBiasFieldMul = NULL;

const char* FNameBiasFieldAdd = NULL;
const char* FNameBiasFieldMul = NULL;

float ThresholdForegroundMin = -FLT_MAX;
float ThresholdForegroundMax = FLT_MAX;
bool ThresholdForegroundFlag = false;
const char* FNameMaskImage = NULL;
float SamplingDensity = 0;
size_t NumberOfHistogramBins = 256;

const char* FNameInputImage = NULL;
const char* FNameOutputImage = NULL;
bool OutputFloatImage = false;

int PolynomialDegreeAdd = 0;
int PolynomialDegreeMul = 2;
bool IncrementalPolynomials = false;

cmtk::Types::Coordinate StepMax = 1.0;
cmtk::Types::Coordinate StepMin = 0.1;

bool EstimateDropOff = false;
bool LogIntensities = false;

int
main( const int argc, const char *argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "MR Image Intensity Bias Field Correction" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This program corrects intensity inhomogeneity artifacts in MR images using a bias field estimated via entropy minimization." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "BiasFieldCorrection" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Be verbose" );

    cl.AddSwitch( Key( 'L', "log-intensities" ), &LogIntensities, true, "Use log intensities for entropy estimation." );
    cl.AddOption( Key( 'm', "mask" ), &FNameMaskImage, "Binary mask image filename." );
    cl.AddOption( Key( 't', "thresh-min" ), &ThresholdForegroundMin, "Minimum intensity threshold for image foreground.", &ThresholdForegroundFlag );
    cl.AddOption( Key( 'T', "thresh-max" ), &ThresholdForegroundMax, "Minimum intensity threshold for image foreground.", &ThresholdForegroundFlag );
    cl.AddOption( Key( 's', "sampling-density" ), &SamplingDensity, "Pixel sampling density (default: 1.0; all pixels)" );
    cl.AddOption( Key( 'n', "num-bins" ), &NumberOfHistogramBins, "Number of histogram bins for entropy estimation [default: 256]" );

    cl.AddOption( Key( "step-max" ), &StepMax, "Maximum (initial) search step size." );
    cl.AddOption( Key( "step-min" ), &StepMin, "Minimum (final) search step size." );
    cl.AddOption( Key( 'A', "degree-add" ), &PolynomialDegreeAdd, "Polynomial degree for additive correction." );
    cl.AddOption( Key( 'M', "degree-mul" ), &PolynomialDegreeMul, "Polynomial degree for multiplicative correction." );
    cl.AddSwitch( Key( 'I', "incremental" ), &IncrementalPolynomials, true, "Incrementally increase polynomial degrees." );
    cl.AddSwitch( Key( 'd', "dropoff" ), &EstimateDropOff, true, "Estimate signal drop-off from coil coverage." );

    cl.AddOption( Key( "write-bias-add" ), &FNameBiasFieldAdd, "File name for output of additive bias field." );
    cl.AddOption( Key( "write-bias-mul" ), &FNameBiasFieldMul, "File name for output of multiplicative bias field." );
    cl.AddSwitch( Key( 'F', "write-float" ), &OutputFloatImage, true, "Write output image with floating point pixel data [default: input data type]." );

    cl.AddOption( Key( "import-bias-add" ), &ImportBiasFieldAdd, "Import additive bias field (disables optimization)." );
    cl.AddOption( Key( "import-bias-mul" ), &ImportBiasFieldMul, "Import multiplicative bias field (disables optimization)." );
    
    cl.Parse();

    FNameInputImage = cl.GetNext();
    FNameOutputImage = cl.GetNext();
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  cmtk::UniformVolume::SmartPtr inputImage( cmtk::VolumeIO::ReadOriented( FNameInputImage, Verbose ) );
  if ( ! inputImage || ! inputImage->GetData() )
    {
    cmtk::StdErr << "ERROR: Could not read input image " << FNameInputImage << "\n";
    exit( 1 );
    }

  cmtk::UniformVolume::SmartPtr maskImage;
  if ( ThresholdForegroundFlag )
    {
    maskImage = cmtk::UniformVolume::SmartPtr( inputImage->Clone( true /*copyData*/ ) );
    maskImage->GetData()->SetPaddingValue( 0.0 );
    maskImage->GetData()->ThresholdToPadding( ThresholdForegroundMin, ThresholdForegroundMax );
    }
  
  if ( FNameMaskImage )
    {
    maskImage = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( FNameMaskImage, Verbose ) );
    if ( ! maskImage || ! maskImage->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read mask image " << FNameMaskImage << "\n";
      exit( 1 );
      }
    }

  typedef cmtk::EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr FunctionalPointer;
  FunctionalPointer functional( NULL );
  
  cmtk::UniformVolume::SmartPtr outputImage;
  if ( ImportBiasFieldAdd || ImportBiasFieldMul )
    {
    functional = cmtk::CreateEntropyMinimizationIntensityCorrectionFunctional( 1, 1 );

    if ( SamplingDensity > 0 )
      functional->SetSamplingDensity( SamplingDensity );

    if ( NumberOfHistogramBins )
      functional->SetNumberOfHistogramBins( NumberOfHistogramBins );

    functional->SetInputImage( inputImage );
    
    if ( ImportBiasFieldAdd )
      {
      cmtk::UniformVolume::SmartPtr biasAdd( cmtk::VolumeIO::ReadOriented( ImportBiasFieldAdd, Verbose ) );
      if ( ! biasAdd || ! biasAdd->GetData() )
	{
	cmtk::StdErr << "ERROR: Could not read additive bias field image " << ImportBiasFieldAdd << "\n";
	exit( 1 );
	}
      functional->SetBiasFieldAdd( biasAdd );
      }

    if ( ImportBiasFieldMul )
      {
      cmtk::UniformVolume::SmartPtr biasMul( cmtk::VolumeIO::ReadOriented( ImportBiasFieldMul, Verbose ) );
      if ( ! biasMul || ! biasMul->GetData() )
	{
	cmtk::StdErr << "ERROR: Could not read multiplicative bias field image " << ImportBiasFieldMul << "\n";
	exit( 1 );
	}
      functional->SetBiasFieldMul( biasMul );
      }

    outputImage = functional->GetOutputImage( true /*update*/ );
    }
  else
    {
    cmtk::CoordinateVector v;
    size_t polynomialDegreeFrom = std::max( PolynomialDegreeAdd, PolynomialDegreeMul );
    const size_t polynomialDegreeTo = polynomialDegreeFrom + 1;
    if ( IncrementalPolynomials )
      polynomialDegreeFrom = 1;
    
    for ( size_t polynomialDegree = polynomialDegreeFrom; polynomialDegree < polynomialDegreeTo; ++polynomialDegree )
      {
      const size_t degreeAdd = std::min<size_t>( polynomialDegree, PolynomialDegreeAdd );
      const size_t degreeMul = std::min<size_t>( polynomialDegree, PolynomialDegreeMul );

      functional = cmtk::CreateEntropyMinimizationIntensityCorrectionFunctional( degreeAdd, degreeMul, functional );
      functional->SetUseLogIntensities( LogIntensities );

      if ( SamplingDensity > 0 )
	functional->SetSamplingDensity( SamplingDensity );
      
      if ( NumberOfHistogramBins )
	functional->SetNumberOfHistogramBins( NumberOfHistogramBins );
      
      functional->SetInputImage( inputImage );
      if ( maskImage && maskImage->GetData() )
	{
	functional->SetForegroundMask( maskImage );
	}
      else
	{
	cmtk::StdErr << "ERROR: please use a mask image. Seriously.\n";
	exit( 1 );
	}
      functional->GetParamVector( v );
      
      if ( Verbose )
	{
	cmtk::StdErr.printf( "Estimating bias field with order %d multiplicative / %d additive polynomials.\nNumber of parameters: %d\n", degreeMul, degreeAdd, v.Dim );
	}

      if ( (PolynomialDegreeAdd > 0) || (PolynomialDegreeMul > 0) )
	{
	try
	  {
	  cmtk::Optimizer::SmartPtr optimizer;
	  optimizer = cmtk::Optimizer::SmartPtr( new cmtk::BestDirectionOptimizer );
	  optimizer->SetFunctional( functional );
	  
	  optimizer->Optimize( v, StepMax, StepMin );
	  }
	catch ( const char* cp )
	  {
	  cmtk::StdErr << "EXCEPTION: " << cp << "\n";
	  }
	v.Print();
	}
      }

    // make sure everything is according to optimum parameters
    outputImage = functional->GetOutputImage( v, false/*foregroundOnly*/ );
    }
  
  if ( EstimateDropOff )
    {
    std::vector<float> edgeImageMean( outputImage->GetDims( 2 ) );
    for ( int i = 0; i < outputImage->GetDims( 2 ); ++i )
      {
      cmtk::ScalarImage::SmartPtr slice( outputImage->GetOrthoSlice( cmtk::AXIS_Z, i ) );
      slice->ApplySobel2DFilter();
      
      cmtk::Types::DataItem mean, variance;
      slice->GetPixelData()->GetStatistics( mean, variance );

      edgeImageMean[i] = mean;
      }

    const float maxOfMeans = *(std::max_element( edgeImageMean.begin(), edgeImageMean.end() ) );

    for ( int i = 0; i < outputImage->GetDims( 2 ); ++i )
      {
      const float scale = maxOfMeans / edgeImageMean[i];
      for ( int y = 0; y < outputImage->GetDims( 1 ); ++y )
	for ( int x = 0; x < outputImage->GetDims( 0 ); ++x )
	  {
	  const size_t offset = outputImage->GetOffsetFromIndex( x, y, i );
	  outputImage->SetDataAt( scale * outputImage->GetDataAt( offset ), offset );
	  }
      }
    }
  
  if ( FNameOutputImage )
    {
    if ( ! OutputFloatImage )
      {
      outputImage->GetData()->ReplacePaddingData( 0.0 );
      cmtk::TypedArray::SmartPtr convertedData( outputImage->GetData()->Convert( inputImage->GetData()->GetType() ) );
      outputImage->SetData( convertedData );
      }
    cmtk::VolumeIO::Write( outputImage, FNameOutputImage, Verbose );
    }

  if ( FNameBiasFieldAdd && PolynomialDegreeAdd )
    {
    cmtk::UniformVolume::SmartPtr biasField( functional->GetBiasFieldAdd( true /*completeImage*/ ) );
    cmtk::VolumeIO::Write( biasField, FNameBiasFieldAdd, Verbose );
    }

  if ( FNameBiasFieldMul && PolynomialDegreeMul )
    {
    cmtk::UniformVolume::SmartPtr biasField( functional->GetBiasFieldMul( true /*completeImage*/ ) );
    cmtk::VolumeIO::Write( biasField, FNameBiasFieldMul, Verbose );
    }

  return 0;
}
