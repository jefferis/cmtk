/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
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
#include <cmtkProgressConsole.h>

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>

#include <cmtkBestDirectionOptimizer.h>
#include <GPU/cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice.h>

#include <math.h>
#include <vector>
#include <algorithm>

#ifdef CMTK_USE_SQLITE
#  include <cmtkImageXformDB.h>
#endif

bool Verbose = false;

const char* FNameBiasFieldAdd = NULL;
const char* FNameBiasFieldMul = NULL;

float ThresholdForegroundMin = -FLT_MAX;
float ThresholdForegroundMax = FLT_MAX;
bool ThresholdForegroundFlag = false;
const char* FNameMaskImage = NULL;

bool UseSamplingDensity = false;
float SamplingDensity = 1.0;

unsigned int NumberOfHistogramBins = 256;

const char* FNameInputImage = NULL;
const char* FNameOutputImage = NULL;
bool OutputFloatImage = false;

int PolynomialDegreeAdd = 0;
int PolynomialDegreeMul = 2;
bool IncrementalPolynomials = false;

cmtk::Types::Coordinate StepMax = 1.0;
cmtk::Types::Coordinate StepMin = 0.1;

bool LogIntensities = false;

#ifdef CMTK_USE_SQLITE
const char* updateDB = NULL;
#endif

int
main( const int argc, const char *argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv, cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "MR Image Intensity Bias Field Correction on GPU using CUDA" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This program corrects intensity inhomogeneity artifacts in MR images using a bias field estimated via entropy minimization using GPU-accelerated computation." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Artifact Correction" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Be verbose" )->SetProperties( cmtk::CommandLine::PROPS_NOXML );
    cl.AddSwitch( Key( 'F', "write-float" ), &OutputFloatImage, true, "Write output image with floating point pixel data. If this is not given, the input data type is used." );

    cl.BeginGroup( "Bias Field", "Bias Field Parameterization" );
    cl.AddOption( Key( 'A', "degree-add" ), &PolynomialDegreeAdd, "Polynomial degree for additive correction." );
    cl.AddOption( Key( 'M', "degree-mul" ), &PolynomialDegreeMul, "Polynomial degree for multiplicative correction." );
    cl.AddSwitch( Key( 'I', "incremental" ), &IncrementalPolynomials, true, "Incrementally increase polynomial degrees." );
    cl.EndGroup();

    cl.BeginGroup( "Preprocessing", "Input Image Preprocessing" );
    cl.AddOption( Key( 'm', "mask" ), &FNameMaskImage, "Binary mask image filename." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_LABELS );
    cl.AddOption( Key( 't', "thresh-min" ), &ThresholdForegroundMin, "Minimum intensity threshold for image foreground.", &ThresholdForegroundFlag );
    cl.AddOption( Key( 'T', "thresh-max" ), &ThresholdForegroundMax, "Minimum intensity threshold for image foreground.", &ThresholdForegroundFlag );
    cl.EndGroup();

    cl.BeginGroup( "Entropy Estimation", "Entropy Estimation Settings" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddSwitch( Key( 'L', "log-intensities" ), &LogIntensities, true, "Use log intensities for entropy estimation." );
    cl.AddOption( Key( 's', "sampling-density" ), &SamplingDensity, "Image sampling density to use only subset of image pixels", &UseSamplingDensity );
    cl.AddOption( Key( 'n', "num-bins" ), &NumberOfHistogramBins, "Number of histogram bins for entropy estimation" );
    cl.EndGroup();

    cl.BeginGroup( "Optimization", "Optimization Algorithm Settings" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );;
    cl.AddOption( Key( "step-max" ), &StepMax, "Maximum (initial) search step size." );
    cl.AddOption( Key( "step-min" ), &StepMin, "Minimum (final) search step size." );
    cl.EndGroup();

    cl.BeginGroup( "Bias Field I/O", "Output of Bias Fields" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );;
    cl.AddOption( Key( "write-bias-add" ), &FNameBiasFieldAdd, "File name for output of additive bias field." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-bias-mul" ), &FNameBiasFieldMul, "File name for output of multiplicative bias field." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.EndGroup();
    
#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &updateDB, "Path to image/transformation database that should be updated with the newly created image." );
    cl.EndGroup();
#endif

    cl.AddParameter( &FNameInputImage, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &FNameOutputImage, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse();
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  // Instantiate programm progress indicator.
  cmtk::ProgressConsole progressIndicator( "Intensity Bias Field Correction" );

  cmtk::UniformVolume::SmartPtr inputImage( cmtk::VolumeIO::ReadOriented( FNameInputImage, Verbose ) );
  if ( ! inputImage || ! inputImage->GetData() )
    {
    cmtk::StdErr << "ERROR: Could not read input image " << FNameInputImage << "\n";
    exit( 1 );
    }

  cmtk::UniformVolume::SmartPtr maskImage;
  if ( FNameMaskImage )
    {
    maskImage = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( FNameMaskImage, Verbose ) );
    if ( ! maskImage || ! maskImage->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read mask image " << FNameMaskImage << "\n";
      exit( 1 );
      }
    }
  else
    {
    if ( ThresholdForegroundFlag )
      {
      maskImage = cmtk::UniformVolume::SmartPtr( inputImage->Clone( true /*copyData*/ ) );
      maskImage->GetData()->SetPaddingValue( 0.0 );
      maskImage->GetData()->ThresholdToPadding( cmtk::Types::DataItemRange( ThresholdForegroundMin, ThresholdForegroundMax ) );
      }  
    }

  typedef cmtk::EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr FunctionalPointer;
  FunctionalPointer functional( NULL );
  
  cmtk::UniformVolume::SmartPtr outputImage;

  cmtk::CoordinateVector v;
  size_t polynomialDegreeFrom = std::max( PolynomialDegreeAdd, PolynomialDegreeMul );
  const size_t polynomialDegreeTo = polynomialDegreeFrom + 1;
  if ( IncrementalPolynomials )
    polynomialDegreeFrom = 1;
  
  for ( size_t polynomialDegree = polynomialDegreeFrom; polynomialDegree < polynomialDegreeTo; ++polynomialDegree )
    {
    const size_t degreeAdd = std::min<size_t>( polynomialDegree, PolynomialDegreeAdd );
    const size_t degreeMul = std::min<size_t>( polynomialDegree, PolynomialDegreeMul );
    
    functional = cmtk::CreateEntropyMinimizationIntensityCorrectionFunctionalDevice( degreeAdd, degreeMul, functional );
    functional->SetUseLogIntensities( LogIntensities );
    
    if ( SamplingDensity > 0 )
      functional->SetSamplingDensity( SamplingDensity );
    
    if ( NumberOfHistogramBins )
      functional->SetNumberOfHistogramBins( NumberOfHistogramBins );
    
    functional->SetInputImage( inputImage );
    if ( maskImage && maskImage->GetData() )
      {
      functional->SetForegroundMask( *maskImage );
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
  
  if ( FNameOutputImage )
    {
    if ( ! OutputFloatImage )
      {
      outputImage->GetData()->ReplacePaddingData( 0.0 );
      cmtk::TypedArray::SmartPtr convertedData( outputImage->GetData()->Convert( inputImage->GetData()->GetType() ) );
      outputImage->SetData( convertedData );
      }
    cmtk::VolumeIO::Write( *outputImage, FNameOutputImage, Verbose );

#ifdef CMTK_USE_SQLITE
    if ( updateDB )
      {
      cmtk::ImageXformDB db( updateDB );
      db.AddImage( FNameOutputImage, FNameInputImage );
      }
#endif
    }

  if ( FNameBiasFieldAdd && PolynomialDegreeAdd )
    {
    cmtk::UniformVolume::SmartPtr biasField( functional->GetBiasFieldAdd( true /*completeImage*/ ) );
    cmtk::VolumeIO::Write( *biasField, FNameBiasFieldAdd, Verbose );

#ifdef CMTK_USE_SQLITE
    if ( updateDB )
      {
      cmtk::ImageXformDB db( updateDB );
      db.AddImage( FNameBiasFieldAdd, FNameInputImage );
      }
#endif
    }

  if ( FNameBiasFieldMul && PolynomialDegreeMul )
    {
    cmtk::UniformVolume::SmartPtr biasField( functional->GetBiasFieldMul( true /*completeImage*/ ) );
    cmtk::VolumeIO::Write( *biasField, FNameBiasFieldMul, Verbose );

#ifdef CMTK_USE_SQLITE
    if ( updateDB )
      {
      cmtk::ImageXformDB db( updateDB );
      db.AddImage( FNameBiasFieldMul, FNameInputImage );
      }
#endif
    }

  return 0;
}
