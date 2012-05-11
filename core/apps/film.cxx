/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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
#include <System/cmtkProgressConsole.h>

#include <Registration/cmtkAffineRegistration.h>
#include <Registration/cmtkProtocolCallback.h>

#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkCubicInterpolator.h>
#include <Base/cmtkSincInterpolator.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkVector3D.h>

#include <Recon/cmtkInverseInterpolationVolumeReconstruction.h>

#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkClassStreamAffineXform.h>
#include <IO/cmtkVolumeIO.h>

#include <algorithm>
#include <map>
#include <vector>

const char* InputFilePath = NULL;
const char* OutputFilePath = NULL;

cmtk::Types::DataItem PaddingValue = 0;
bool PaddingFlag = false;

int InterleaveAxis = -1;
unsigned int NumberOfPasses = 2;

int RegistrationMetric = 4; // MSD

double InjectionKernelSigma = 0.5;
double InjectionKernelRadius = 2;

bool FourthOrderError = false;

int InverseInterpolationKernel = cmtk::Interpolators::CUBIC;
int NumberOfIterations = 20;
bool RegionalIntensityTruncation = true;
double ConstraintWeightLNorm = 0;

const char* ReferenceImagePath = NULL;
const char* InjectedImagePath = NULL;

const char* ExportXformsPath = NULL;
const char* ImportXformsPath = NULL;

std::map<size_t,float> PassWeights;

bool WriteImagesAsFloat = false;

void
CallbackSetPassWeight( const char* argv )
{
  int pass = 0;
  float weight = 1.0;
  if ( 2 == sscanf( argv, "%d:%f", &pass, &weight ) )
    {
    PassWeights[pass] = weight;
    }
  else
    {
    cmtk::StdErr << "ERROR: pass weights must be given as 'pass:weight', where 'pass' is an integer and 'weight' is a number between 0 and 1.\n"
	      << "       Parameter provided was '" << argv << "'\n";
    throw cmtk::ExitException( 1 );
    }
}

template<class TInterpolator>
cmtk::UniformVolume::SmartPtr
GetReconstructedImage( cmtk::UniformVolume::SmartPtr& volume, cmtk::UniformVolume::SmartPtr& refImage, std::vector<cmtk::Xform::SmartPtr>& xformsToPassImages )
{
  if ( InterleaveAxis < 0 )
    InterleaveAxis = cmtk::VolumeInjectionReconstruction::GuessInterleaveAxis( volume );
  
  cmtk::InverseInterpolationVolumeReconstruction<TInterpolator> volRecon( volume, NumberOfPasses, InterleaveAxis );
  
  for ( std::map<size_t,float>::const_iterator it = PassWeights.begin(); it != PassWeights.end(); ++it )
    {
    volRecon.SetPassWeight( it->first, it->second );
    }
  
  if ( refImage )
    volRecon.SetReferenceImage( refImage );

  volRecon.SetUseRegionalIntensityTruncation( RegionalIntensityTruncation );
  volRecon.SetUseFourthOrderError( FourthOrderError );
  volRecon.SetConstraintWeightLNorm( ConstraintWeightLNorm );

  if ( xformsToPassImages.size() == NumberOfPasses )
    {
    volRecon.SetTransformationsToPassImages( xformsToPassImages );
    }
  else
    {
    cmtk::DebugOutput( 2 ) << "Computing transformations between passes...\n";
    volRecon.ComputeTransformationsToPassImages( RegistrationMetric );
    xformsToPassImages = volRecon.GetTransformationsToPassImages();
    }
  
  if ( ExportXformsPath )
    {
    cmtk::ClassStreamOutput stream( ExportXformsPath, cmtk::ClassStreamOutput::MODE_WRITE );
    if ( stream.IsValid() )
      {
      cmtk::DebugOutput( 2 ) << "Exporting transformations between passes to " << ExportXformsPath << "\n";
      for ( unsigned int pass = 0; pass < NumberOfPasses; ++pass )
	{
	stream << dynamic_cast<cmtk::AffineXform&>( *xformsToPassImages[pass] );
	}
      }
    else
      {
      cmtk::StdErr << "ERROR: Could not open transformation file" << ExportXformsPath << "\n";
      }
    }

  cmtk::DebugOutput( 2 ) << "Volume injection...\n";
  volRecon.VolumeInjectionAnisotropic( InjectionKernelSigma, InjectionKernelRadius );
  if ( InjectedImagePath )
    {
    cmtk::UniformVolume::SmartPtr outputImage = volRecon.GetCorrectedImage();
    if ( !WriteImagesAsFloat && outputImage->GetData()->GetType() != volume->GetData()->GetType() )
      {
      outputImage = cmtk::UniformVolume::SmartPtr( outputImage->CloneGrid() );
      outputImage->SetData( cmtk::TypedArray::SmartPtr( volRecon.GetCorrectedImage()->GetData()->Convert( volume->GetData()->GetType() ) ) );
      }
    cmtk::VolumeIO::Write( *outputImage, InjectedImagePath );
    }

  volRecon.Optimize( NumberOfIterations );
  return volRecon.GetCorrectedImage();
}  

int
doMain( const int argc, const char* argv[] )
{
  /*
  // Parse command line
  */
  try
    {
    cmtk::CommandLine cl( cmtk::CommandLine::PROPS_XML );
    typedef cmtk::CommandLine::Key Key;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Fix interleaved motion using inverse interpolation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool splits an interleaved input image into the pass images, co-registers them, and reconstructs a motion-corrected image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Artifact Correction" );

    cl.BeginGroup( "input", "Input Options" );
    cl.AddOption( Key( "padding-value" ), &PaddingValue, "Set padding value for input image. Pixels with this value will be ignored.", &PaddingFlag );
    cl.EndGroup();

    cl.BeginGroup( "interleave", "Interleaving Options" );
    cmtk::CommandLine::EnumGroup<int>::SmartPtr
      interleaveGroup = cl.AddEnum( "interleave-axis", &InterleaveAxis, "Define interleave axis: this is the through-slice direction of the acquisition." );
    interleaveGroup->AddSwitch( Key( "guess-from-input" ), -1, "Guess from input image" );
    interleaveGroup->AddSwitch( Key( 'a', "axial" ), (int)cmtk::AXIS_Z, "Interleaved axial images" );
    interleaveGroup->AddSwitch( Key( 's', "sagittal" ),(int)cmtk::AXIS_X, "Interleaved sagittal images" );
    interleaveGroup->AddSwitch( Key( 'c', "coronal" ), (int)cmtk::AXIS_Y, "Interleaved coronal images" );
    interleaveGroup->AddSwitch( Key( 'x', "interleave-x" ), (int)cmtk::AXIS_X, "Interleaved along x axis" );
    interleaveGroup->AddSwitch( Key( 'y', "interleave-y" ), (int)cmtk::AXIS_Y, "Interleaved along y axis" );
    interleaveGroup->AddSwitch( Key( 'z', "interleave-z" ), (int)cmtk::AXIS_Z, "Interleaved along z axis" );

    cl.AddOption( Key( 'p', "passes" ), &NumberOfPasses, "Number of interleaved passes" );
    cl.AddCallback( Key( 'W', "pass-weight" ), CallbackSetPassWeight, "Set contribution weight for a pass in the form 'pass:weight'" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
;
    cl.EndGroup();

    cl.BeginGroup( "motion", "Motion Correction / Registration Options" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( 'R', "reference-image" ), &ReferenceImagePath, "Use a separate high-resolution reference image for registration" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

    cmtk::CommandLine::EnumGroup<int>::SmartPtr
      metricGroup = cl.AddEnum( "registration-metric", &RegistrationMetric, "Registration metric for motion estimation by image-to-image registration." );
    metricGroup->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    metricGroup->AddSwitch( Key( "nmi" ), 0, "Use Normalized Mutual Information for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "mi" ), 1, "Use standard Mutual Information for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "cr" ), 2, "Use Correlation Ratio for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "msd" ), 4, "Use Mean Squared Differences for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "cc" ), 5, "Use Cross-Correlation for pass-to-refereence registration" );

    cl.AddOption( Key( "import-xforms-path" ), &ImportXformsPath, "Path of file from which to import transformations between passes." )
      ->SetProperties( cmtk::CommandLine::PROPS_FILENAME | cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( "export-xforms-path" ), &ExportXformsPath, "Path of file to which to export transformations between passes." )
      ->SetProperties( cmtk::CommandLine::PROPS_FILENAME | cmtk::CommandLine::PROPS_ADVANCED | cmtk::CommandLine::PROPS_OUTPUT );

    cl.BeginGroup( "inject", "Initial Volume Injection Options" );
    cl.AddOption( Key( 'S', "injection-kernel-sigma" ), &InjectionKernelSigma, "Standard deviation of Gaussian kernel for volume injection in multiples of pixel size in each direction." );
    cl.AddOption( Key( 'r', "injection-kernel-radius" ), &InjectionKernelRadius, "Truncation radius factor of injection kernel. The kernel is truncated at sigma*radius, where sigma is the kernel standard deviation." );
    cl.EndGroup();

    cl.BeginGroup( "invint", "Inverse Interpolation Options" );
    cmtk::CommandLine::EnumGroup<int>::SmartPtr kernelGroup = 
      cl.AddEnum( "inverse-interpolation-kernel", &InverseInterpolationKernel, "Kernel for the inverse interpolation reconstruction" );
    kernelGroup->AddSwitch( Key( 'C', "cubic" ), cmtk::Interpolators::CUBIC, "Tricubic interpolation" );
    kernelGroup->AddSwitch( Key( 'L', "linear" ), cmtk::Interpolators::LINEAR, "Trilinear interpolation (faster but less accurate)" );
    kernelGroup->AddSwitch( Key( 'H', "hamming-sinc" ), cmtk::Interpolators::HAMMING_SINC, "Hamming-windowed sinc interpolation" );
    kernelGroup->AddSwitch( Key( 'O', "cosine-sinc" ), cmtk::Interpolators::COSINE_SINC, "Cosine-windowed sinc interpolation (most accurate but slowest)" );

    cl.AddSwitch( Key( 'f', "fourth-order-error" ), &FourthOrderError, true, "Use fourth-order (rather than second-order) error for optimization." )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( 'n', "num-iterations" ), &NumberOfIterations, "Maximum number of inverse interpolation iterations" );
    cl.EndGroup();

    cl.BeginGroup( "regularize", "Reconstruction Regularization Options" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( "l-norm-weight" ), &ConstraintWeightLNorm, "Set constraint weight for Tikhonov-type L-Norm regularization (0 disables constraint)" );
    cl.AddSwitch( Key( 'T', "no-truncation" ), &RegionalIntensityTruncation, false, "Turn off regional intensity truncatrion" );
    cl.EndGroup();
    
    cl.BeginGroup( "output", "Output Options" );
    cl.AddOption( Key( "write-injected-image" ), &InjectedImagePath, "Write initial volume injection image to path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddSwitch( Key( 'F', "write-images-as-float" ), &WriteImagesAsFloat, true, "Write output images as floating point [default: same as input]" );
    cl.EndGroup();

    cl.AddParameter( &InputFilePath, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &OutputFilePath, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  // Instantiate programm progress indicator.
  cmtk::ProgressConsole progressIndicator( "Interleaved Image Motion Correction" );

  /*
  // Read input image
  */
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InputFilePath ) );
  if ( ! volume || ! volume->GetData() )
    {
    cmtk::StdErr << "ERROR: Could not read image " << InputFilePath << "\n";
    return 1;
    }

  if ( PaddingFlag )
    {
    volume->GetData()->SetPaddingValue( PaddingValue );
    }

  cmtk::UniformVolume::SmartPtr refImage;
  if ( ReferenceImagePath )
    {
    refImage = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( ReferenceImagePath ) );
    if ( ! refImage || ! refImage->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read image " << ReferenceImagePath << "\n";
      return 1;
      }
    }

  std::vector<cmtk::Xform::SmartPtr> xformsToPassImages;
  if ( ImportXformsPath )
    {
    cmtk::ClassStreamInput stream( ImportXformsPath );
    if ( stream.IsValid() )
      {
      cmtk::DebugOutput( 1 ) << "Importing transformations between passes from " << ImportXformsPath << "\n";

      cmtk::AffineXform xform;
      for ( unsigned int pass = 0; pass < NumberOfPasses; ++pass )
	{
	stream >> xform;
	xformsToPassImages.push_back( cmtk::Xform::SmartPtr( xform.Clone() ) );
	}
      }
    else
      {
      cmtk::StdErr << "ERROR: Could not open transformation file" << ImportXformsPath << "\n";
      }
    }

  cmtk::UniformVolume::SmartPtr correctedVolume;
  switch ( InverseInterpolationKernel )
    {
    default:
    case cmtk::Interpolators::LINEAR: 
      correctedVolume = GetReconstructedImage<cmtk::Interpolators::Linear>( volume, refImage, xformsToPassImages );
      break;
    case cmtk::Interpolators::CUBIC: 
      correctedVolume = GetReconstructedImage<cmtk::Interpolators::Cubic>( volume, refImage, xformsToPassImages );
      break;
    case cmtk::Interpolators::HAMMING_SINC: 
      correctedVolume = GetReconstructedImage< cmtk::Interpolators::HammingSinc<3> >( volume, refImage, xformsToPassImages );
      break;
    case cmtk::Interpolators::COSINE_SINC: 
      correctedVolume = GetReconstructedImage< cmtk::Interpolators::CosineSinc<3> >( volume, refImage, xformsToPassImages );
      break;
    }

  cmtk::UniformVolume::SmartPtr outputImage = correctedVolume;
  if ( !WriteImagesAsFloat && outputImage->GetData()->GetType() != volume->GetData()->GetType() )
    {
    outputImage = cmtk::UniformVolume::SmartPtr( outputImage->CloneGrid() );
    outputImage->SetData( cmtk::TypedArray::SmartPtr( correctedVolume->GetData()->Convert( volume->GetData()->GetType() ) ) );
    }
  cmtk::VolumeIO::Write( *outputImage, OutputFilePath );
  
  return 0;
}

#include "cmtkSafeMain" 
