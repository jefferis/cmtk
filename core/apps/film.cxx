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
#include <cmtkVector3D.h>

#include <cmtkAffineRegistration.h>
#include <cmtkProtocolCallback.h>

#include <cmtkLinearInterpolator.h>
#include <cmtkCubicInterpolator.h>
#include <cmtkSincInterpolator.h>
#include <cmtkInverseInterpolationVolumeReconstruction.h>

#include <cmtkClassStream.h>
#include <cmtkClassStreamXform.h>
#include <cmtkClassStreamAffineXform.h>

#include <algorithm>
#include <map>
#include <vector>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace film
{
#endif
bool Verbose = false;

const char* InputFilePath = NULL;
const char* OutputFilePath = NULL;

int InterleaveAxis = -1;
unsigned int NumberOfPasses = 2;

int RegistrationMetric = 4; // MSD

double InjectionKernelSigma = 0.5;
double InjectionKernelRadius = 2;

bool FourthOrderError = false;
int InverseInterpolationKernel = 1;
int NumberOfIterations = 20;
bool RegionalIntensityTruncation = true;

const char* ReferenceImagePath = NULL;
const char* InjectedImagePath = NULL;

const char* ExportXformsPath = NULL;
const char* ImportXformsPath = NULL;

std::map<size_t,float> PassWeights;

bool WriteImagesAsFloat = false;

const char*
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
    exit( 1 );
    }
  return NULL;
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

  if ( xformsToPassImages.size() == NumberOfPasses )
    {
    volRecon.SetTransformationsToPassImages( xformsToPassImages );
    }
  else
    {
    if ( Verbose )
      {
      cmtk::StdErr << "Computing transformations between passes...\n";
      }
    volRecon.ComputeTransformationsToPassImages( RegistrationMetric );
    xformsToPassImages = volRecon.GetTransformationsToPassImages();
    }
  
  if ( ExportXformsPath )
    {
    cmtk::ClassStream stream( ExportXformsPath, cmtk::ClassStream::WRITE );
    if ( stream.IsValid() )
      {
      if ( Verbose )
	{
	cmtk::StdErr << "Exporting transformations between passes to " << ExportXformsPath << "\n";
	}
      for ( unsigned int pass = 0; pass < NumberOfPasses; ++pass )
	{
	stream << *xformsToPassImages[pass];
	}
      }
    else
      {
      cmtk::StdErr << "ERROR: Could not open transformation file" << ExportXformsPath << "\n";
      }
    }

  if ( Verbose )
    {
    cmtk::StdErr << "Volume injection...\n";
    }
  volRecon.VolumeInjectionAnisotropic( InjectionKernelSigma, InjectionKernelRadius );
  if ( InjectedImagePath )
    {
    cmtk::UniformVolume::SmartPtr outputImage = volRecon.GetCorrectedImage();
    if ( !WriteImagesAsFloat && outputImage->GetData()->GetType() != volume->GetData()->GetType() )
      {
      outputImage = cmtk::UniformVolume::SmartPtr( outputImage->CloneGrid() );
      outputImage->SetData( cmtk::TypedArray::SmartPtr( volRecon.GetCorrectedImage()->GetData()->Convert( volume->GetData()->GetType() ) ) );
      }
    cmtk::VolumeIO::Write( outputImage, InjectedImagePath, Verbose );
    }

  volRecon.Optimize( NumberOfIterations );
  return volRecon.GetCorrectedImage();
}  
  

int
main( int argc, char* argv[] )
{
  /*
  // Parse command line
  */
  try
    {
    cmtk::CommandLine cl( argc, argv );
    typedef cmtk::CommandLine::Key Key;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Fix interleaved motion using inverse interpolation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool splits an interleaved input image into the pass images, co-registers them, and reconstructs a motion-corrected image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] inImage outImage" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Artifact Correction" );

    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose operation" );

    cl.BeginGroup( "interleave", "Interleaving Options" );
    cmtk::CommandLine::EnumGroup::SmartPtr
      interleaveGroup = cl.AddEnum( "interleave-axis", &InterleaveAxis, "Define interleave axis: this is the through-slice direction of the acquisition." );
    interleaveGroup->AddSwitch( Key( "guess-from-input" ), -1, "Guess from input image" );
    interleaveGroup->AddSwitch( Key( 'a', "axial" ), (int)cmtk::AXIS_Z, "Interleaved axial images" );
    interleaveGroup->AddSwitch( Key( 's', "sagittal" ),(int)cmtk::AXIS_X, "Interleaved sagittal images" );
    interleaveGroup->AddSwitch( Key( 'c', "coronal" ), (int)cmtk::AXIS_Y, "Interleaved coronal images" );
    interleaveGroup->AddSwitch( Key( 'x', "interleave-x" ), (int)cmtk::AXIS_X, "Interleaved along x axis" );
    interleaveGroup->AddSwitch( Key( 'y', "interleave-y" ), (int)cmtk::AXIS_Y, "Interleaved along y axis" );
    interleaveGroup->AddSwitch( Key( 'z', "interleave-z" ), (int)cmtk::AXIS_Z, "Interleaved along z axis" );

    cl.AddOption( Key( 'p', "passes" ), &NumberOfPasses, "Number of interleaved passes" );
    cl.AddCallback( Key( 'W', "pass-weight" ), CallbackSetPassWeight, "Set contribution weight for a pass in the form 'pass:weight'" );

    cl.BeginGroup( "motion", "Motion Correction / Registration Options" );
    cl.AddOption( Key( 'R', "reference-image" ), &ReferenceImagePath, "Use a separate high-resolution reference image for registration" );

    cmtk::CommandLine::EnumGroup::SmartPtr
      metricGroup = cl.AddEnum( "registration-metric", &RegistrationMetric, "Registration metric for motion estimation by image-to-image registration." );
    metricGroup->AddSwitch( Key( "nmi" ), 0, "Use Normalized Mutual Information for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "mi" ), 1, "Use standard Mutual Information for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "cr" ), 2, "Use Correlation Ratio for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "msd" ), 4, "Use Mean Squared Differences for pass-to-refereence registration" );
    metricGroup->AddSwitch( Key( "cc" ), 5, "Use Cross-Correlation for pass-to-refereence registration" );

    cl.AddOption( Key( "import-xforms-path" ), &ImportXformsPath, "Path of file from which to import transformations between passes." );
    cl.AddOption( Key( "export-xforms-path" ), &ExportXformsPath, "Path of file to which to export transformations between passes." );

    cl.AddOption( Key( 'S', "injection-kernel-sigma" ), &InjectionKernelSigma, "Standard deviation of Gaussian kernel for volume injection in multiples of pixel size in each direction." );
    cl.AddOption( Key( 'r', "injection-kernel-radius" ), &InjectionKernelRadius, "Truncation radius factor of injection kernel. The kernel is truncated at sigma*radius, where sigma is the kernel standard deviation." );

    cl.BeginGroup( "invint", "Inverse Interpolation Options" );
    cmtk::CommandLine::EnumGroup::SmartPtr kernelGroup = 
      cl.AddEnum( "inverse-interpolation-kernel", &InverseInterpolationKernel, "Kernel for the inverse interpolation reconstruction" );
    kernelGroup->AddSwitch( Key( 'C', "cubic" ), 1, "Tricubic interpolation" );
    kernelGroup->AddSwitch( Key( 'L', "linear" ), 0, "Trilinear interpolation (faster but less accurate)" );
    kernelGroup->AddSwitch( Key( 'H', "hamming-sinc" ), 2, "Hamming-windowed sinc interpolation" );
    kernelGroup->AddSwitch( Key( 'O', "cosine-sinc" ), 3, "Cosine-windowed sinc interpolation (most accurate but slowest)" );

    cl.AddSwitch( Key( 'f', "fourth-order-error" ), &FourthOrderError, true, "Use fourth-order (rather than second-order) error for optimization." );
    cl.AddOption( Key( 'n', "num-iterations" ), &NumberOfIterations, "Maximum number of inverse interpolation iterations" );
    cl.AddSwitch( Key( 'T', "no-truncation" ), &RegionalIntensityTruncation, false, "Turn off regional intensity truncatrion" );
    
    cl.BeginGroup( "output", "Output Options" );
    cl.AddOption( Key( "write-injected-image" ), &InjectedImagePath, "Write initial volume injection image to path" );
    cl.AddSwitch( Key( 'F', "write-images-as-float" ), &WriteImagesAsFloat, true, "Write output images as floating point [default: same as input]" );

    cl.AddParameter( &InputFilePath, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &OutputFilePath, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse();
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  /*
  // Read input image
  */
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InputFilePath, Verbose ) );
  if ( ! volume || ! volume->GetData() )
    {
    cmtk::StdErr << "ERROR: Could not read image " << InputFilePath << "\n";
    return 1;
    }

  cmtk::UniformVolume::SmartPtr refImage;
  if ( ReferenceImagePath )
    {
    refImage = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( ReferenceImagePath, Verbose ) );
    if ( ! refImage || ! refImage->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read image " << ReferenceImagePath << "\n";
      return 1;
      }
    }

  std::vector<cmtk::Xform::SmartPtr> xformsToPassImages;
  if ( ImportXformsPath )
    {
    cmtk::ClassStream stream( ImportXformsPath, cmtk::ClassStream::READ );
    if ( stream.IsValid() )
      {
      if ( Verbose )
	{
	cmtk::StdErr << "Importing transformations between passes from " << ImportXformsPath << "\n";
	}
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
    case 0 : 
      correctedVolume = GetReconstructedImage<cmtk::Interpolators::Linear>( volume, refImage, xformsToPassImages );
      break;
    case 1 : 
      correctedVolume = GetReconstructedImage<cmtk::Interpolators::Cubic>( volume, refImage, xformsToPassImages );
      break;
    case 2 : 
      correctedVolume = GetReconstructedImage< cmtk::Interpolators::HammingSinc<3> >( volume, refImage, xformsToPassImages );
      break;
    case 3 : 
      correctedVolume = GetReconstructedImage< cmtk::Interpolators::CosineSinc<3> >( volume, refImage, xformsToPassImages );
      break;
    }

  cmtk::UniformVolume::SmartPtr outputImage = correctedVolume;
  if ( !WriteImagesAsFloat && outputImage->GetData()->GetType() != volume->GetData()->GetType() )
    {
    outputImage = cmtk::UniformVolume::SmartPtr( outputImage->CloneGrid() );
    outputImage->SetData( cmtk::TypedArray::SmartPtr( correctedVolume->GetData()->Convert( volume->GetData()->GetType() ) ) );
    }
  cmtk::VolumeIO::Write( outputImage, OutputFilePath, Verbose);
  
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace film
} // namespace apps
} // namespace cmtk
#endif

