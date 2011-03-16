/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkThreads.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkVector3D.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkCubicInterpolator.h>
#include <Base/cmtkSincInterpolator.h>
#include <Base/cmtkTransformChangeFromSpaceAffine.h>

#include <Registration/cmtkAffineRegistration.h>
#include <Registration/cmtkProtocolCallback.h>

#include <Recon/cmtkInverseInterpolationVolumeReconstruction.h>
#include <Recon/cmtkPointSpreadFunctionBox.h>
#include <Recon/cmtkPointSpreadFunctionGaussian.h>
#include <Recon/cmtkDeblurringVolumeReconstruction.h>

#include <IO/cmtkXformIO.h>
#include <IO/cmtkVolumeIO.h>

#include <algorithm>
#include <map>
#include <vector>

const char* ReconstructionGridPath = NULL;
bool ExcludeFirstImage = false;

std::vector<const char*> XformPaths;
std::vector<const char*> ImagePaths;

std::vector<cmtk::Xform::SmartPtr> Xforms;
std::vector<cmtk::UniformVolume::SmartPtr> Images;

const char* OutputImagePath = "reconstructed.nii";
const char* LowestMaxErrorImagePath = NULL;

bool VolumeInjectionIsotropic = false;
double InjectionKernelSigma = 1;
double InjectionKernelRadius = 2;

bool FourthOrderError = false;
double ConstraintWeightLNorm = 0;

int InverseInterpolationKernel = cmtk::Interpolators::CUBIC;

enum {
  DEBLURRING_BOX = 1,
  DEBLURRING_GAUSSIAN = 2
};
int DeblurringKernel = 0;

cmtk::Vector3D PointSpreadFunction;
bool PointSpreadFunctionSet = false;
cmtk::Types::Coordinate PointSpreadFunctionScale = 1.0;

void
CallbackSetPSF( const char* arg )
{
  float xyz[3];
  if ( 3 != sscanf( arg, "%f,%f,%f", xyz, xyz+1, xyz+2 ) )
    {
    throw "ERROR: point spread function size must be given as three comma-separated real values: x,y,z\n";
    }
  PointSpreadFunction = cmtk::Vector3D( xyz );
  PointSpreadFunctionSet = true;
}

int NumberOfIterations = 20;
bool RegionalIntensityTruncation = true;

const char* SplattedImagePath = NULL;
bool WriteImagesAsFloat = false;

std::map<size_t,float> PassWeights;

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

bool UseCropRegion = false;
cmtk::DataGrid::RegionType CropRegion;

void
CallbackCrop( const char* arg )
{
  int cropFrom[3], cropTo[3];
  UseCropRegion = (6 == sscanf( arg, "%d,%d,%d,%d,%d,%d", cropFrom, cropFrom+1, cropFrom+2, cropTo,cropTo+1,cropTo+2 ) );

  if ( UseCropRegion )
    {
    CropRegion = cmtk::DataGrid::RegionType( cmtk::DataGrid::IndexType( cropFrom ), cmtk::DataGrid::IndexType( cropTo ) );
    }
  else
    {
    cmtk::StdErr.printf( "ERROR: string '%s' does not describe a valid crop region\n", arg );
    throw cmtk::ExitException( 1 );
    }
}

cmtk::UniformVolume::SmartPtr ReconGrid( NULL );

void
CallbackReconGrid( const char* arg )
{
  int gridDims[3] = { 0, 0, 0 };
  float gridDelta[3] = { 0, 0, 0 };
  float gridOffset[3] = { 0, 0, 0 };

  const size_t numArgs = sscanf( arg, "%d,%d,%d:%f,%f,%f:%f,%f,%f", gridDims, gridDims+1, gridDims+2, gridDelta, gridDelta+1, gridDelta+2, gridOffset, gridOffset+1, gridOffset+2 );
  if ( (numArgs != 6) && (numArgs != 9) )
    {
    cmtk::StdErr.printf( "ERROR: reconstruction volume definition must be int,int,int:float,float,float or int,int,int:float,float,float:float,float,float\n", arg );
    throw cmtk::ExitException( 1 );
    }
  
  ReconGrid = cmtk::UniformVolume::SmartPtr( new cmtk::UniformVolume( cmtk::UniformVolume::IndexType( gridDims ), gridDelta[0], gridDelta[1], gridDelta[2] ) );
  ReconGrid->m_MetaInformation[cmtk::META_SPACE] = ReconGrid->m_MetaInformation[cmtk::META_SPACE_ORIGINAL] = cmtk::AnatomicalOrientation::ORIENTATION_STANDARD;

  if ( numArgs == 9 )
    {
    ReconGrid->SetOffset( cmtk::Vector3D( gridOffset ) );
    }
}

void
WriteOutputImage( cmtk::UniformVolume::SmartPtr& image, const char* path )
{
  cmtk::UniformVolume::SmartPtr outputImage = image;

  const cmtk::ScalarDataType type = Images[0]->GetData()->GetType();
  if ( !WriteImagesAsFloat && (outputImage->GetData()->GetType() != type) )
    {
    outputImage = cmtk::UniformVolume::SmartPtr( outputImage->CloneGrid() );
    outputImage->SetData( cmtk::TypedArray::SmartPtr( image->GetData()->Convert( type ) ) );
    }
  cmtk::VolumeIO::Write( *outputImage, path );
}

template<class TRecon>
cmtk::UniformVolume::SmartPtr
ReconstructVolume()
{
  TRecon volRecon( ReconGrid, Images );
  volRecon.SetTransformationsToPassImages( Xforms );

  for ( std::map<size_t,float>::const_iterator it = PassWeights.begin(); it != PassWeights.end(); ++it )
    {
    volRecon.SetPassWeight( it->first, it->second );
    }
  
  volRecon.SetUseRegionalIntensityTruncation( RegionalIntensityTruncation );
  volRecon.SetUseFourthOrderError( FourthOrderError );
  volRecon.SetConstraintWeightLNorm( ConstraintWeightLNorm );

  cmtk::DebugOutput( 1 ) << "Volume injection...\n";
  if ( VolumeInjectionIsotropic )
    volRecon.VolumeInjectionIsotropic( InjectionKernelSigma, InjectionKernelRadius );
  else
    volRecon.VolumeInjectionAnisotropic( InjectionKernelSigma, InjectionKernelRadius );

  if ( SplattedImagePath )
    {
    WriteOutputImage( volRecon.GetCorrectedImage(), SplattedImagePath );
    }

  const double timeBaseline = cmtk::Timers::GetTimeProcess();
  
  if ( NumberOfIterations )
    {
    volRecon.Optimize( NumberOfIterations );
    }

  cmtk::DebugOutput( 1 ) << "OPT_TIME\t" << cmtk::Timers::GetTimeProcess() - timeBaseline << "\n";

  if ( LowestMaxErrorImagePath )
    {
    WriteOutputImage( volRecon.GetLowestMaxErrorImage(), LowestMaxErrorImagePath );
    }
  
  return volRecon.GetCorrectedImage();
}

template<class TRecon>
cmtk::UniformVolume::SmartPtr
ReconstructVolumeDeblurring()
{
  cmtk::Vector3D psf = PointSpreadFunctionScale * PointSpreadFunction;
  TRecon volRecon( ReconGrid, Images, psf );
  volRecon.SetTransformationsToPassImages( Xforms );

  for ( std::map<size_t,float>::const_iterator it = PassWeights.begin(); it != PassWeights.end(); ++it )
    {
    volRecon.SetPassWeight( it->first, it->second );
    }
  
  volRecon.SetUseRegionalIntensityTruncation( RegionalIntensityTruncation );
  volRecon.SetUseFourthOrderError( FourthOrderError );
  volRecon.SetConstraintWeightLNorm( ConstraintWeightLNorm );

  cmtk::DebugOutput( 1 ) << "Volume injection...\n";
  if ( VolumeInjectionIsotropic )
    volRecon.VolumeInjectionIsotropic( InjectionKernelSigma, InjectionKernelRadius );
  else
    volRecon.VolumeInjectionAnisotropic( InjectionKernelSigma, InjectionKernelRadius );

  if ( SplattedImagePath )
    {
    WriteOutputImage( volRecon.GetCorrectedImage(), SplattedImagePath );
    }

  const double timeBaseline = cmtk::Timers::GetTimeProcess();
  
  if ( NumberOfIterations )
    {
    volRecon.Optimize( NumberOfIterations );
    }

  cmtk::DebugOutput( 1 ) << "OPT_TIME\t" << cmtk::Timers::GetTimeProcess() - timeBaseline << "\n";

  if ( LowestMaxErrorImagePath )
    {
    WriteOutputImage( volRecon.GetLowestMaxErrorImage(), LowestMaxErrorImagePath );
    }
  
  return volRecon.GetCorrectedImage();
}

int
doMain( const int argc, const char* argv[] )
{
  cmtk::Threads::CheckEnvironment(); // need this to check for "CMTK_NUM_THREADS" and constrain OpenMP accordingly

  /*
  // Parse command line
  */
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Volume reconstruction" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Iterative volume reconstruction from co-registered images using inverse interpolation or joint deblurring" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] refImage xform0 inImage0 [xform1 inImage1 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Input", "Input Options" );
    cl.AddSwitch( Key( 'x', "exclude-first-image" ), &ExcludeFirstImage, true, "Exclude first image from reconstruction as a separate registration target image)" );
    cl.AddCallback( Key( "crop" ), CallbackCrop, "Crop reference to pixel region x0,y0,z1:x1,y1,z1" );
    cl.AddCallback( Key( 'W', "pass-weight" ), CallbackSetPassWeight, "Set contribution weight for a pass in the form 'pass:weight'" );
    cl.EndGroup();

    cl.BeginGroup( "ReconGrid", "Reconstruction Grid" );
    cl.AddCallback( Key( "recon-grid" ), CallbackReconGrid, "Define reconstruction grid as Nx,Ny,Nz:dX,dY,dZ[:Ox,Oy,Oz] (dims:pixel:offset)" );
    cl.AddOption( Key( 'R', "recon-grid-path" ), &ReconstructionGridPath, "Give path to grid that defines reconstructed image grid [including offset]" );
    cl.EndGroup();

    cl.BeginGroup( "Injection", "Initial Volume Injection Parameters" );
    cl.AddSwitch( Key( "isotropic-injection" ), &VolumeInjectionIsotropic, true, "Use isotropic volume injection [otherwise: scaled with pass image pixel size per dimension]" );
    cl.AddOption( Key( 'S', "injection-kernel-sigma" ), &InjectionKernelSigma, "Standard deviation of Gaussian kernel for volume injection in multiples of pixel size in each direction." );
    cl.AddOption( Key( 'r', "injection-kernel-radius" ), &InjectionKernelRadius, "Truncation radius factor of injection kernel. The kernel is truncated at sigma*radius, where sigma is the kernel standard deviation." );
    cl.EndGroup();
    
    cl.BeginGroup( "Reconstruction", "Volume Reconstruction Options" );
    cmtk::CommandLine::EnumGroup<int>::SmartPtr kernelGroup = 
      cl.AddEnum( "inverse-interpolation-kernel", &InverseInterpolationKernel, "Kernel for the inverse interpolation reconstruction" );
    kernelGroup->AddSwitch( Key( 'C', "cubic" ), cmtk::Interpolators::CUBIC, "Tricubic interpolation" );
    kernelGroup->AddSwitch( Key( 'L', "linear" ), cmtk::Interpolators::LINEAR, "Trilinear interpolation (faster but less accurate)" );
    kernelGroup->AddSwitch( Key( 'H', "hamming-sinc" ), cmtk::Interpolators::HAMMING_SINC, "Hamming-windowed sinc interpolation" );
    kernelGroup->AddSwitch( Key( 'O', "cosine-sinc" ), cmtk::Interpolators::COSINE_SINC, "Cosine-windowed sinc interpolation (most accurate but slowest)" );

    cmtk::CommandLine::EnumGroup<int>::SmartPtr deblurGroup = 
      cl.AddEnum( "deblurring", &DeblurringKernel, "Kernel shape to approximate the point spread function for joint deblurring reconstruction (selecting one of these disables inverse interpolation reconstruction)" );
    deblurGroup->AddSwitch( Key( "box" ), (int)DEBLURRING_BOX, "Box-shaped kernel" );
    deblurGroup->AddSwitch( Key( "gaussian" ), (int)DEBLURRING_GAUSSIAN, "Gaussian kernel" );

    cl.AddCallback( Key( "psf" ), CallbackSetPSF, "Explicitly set point spread function size as x,y,z. Use with 'deblurring' kernel reconstrunction." );
    cl.AddOption( Key( "psf-scale" ), &PointSpreadFunctionScale, "Scale point spread function size by this value. Use with 'deblurring' kernel reconstrunction." );
    cl.EndGroup();
    
    cl.BeginGroup( "Optimization", "Optimization Parameters" );
    cl.AddOption( Key( 'n', "num-iterations" ), &NumberOfIterations, "Maximum number of inverse interpolation iterations" );
    cl.AddSwitch( Key( 'f', "fourth-order-error" ), &FourthOrderError, true, "Use fourth-order (rather than second-order) error for optimization." );
    cl.EndGroup();

    cl.BeginGroup( "Regularization", "Regularization Parameters" );
    cl.AddOption( Key( "l-norm-weight" ), &ConstraintWeightLNorm, "Set constraint weight for Tikhonov-type L-Norm regularization (0 disables constraint)" );
    cl.AddSwitch( Key( 'T', "no-truncation" ), &RegionalIntensityTruncation, false, "Turn off non-linear regional intensity truncation" );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'o', "output" ), &OutputImagePath, "Output path for final reconstructed image" );
    cl.AddOption( Key( "write-injected-image" ), &SplattedImagePath, "Write initial volume-injected image to this path" );
    cl.AddOption( Key( "write-lowest-max-error-image" ), &LowestMaxErrorImagePath, "Optional path to write reconstructed image with lowest MAXIMUM error." );
    cl.AddSwitch( Key( 'F', "write-images-as-float" ), &WriteImagesAsFloat, true, "Write output images as floating point" );
    cl.EndGroup();

    cl.Parse( argc, argv );

    ImagePaths.push_back( cl.GetNext() );
    XformPaths.push_back( NULL );
    
    const char* nextXform = cl.GetNext();
    const char* nextImage = cl.GetNext();
    while ( nextXform && nextImage )
      {
      XformPaths.push_back( nextXform );
      ImagePaths.push_back( nextImage );

      nextXform = cl.GetNextOptional();
      nextImage = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  if ( ExcludeFirstImage )
    ReconstructionGridPath = ImagePaths[0];

  if ( ReconstructionGridPath )
    {
    ReconGrid = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( ReconstructionGridPath ) );
    if ( ! ReconGrid )
      {
      cmtk::StdErr << "ERROR: Could not read reconstruction grid from image " << ReconstructionGridPath << "\n";
      throw cmtk::ExitException( 1 );
      }
    }
  
  for ( size_t idx = (ExcludeFirstImage?1:0); idx < ImagePaths.size(); ++idx )
    {
    cmtk::UniformVolume::SmartPtr image( cmtk::VolumeIO::ReadOriented( ImagePaths[idx] ) );
    if ( ! image || ! image->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read image " << ImagePaths[idx] << "\n";
      throw cmtk::ExitException( 1 );
      }

    cmtk::AffineXform::SmartPtr affineXform( new cmtk::AffineXform );
    if ( XformPaths[idx] && strcmp( XformPaths[idx], "--" ) )
      {
      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( XformPaths[idx] ) );
      if ( ! xform )
	{
	cmtk::StdErr << "ERROR: Could read affine transformation from file" << XformPaths[idx] << "\n";
	}
      affineXform = cmtk::AffineXform::SmartPtr::DynamicCastFrom( xform );
      if ( ! affineXform )
	{
	cmtk::StdErr << "ERROR: transformation " << XformPaths[idx] << " is not affine\n";
	}
      }

    if ( affineXform->m_MetaInformation[cmtk::META_SPACE] != cmtk::AnatomicalOrientation::ORIENTATION_STANDARD )
      {
      cmtk::TransformChangeFromSpaceAffine toStandardSpace( *affineXform, *ReconGrid, *image );
      *affineXform = toStandardSpace.GetTransformation();
      affineXform->m_MetaInformation[cmtk::META_SPACE] = cmtk::AnatomicalOrientation::ORIENTATION_STANDARD;
      }

    Images.push_back( image );
    Xforms.push_back( affineXform );
    }

  if ( ! ReconGrid )
    {
    // No recon grid from command line: use first input image.
    ReconGrid = Images[0];
    }
  else
    {
    // If we have a pre-defined reconstruction grid, make its physical coordinates match the first input image
    // First, get the recon grid offset
    const cmtk::UniformVolume::CoordinateVectorType offset = ReconGrid->m_Offset;
    // Convert offset to input image index coordinates
    const cmtk::UniformVolume::CoordinateVectorType indexOffset = cmtk::ComponentDivide( offset, Images[0]->m_Delta );
    // New offset is the index grid offset transformed to physical space
    const cmtk::UniformVolume::CoordinateVectorType newOffset = Images[0]->IndexToPhysical( indexOffset );
    // Copy image-to-physical matrix from input to recon image
    ReconGrid->SetImageToPhysicalMatrix( Images[0]->GetImageToPhysicalMatrix() );
    // Finally, copy new offset into recon image-to-physical matrix.
    for ( int i = 0; i < 3; ++i )
      {
      ReconGrid->m_IndexToPhysicalMatrix[3][i] = newOffset[i];
      }
    }
  
  if ( UseCropRegion )
    {
    ReconGrid->CropRegion() = CropRegion;
    ReconGrid = cmtk::UniformVolume::SmartPtr( ReconGrid->GetCroppedVolume() );
    }

  cmtk::DebugOutput( 1 ).GetStream().printf( "Reconstruction grid: %dx%dx%d pixels, %fx%fx%f pixel size, offset=%f,%f,%f\n",
					     ReconGrid->m_Dims[0], ReconGrid->m_Dims[1], ReconGrid->m_Dims[2], (float)ReconGrid->m_Delta[0], (float)ReconGrid->m_Delta[1], (float)ReconGrid->m_Delta[2],
					     (float)ReconGrid->m_Offset[0], (float)ReconGrid->m_Offset[1], (float)ReconGrid->m_Offset[2] );
  
  cmtk::UniformVolume::SmartPtr correctedVolume;
  if ( !DeblurringKernel )
    {
    switch ( InverseInterpolationKernel )
      {
      case cmtk::Interpolators::LINEAR:
      default:
	correctedVolume = ReconstructVolume< cmtk::InverseInterpolationVolumeReconstruction<cmtk::Interpolators::Linear> >();
	break;
      case cmtk::Interpolators::CUBIC:
	correctedVolume = ReconstructVolume< cmtk::InverseInterpolationVolumeReconstruction<cmtk::Interpolators::Cubic> >();
	break;
      case cmtk::Interpolators::HAMMING_SINC:
	correctedVolume = ReconstructVolume< cmtk::InverseInterpolationVolumeReconstruction<cmtk::Interpolators::HammingSinc<3> > >();
	break;
      case cmtk::Interpolators::COSINE_SINC:
	correctedVolume = ReconstructVolume< cmtk::InverseInterpolationVolumeReconstruction<cmtk::Interpolators::CosineSinc<3> > >();
	break;
      }
    }
  else
    {
    if ( ! PointSpreadFunctionSet )
      {
      cmtk::StdErr << "ERROR: must set point spread function size for deblurring reconstruction\n";
      throw cmtk::ExitException( 1 );
      }

    switch ( DeblurringKernel )
      {
      case DEBLURRING_BOX:
      default:
	correctedVolume = ReconstructVolumeDeblurring< cmtk::DeblurringVolumeReconstruction<cmtk::PointSpreadFunctions::Box> >();
	break;
      case DEBLURRING_GAUSSIAN:
	correctedVolume = ReconstructVolumeDeblurring< cmtk::DeblurringVolumeReconstruction<cmtk::PointSpreadFunctions::Gaussian> >();
	break;
      }
    }
  
  WriteOutputImage( correctedVolume, OutputImagePath );
  
  return 0;
}

#include "cmtkSafeMain"
