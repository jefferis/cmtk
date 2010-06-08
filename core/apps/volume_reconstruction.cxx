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

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>
#include <cmtkVector3D.h>

#include <cmtkAffineRegistration.h>
#include <cmtkProtocolCallback.h>

#include <cmtkLinearInterpolator.h>
#include <cmtkCubicInterpolator.h>
#include <cmtkSincInterpolator.h>
#include <cmtkInverseInterpolationVolumeReconstruction.h>

#include <cmtkPointSpreadFunctionBox.h>
#include <cmtkPointSpreadFunctionGaussian.h>
#include <cmtkDeblurringVolumeReconstruction.h>

#include <cmtkXformIO.h>
#include <cmtkTransformChangeFromSpaceAffine.h>

#include <algorithm>
#include <map>
#include <vector>

bool Verbose = false;

const char* ReconstructionGridPath = NULL;
bool ExcludeFirstImage = false;

std::vector<const char*> XformPaths;
std::vector<const char*> ImagePaths;

std::vector<cmtk::Xform::SmartPtr> Xforms;
std::vector<cmtk::UniformVolume::SmartPtr> Images;

const char* OutputImagePath = "reconstructed.hdr";
const char* LowestMaxErrorImagePath = NULL;

bool VolumeInjectionIsotropic = false;
double VolumeInjectionSigma = 1;
double VolumeInjectionRadius = 2;

bool FourthOrderError = false;
double ConstraintWeightLNorm = 0;

cmtk::Interpolators::InterpolationEnum Interpolation = cmtk::Interpolators::LINEAR;

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
    exit( 1 );
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
    exit( 1 );
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
    exit( 1 );
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
  cmtk::VolumeIO::Write( *outputImage, path, Verbose );
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

  if ( Verbose )
    {
    cmtk::StdErr << "Volume injection...\n";
    }
  if ( VolumeInjectionIsotropic )
    volRecon.VolumeInjectionIsotropic( VolumeInjectionSigma, VolumeInjectionRadius );
  else
    volRecon.VolumeInjectionAnisotropic( VolumeInjectionSigma, VolumeInjectionRadius );

  if ( SplattedImagePath )
    {
    WriteOutputImage( volRecon.GetCorrectedImage(), SplattedImagePath );
    }

  const double timeBaseline = cmtk::Timers::GetTimeProcess();
  
  if ( NumberOfIterations )
    {
    volRecon.Optimize( NumberOfIterations );
    }

  if ( Verbose )
    {
    cmtk::StdErr << "OPT_TIME\t" << cmtk::Timers::GetTimeProcess() - timeBaseline << "\n";
    }

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

  if ( Verbose )
    {
    cmtk::StdErr << "Volume injection...\n";
    }
  if ( VolumeInjectionIsotropic )
    volRecon.VolumeInjectionIsotropic( VolumeInjectionSigma, VolumeInjectionRadius );
  else
    volRecon.VolumeInjectionAnisotropic( VolumeInjectionSigma, VolumeInjectionRadius );

  if ( SplattedImagePath )
    {
    WriteOutputImage( volRecon.GetCorrectedImage(), SplattedImagePath );
    }

  const double timeBaseline = cmtk::Timers::GetTimeProcess();
  
  if ( NumberOfIterations )
    {
    volRecon.Optimize( NumberOfIterations );
    }

  if ( Verbose )
    {
    cmtk::StdErr << "OPT_TIME\t" << cmtk::Timers::GetTimeProcess() - timeBaseline << "\n";
    }

  if ( LowestMaxErrorImagePath )
    {
    WriteOutputImage( volRecon.GetLowestMaxErrorImage(), LowestMaxErrorImagePath );
    }
  
  return volRecon.GetCorrectedImage();
}

int
main( const int argc, const char* argv[] )
{
  /*
  // Parse command line
  */
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Volume reconstruction" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Iterative volume reconstruction from co-registered images using inverse interpolation or joint deblurring" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] refImage xform0 inImage0 [xform1 inImage1 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose operation" );

    cl.AddSwitch( Key( 'x', "exclude-first-image" ), &ExcludeFirstImage, true, "Exclude first image from reconstruction as a separate registration target image)" );
    cl.AddCallback( Key( "recon-grid" ), CallbackReconGrid, "Define reconstruction grid as Nx,Ny,Nz:dX,dY,dZ[:Ox,Oy,Oz] (dims:pixel:offset)" );
    cl.AddOption( Key( 'R', "recon-grid-path" ), &ReconstructionGridPath, "Give path to grid that defines reconstructed image grid [including offset]" );
    cl.AddCallback( Key( "crop" ), CallbackCrop, "Crop reference to pixel region x0,y0,z1:x1,y1,z1" );
    cl.AddOption( Key( 'o', "output" ), &OutputImagePath, "Output image path" );

    cl.AddCallback( Key( 'W', "pass-weight" ), CallbackSetPassWeight, "Set contribution weight for a pass in the form 'pass:weight'" );

    cl.AddSwitch( Key( "isotropic-injection" ), &VolumeInjectionIsotropic, true, "Use isotropic volume injection [otherwise: scaled with pass image pixel size per dimension]" );
    cl.AddOption( Key( 'S', "gauss-sigma" ), &VolumeInjectionSigma, "Volume injection Gaussian kernel width factor, taken times pixel size." );
    cl.AddOption( Key( 'r', "radius" ), &VolumeInjectionRadius, "Volume injection kernel truncation factor, taken times pixel size." );
    
    cl.AddSwitch( Key( 'L', "linear" ), &Interpolation, cmtk::Interpolators::LINEAR, "Trilinear interpolation" );
    cl.AddSwitch( Key( 'C', "cubic" ), &Interpolation, cmtk::Interpolators::CUBIC, "Tricubic interpolation" );
    cl.AddSwitch( Key( 'H', "hamming-sinc" ), &Interpolation, cmtk::Interpolators::HAMMING_SINC, "Hamming-windowed sinc interpolation" );
    cl.AddSwitch( Key( 'O', "cosine-sinc" ), &Interpolation, cmtk::Interpolators::COSINE_SINC, "Cosine-windowed sinc interpolation" );

    cl.AddSwitch( Key( "deblurring-box" ), &DeblurringKernel, (int)DEBLURRING_BOX, "Deblurring reconstruction [box PSF]" );
    cl.AddSwitch( Key( "deblurring-gaussian" ), &DeblurringKernel, (int)DEBLURRING_GAUSSIAN, "Deblurring reconstruction [Gaussian PSF]" );
    cl.AddCallback( Key( "psf" ), CallbackSetPSF, "Set point spread function as x,y,z" );
    cl.AddOption( Key( "psf-scale" ), &PointSpreadFunctionScale, "Set point spread function global scale as real value" );

    cl.AddSwitch( Key( 'f', "fourth-order-error" ), &FourthOrderError, true, "Use fourth-order (rather than second-order) error for optimization." );
    cl.AddOption( Key( "l-norm-weight" ), &ConstraintWeightLNorm, "Set constraint weight for L-Norm regularization (values <= 0 disable constraint)" );
    cl.AddOption( Key( 'n', "num-iterations" ), &NumberOfIterations, "Maximum number of inverse interpolation iterations" );
    cl.AddSwitch( Key( 'T', "no-truncation" ), &RegionalIntensityTruncation, false, "Turn off regional intensity truncatrion" );
    
    cl.AddOption( Key( "write-splatted-image" ), &SplattedImagePath, "Write initial Gaussian-splatted image to path" );
    cl.AddOption( Key( "write-lowest-max-error-image" ), &LowestMaxErrorImagePath, "Optional path to write reconstructed image with lowest MAXIMUM error." );
    
    cl.AddSwitch( Key( 'F', "write-images-as-float" ), &WriteImagesAsFloat, true, "Write output images as floating point" );

    cl.Parse();

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
      exit( 1 );
      }
    }
  
  for ( size_t idx = (ExcludeFirstImage?1:0); idx < ImagePaths.size(); ++idx )
    {
    cmtk::UniformVolume::SmartPtr image( cmtk::VolumeIO::ReadOriented( ImagePaths[idx], Verbose ) );
    if ( ! image || ! image->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read image " << ImagePaths[idx] << "\n";
      exit( 1 );
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
    ReconGrid = Images[0];
  
  if ( UseCropRegion )
    {
    ReconGrid->CropRegion() = CropRegion;
    ReconGrid = cmtk::UniformVolume::SmartPtr( ReconGrid->GetCroppedVolume() );
    }

  if ( Verbose )
    {
    cmtk::StdErr.printf( "Reconstruction grid: %dx%dx%d pixels, %fx%fx%f pixel size, offset=%f,%f,%f\n",
		      ReconGrid->m_Dims[0], ReconGrid->m_Dims[1], ReconGrid->m_Dims[2], (float)ReconGrid->m_Delta[0], (float)ReconGrid->m_Delta[1], (float)ReconGrid->m_Delta[2],
		      (float)ReconGrid->m_Offset[0], (float)ReconGrid->m_Offset[1], (float)ReconGrid->m_Offset[2] );
    }
  
  cmtk::UniformVolume::SmartPtr correctedVolume;
  if ( !DeblurringKernel )
    {
    switch ( Interpolation )
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
      exit( 1 );
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


