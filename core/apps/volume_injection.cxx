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

#include <cmtkVolumeInjectionReconstruction.h>

#include <cmtkXformIO.h>

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

const char* OutputImagePath = "volume_injection.hdr";
bool WriteImagesAsFloat = false;

bool VolumeInjectionIsotropic = false;
double VolumeInjectionSigma = 1;
int VolumeInjectionRadius = 0;

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

int CropFromIndex[3], CropToIndex[3];
bool UseCropRegion = false;

void
CallbackCrop( const char* arg )
{
  if ( 6 == sscanf( arg, "%d,%d,%d:%d,%d,%d", CropFromIndex, CropFromIndex+1, CropFromIndex+2, CropToIndex, CropToIndex+1, CropToIndex+2 ) )
    UseCropRegion = true;
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
  float gridOrigin[3] = { 0, 0, 0 };

  const size_t numArgs = sscanf( arg, "%d,%d,%d:%f,%f,%f:%f,%f,%f", gridDims, gridDims+1, gridDims+2, gridDelta, gridDelta+1, gridDelta+2, gridOrigin, gridOrigin+1, gridOrigin+2 );
  if ( (numArgs != 6) && (numArgs != 9) )
    {
    cmtk::StdErr.printf( "ERROR: reconstruction volume definition must be int,int,int:float,float,float or int,int,int:float,float,float:float,float,float\n", arg );
    exit( 1 );
    }
  
  ReconGrid = cmtk::UniformVolume::SmartPtr( new cmtk::UniformVolume( gridDims, gridDelta[0], gridDelta[1], gridDelta[2] ) );
  
  if ( numArgs == 9 )
    {
    ReconGrid->SetOrigin( cmtk::Vector3D( gridOrigin[0], gridOrigin[1], gridOrigin[2] ) );
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
  cmtk::VolumeIO::Write( outputImage, path, Verbose );
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
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Volume injection" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Reconstruction a high-resolution volume from multiple co-registered (low-resolution) images using forward volume injection" ); 
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] refImage xform0 inImage0 [xform1 inImage1 ...]" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose operation" );

    cl.AddSwitch( Key( 'x', "exclude-first-image" ), &ExcludeFirstImage, true, "Exclude first image from reconstruction as a separate registration target image)" );
    cl.AddCallback( Key( "recon-grid" ), CallbackReconGrid, "Define reconstruction grid as Nx,Ny,Nz:dX,dY,dZ[:Ox,Oy,Oz] (dims:pixel:origin)" );
    cl.AddOption( Key( 'R', "recon-grid-path" ), &ReconstructionGridPath, "Give path to grid that defines reconstructed image grid [including offset]" );
    cl.AddCallback( Key( "crop" ), CallbackCrop, "Crop reference to pixel region x0,y0,z1:x1,y1,z1" );
    cl.AddOption( Key( 'o', "output" ), &OutputImagePath, "Output image path [default: reconstructed.hdr]" );

    cl.AddCallback( Key( 'W', "pass-weight" ), CallbackSetPassWeight, "Set contribution weight for a pass in the form 'pass:weight'" );

    cl.AddSwitch( Key( "isotropic-injection" ), &VolumeInjectionIsotropic, true, "Use isotropic volume injection [default: scaled with pass image pixel size per dimension]" );
    cl.AddOption( Key( 'S', "gauss-sigma" ), &VolumeInjectionSigma, "Gauss contribution [default: 1]" );
    cl.AddOption( Key( 'r', "radius" ), &VolumeInjectionRadius, "VolumeInjectionRadius of affected pixel [default: 0]" );

    cl.AddSwitch( Key( 'F', "write-images-as-float" ), &WriteImagesAsFloat, true, "Write output images as floating point [default: same as input]" );

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
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  if ( ExcludeFirstImage )
    ReconstructionGridPath = ImagePaths[0];

  for ( size_t idx = (ExcludeFirstImage?1:0); idx < ImagePaths.size(); ++idx )
    {
    cmtk::UniformVolume::SmartPtr image( cmtk::VolumeIO::ReadOriented( ImagePaths[idx], Verbose ) );
    if ( ! image || ! image->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read image " << ImagePaths[idx] << "\n";
      exit( 1 );
      }

    cmtk::Xform::SmartPtr xform( new cmtk::AffineXform );
    if ( XformPaths[idx] && strcmp( XformPaths[idx], "--" ) )
      {
      xform = cmtk::Xform::SmartPtr( cmtk::XformIO::Read( XformPaths[idx] ) );
      if ( ! xform )
	{
	cmtk::StdErr << "ERROR: Could read transformation from file" << XformPaths[idx] << "\n";
	}
      }

    Images.push_back( image );
    Xforms.push_back( xform );
    }

  if ( ! ReconGrid )
    ReconGrid = Images[0];
  
  if ( ReconstructionGridPath )
    {
    ReconGrid = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( ReconstructionGridPath ) );
    if ( ! ReconGrid )
      {
      cmtk::StdErr << "ERROR: Could not read reconstruction grid from image " << ReconstructionGridPath << "\n";
      exit( 1 );
      }
    }
  
  if ( UseCropRegion )
    {
    ReconGrid->SetCropRegion( CropFromIndex, CropToIndex );
    ReconGrid = cmtk::UniformVolume::SmartPtr( ReconGrid->GetCroppedVolume() );
    }

  if ( Verbose )
    {
    cmtk::StdErr.printf( "Reconstruction grid: %dx%dx%d pixels, %fx%fx%f pixel size, origin=%f,%f,%f\n",
		      ReconGrid->m_Dims[0], ReconGrid->m_Dims[1], ReconGrid->m_Dims[2], (float)ReconGrid->m_Delta[0], (float)ReconGrid->m_Delta[1], (float)ReconGrid->m_Delta[2],
		      (float)ReconGrid->m_Origin[0], (float)ReconGrid->m_Origin[1], (float)ReconGrid->m_Origin[2] );
    }
  
  cmtk::VolumeInjectionReconstruction injection( ReconGrid, Images );
  injection.SetTransformationsToPassImages( Xforms );
  
  for ( std::map<size_t,float>::const_iterator it = PassWeights.begin(); it != PassWeights.end(); ++it )
    {
    injection.SetPassWeight( it->first, it->second );
    }
  
  if ( VolumeInjectionIsotropic )
    injection.VolumeInjectionIsotropic( VolumeInjectionSigma, VolumeInjectionRadius );
  else
    injection.VolumeInjectionAnisotropic( VolumeInjectionSigma, VolumeInjectionRadius );
  
  if ( OutputImagePath )
    {
    WriteOutputImage( injection.GetCorrectedImage(), OutputImagePath );
    }
  
  return 0;
}

