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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkconfig.h>

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>
#include <cmtkProgress.h>

#include <cmtkXform.h>
#include <cmtkUniformVolume.h>
#include <cmtkTypedArray.h>

#include <cmtkXformIO.h>
#include <cmtkVolumeIO.h>

#include <cmtkUniformVolumeInterpolatorPartialVolume.h>
#include <cmtkSincInterpolator.h>
#include <cmtkLinearInterpolator.h>
#include <cmtkCubicInterpolator.h>
#include <cmtkNearestNeighborInterpolator.h>

#include <cmtkReformatVolume.h>
#include <cmtkReformatVolumeJacobian.cxx>
#include <cmtkReformatVolumePlain.cxx>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace reformatx
{
#endif
bool Verbose = false;

bool PushMode = false;

bool AutoScaleReference = false;
bool AutoScaleFloating = false;

bool JacobianCorrectGlobal = true;

bool TargetMask = false;

const char* TargetVolumeName = NULL;
const char* ReferenceVolumeName = NULL;
const char* FloatingVolumeName = NULL;

cmtk::UniformVolume::SmartPtr UserDefinedTargetVolume;

const char*
CallbackTargetVolume( const char* arg )
{
  int gridDims[3] = { 0, 0, 0 };
  float gridDelta[3] = { 0, 0, 0 };
  float gridOrigin[3] = { 0, 0, 0 };

  const size_t numArgs = 
    sscanf( arg, "%d,%d,%d:%f,%f,%f:%f,%f,%f", gridDims, gridDims+1, gridDims+2, gridDelta, gridDelta+1, gridDelta+2, gridOrigin, gridOrigin+1, gridOrigin+2 );
  if ( (numArgs != 6) && (numArgs != 9) )
    {
    cmtk::StdErr.printf( "ERROR: target volume definition must be int,int,int:float,float,float or int,int,int:float,float,float:float,float,float\n", arg );
    exit( 1 );
    }
  
  UserDefinedTargetVolume = cmtk::UniformVolume::SmartPtr( new cmtk::UniformVolume( gridDims, gridDelta[0], gridDelta[1], gridDelta[2] ) );

  if ( numArgs == 9 )
    {
    UserDefinedTargetVolume->SetOrigin( cmtk::Vector3D( gridOrigin[0], gridOrigin[1], gridOrigin[2] ) );
    }
  return NULL;
}


cmtk::XformList TargetToReference;
cmtk::XformList ReferenceToFloating;

cmtk::ReformatVolume::Mode Mode = cmtk::ReformatVolume::REFORMAT_PLAIN;
cmtk::Interpolators::InterpolationEnum Interpolation = cmtk::Interpolators::LINEAR;
size_t InterpolatorWindowRadius = 3;

// default to auto-selection of data type
cmtk::ScalarDataType DataType = cmtk::TYPE_NONE;

cmtk::Types::DataItem OutPaddingValue = 0;
bool OutPaddingValueFlag = true;

const char* OutputImageName = "reformat.hdr";

bool CropImages = false;
int CropImagesRegionFrom[3] = { 0,0,0 };
int CropImagesRegionTo[3] = { 0,0,0 };

bool TargetImageOffsetReal = false;
bool TargetImageOffsetPixels = false;
cmtk::Vector3D TargetImageOffset( 0, 0, 0 );

const char*
CallbackCropImages( const char* arg )
{
  CropImages = 
    (6 == sscanf( arg, "%d,%d,%d,%d,%d,%d",
		  &CropImagesRegionFrom[0], &CropImagesRegionFrom[1], &CropImagesRegionFrom[2],
		  &CropImagesRegionTo[0], &CropImagesRegionTo[1], &CropImagesRegionTo[2] ) );
  return NULL;
}

const char*
CallbackTargetImageOffset( const char* arg )
{
  float x, y, z;
  if ( 3 != sscanf( arg, "%f,%f,%f", &x, &y, &z ) )
    {
    x = y = z = 0;
    }

  TargetImageOffset[0] = x;
  TargetImageOffset[1] = y;
  TargetImageOffset[2] = z;

  TargetImageOffsetReal = true;
  TargetImageOffsetPixels = false;

  return NULL;
}

const char*
CallbackTargetImageOffsetPixels( const char* arg )
{
  float x, y, z;
  if ( 3 != sscanf( arg, "%f,%f,%f", &x, &y, &z ) )
    {
    x = y = z = 0;
    }

  TargetImageOffset[0] = x;
  TargetImageOffset[1] = y;
  TargetImageOffset[2] = z;

  TargetImageOffsetReal = false;
  TargetImageOffsetPixels = true;

  return NULL;
}

void
ReformatPushforward()
{
  cmtk::UniformVolume::SmartPtr targetVolume = UserDefinedTargetVolume;
  if ( TargetVolumeName )
    {
    targetVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( TargetVolumeName, Verbose ) );
    }
  if ( ! targetVolume ) 
    {
    cmtk::StdErr << "ERROR: a target volume must be given\n";
    exit( 1 );
    }
  if ( CropImages )
    {
    targetVolume->SetCropRegion( CropImagesRegionFrom, CropImagesRegionTo );
    targetVolume = cmtk::UniformVolume::SmartPtr( targetVolume->GetCroppedVolume() );
    }
  
  TargetToReference.SetEpsilon( 0.1 * targetVolume->GetMinDelta() );

  cmtk::UniformVolume::SmartPtr floatingVolume;
  if ( FloatingVolumeName )
    {
    floatingVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( FloatingVolumeName, Verbose ) );
    if ( ! floatingVolume ) 
      {
      cmtk::StdErr << "ERROR: floating volume " << FloatingVolumeName << " could not be read\n";
      exit( 1 );
      }
    if ( ! floatingVolume->GetData() ) 
      {
      cmtk::StdErr << "ERROR: floating volume data must exist in PushForward mode\n";
      exit( 1 );
      }    
    }
  
  cmtk::ConsoleProgress progressIndicator;

  cmtk::TypedArray::SmartPtr reformatData;
  reformatData = cmtk::TypedArray::SmartPtr( cmtk::ReformatVolume::ReformatPushForwardAccumulate( floatingVolume, TargetToReference, targetVolume ) );

  if ( reformatData ) 
    {
    targetVolume->SetData( reformatData );
    cmtk::VolumeIO::Write( targetVolume, OutputImageName, Verbose );
    } 
  else 
    {
    cmtk::StdErr << "ERROR: no reformatted data was generated.";
    }
}

template<class TInterpolator>
void InitializeReformatVolume( cmtk::TypedArray::SmartPtr& reformatData, cmtk::UniformVolume::SmartPtr& targetVolume, cmtk::UniformVolume::SmartPtr& referenceVolume, cmtk::UniformVolume::SmartPtr& floatingVolume )
{
  switch ( Mode ) 
    {
    default:
    case cmtk::ReformatVolume::REFORMAT_PLAIN: 
    {
    if ( ! floatingVolume )
      {
      cmtk::StdErr << "ERROR: floating volume must be defined in plain reformat mode; use '--floating' command line option\n";
      exit( 1 );
      }

    cmtk::ReformatVolume::Plain plain( DataType );
    typename TInterpolator::SmartPtr interpolator ( new TInterpolator (floatingVolume) );
    if ( OutPaddingValueFlag )
      plain.SetPaddingValue( OutPaddingValue );
    reformatData = cmtk::TypedArray::SmartPtr( cmtk::ReformatVolume::Reformat( targetVolume, TargetToReference, referenceVolume, ReferenceToFloating, interpolator, plain ) );
    break;
    }
    case cmtk::ReformatVolume::REFORMAT_JACOBIAN: 
    {
    cmtk::ReformatVolume::Jacobian jacobian( DataType, JacobianCorrectGlobal );
    typename TInterpolator::SmartPtr interpolator ( new TInterpolator (floatingVolume) );
    reformatData = cmtk::TypedArray::SmartPtr( cmtk::ReformatVolume::Reformat( targetVolume, TargetToReference, referenceVolume, ReferenceToFloating, interpolator, jacobian ) );
    break;
    }
    }
}

void
ReformatPullback()
{
  cmtk::UniformVolume::SmartPtr targetVolume = UserDefinedTargetVolume;
  if ( TargetVolumeName )
    {
    if ( TargetMask )
      targetVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( TargetVolumeName, Verbose ) );
    else
      targetVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( TargetVolumeName, Verbose ) );
    }

  if ( ! targetVolume ) 
    {
    cmtk::StdErr << "ERROR: could not read target volume " << TargetVolumeName << "\n";
    exit( 1 );
    }

  if ( CropImages )
    {
    targetVolume->SetCropRegion( CropImagesRegionFrom, CropImagesRegionTo );
    targetVolume = cmtk::UniformVolume::SmartPtr( targetVolume->GetCroppedVolume() );
    }
  
  TargetToReference.SetEpsilon( 0.1 * targetVolume->GetMinDelta() );
  ReferenceToFloating.SetEpsilon( 0.1 * targetVolume->GetMinDelta() );

  if ( TargetImageOffsetPixels )
    {
    for ( int dim = 0; dim < 3; ++dim )
      TargetImageOffset[dim] *= targetVolume->Delta[dim];
    targetVolume->SetOrigin( TargetImageOffset );

    if ( Verbose )
      {
      cmtk::StdErr << "INFO: setting reference image offset to " << TargetImageOffset[0] << "/" << TargetImageOffset[1] << "/" << TargetImageOffset[2] << "\n";
      }
    }
  
  if ( TargetImageOffsetReal )
    {
    targetVolume->SetOrigin( TargetImageOffset );
    
    if ( Verbose )
      {
      cmtk::StdErr << "INFO: setting reference image offset to " << TargetImageOffset[0] << "/" << TargetImageOffset[1] << "/" << TargetImageOffset[2] << "\n";
      }
    }
  
  cmtk::UniformVolume::SmartPtr referenceVolume;
  if ( ReferenceVolumeName ) 
    {
    referenceVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( ReferenceVolumeName, Verbose ) );
    if ( ! referenceVolume ) 
      {
      cmtk::StdErr << "ERROR: reference volume " << ReferenceVolumeName << " could not be read\n";
      exit( 1 );
      }
    }

  cmtk::UniformVolume::SmartPtr floatingVolume;
  if ( FloatingVolumeName )
    {
    floatingVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( FloatingVolumeName, Verbose ) );
    if ( ! floatingVolume ) 
      {
      cmtk::StdErr << "ERROR: floating volume " << FloatingVolumeName << " could not be read\n";
      exit( 1 );
      }
    }
  
  cmtk::TypedArray::SmartPtr floatingData( NULL );
  if ( floatingVolume ) 
    {
    floatingData = floatingVolume->GetData();
    if ( ! floatingData ) 
      {
      cmtk::StdErr << "ERROR: floating volume " << FloatingVolumeName << " seems to have no data\n";
      exit( 1 );
      }
    }
  
  if ( AutoScaleReference || AutoScaleFloating ) 
    {    
    // in cases with no reference volume, use target for intensity reference
    cmtk::TypedArray::SmartPtr referenceData = (referenceVolume) ? referenceVolume->GetData() : targetVolume->GetData();
    if ( ! referenceData ) 
      {
      cmtk::StdErr << "ERROR: neither reference nor target volume seem to have pixel data\n";
      exit( 1 );
      }
    
    cmtk::Types::DataItem mean1, var1, mean2, var2;
    referenceData->GetStatistics( mean1, var1 );
    floatingData->GetStatistics( mean2, var2 );
    
    if ( Verbose ) 
      {
      cmtk::StdErr.printf( "Before auto-rescale: [1] %f +/- %f, [2] %f +/- %f\n", mean1, sqrt( var1 ), mean2, sqrt( var2 ) );
      }
    
    // auto rescale, that is, determine scaling factor and offset so that after
    // scaling, the intensities in both images have the same mean and standard
    // deviation. Note that due to the squaring of values in std.dev.
    // computation, the spread will not be exactly identical. For integer data,
    // there is also the inherently limited precision.
    if ( AutoScaleReference ) 
      {
      const cmtk::Types::DataItem factor = sqrt(var2) / sqrt(var1);
      referenceData->Rescale( factor, mean2 - factor * mean1 );
      referenceData->GetStatistics( mean1, var1 );
      } 
    else
      {
      const cmtk::Types::DataItem factor = sqrt(var1) / sqrt(var2);
      floatingData->Rescale( factor, mean1 - factor * mean2 );
      floatingData->GetStatistics( mean2, var2 );
      }
    
    if ( Verbose ) 
      {
      cmtk::StdErr.printf( "After auto-rescale: [1] %f +/- %f, [2] %f +/- %f\n", mean1, sqrt( var1 ), mean2, sqrt( var2 ) );
      }
    }
  
  if ( !TargetMask ) 
    {
    // unless we're in mask mode, remove target pixel data.
    cmtk::TypedArray::SmartPtr nullData( NULL );
    targetVolume->SetData( nullData );
    } 
  else
    {
    if ( Verbose )
      cmtk::StdErr << "INFO: Using target data as binary mask.\n";
    }
  
  cmtk::ConsoleProgress progressIndicator;

  cmtk::TypedArray::SmartPtr reformatData;
  switch ( Interpolation ) 
    {
    default:
    case cmtk::Interpolators::LINEAR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> TInterpolator;
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
    break;
    }
    case cmtk::Interpolators::CUBIC:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Cubic> TInterpolator;
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
    break;
    }
    case cmtk::Interpolators::NEAREST_NEIGHBOR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor> TInterpolator;
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
    break;
    }
    case cmtk::Interpolators::COSINE_SINC:
    {
    switch ( InterpolatorWindowRadius )
      {
      case 2:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<2> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }
      case 3:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<3> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }
      case 4:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<4> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }
      case 5:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<5> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }      
      default:
	cmtk::StdErr.printf( "ERROR: Sinc window radius %d is not supported.\n", (int)InterpolatorWindowRadius );
      }
    }
    break;
    case cmtk::Interpolators::HAMMING_SINC:
    {
    switch ( InterpolatorWindowRadius )
      {
      case 2:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<2> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }
      case 3:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<3> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }
      case 4:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<4> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }
      case 5:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<5> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
      break;
      }
      default:
	cmtk::StdErr.printf( "ERROR: Sinc window radius %d is not supported.\n", (int)InterpolatorWindowRadius );
      }
    }
    break;
    case cmtk::Interpolators::PARTIALVOLUME:
    {
    typedef cmtk::UniformVolumeInterpolatorPartialVolume TInterpolator;
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, referenceVolume, floatingVolume );
    break;
    }
    }
  
  if ( reformatData ) 
    {
    targetVolume->SetData( reformatData );
    cmtk::VolumeIO::Write( targetVolume, OutputImageName, Verbose );
    } 
  else 
    {
    cmtk::StdErr << "ERROR: no reformatted data was generated.\n";
    }
}

int 
main( const int argc, char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Volume reformatter" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Extended volume reformatter tool to compute reformatted images and Jacobian maps from arbitrary sequences of concatenated transformations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] [--push] --floating floatingImg target x0 [x1 ...]\n  OR\n"
		       "[options] target [x0 x1 ...] {-j,--jacobian} xx0 [xx1 ...]\n  OR\n"
		       "WHERE x0 ... xN and xx0 ... xxN is [{-i,--inverse}] transformation##" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddCallback( Key( "target-grid" ), CallbackTargetVolume, "Define target grid for reformating as Nx,Ny,Nz:dX,dY,dZ[:Ox,Oy,Oz] (dims:pixel:origin)" );
    cl.AddSwitch( Key( 'm', "mask" ), &TargetMask, true, "Use target pixel data as binary mask." );

    cl.AddSwitch( Key( "auto-scale-reference" ), &AutoScaleReference, true, "Automatically scale reference image intensities to match floating image" );
    cl.AddSwitch( Key( "auto-scale-floating" ), &AutoScaleFloating, true, "Automatically scale floating image intensities to match reference image" );
    
    cl.AddSwitch( Key( "jacobian-correct-global" ), &JacobianCorrectGlobal, true, "Correct Jacobian maps for global scale." );
    cl.AddSwitch( Key( "no-jacobian-correct-global" ), &JacobianCorrectGlobal, false, "Do not correct Jacobian maps for global scale." );

    cl.AddSwitch( Key( 'c', "char" ), &DataType, cmtk::TYPE_CHAR, "8 bits, signed" );
    cl.AddSwitch( Key( 'b', "byte" ), &DataType, cmtk::TYPE_BYTE, "8 bits, unsigned" );
    cl.AddSwitch( Key( 's', "short" ), &DataType, cmtk::TYPE_SHORT, "16 bits, signed" );
    cl.AddSwitch( Key( 'u', "ushort" ), &DataType, cmtk::TYPE_USHORT, "16 bits, unsigned" );
    cl.AddSwitch( Key( 'i', "int" ), &DataType, cmtk::TYPE_INT, "32 bits signed" );
    cl.AddSwitch( Key( 'f', "float" ), &DataType, cmtk::TYPE_FLOAT, "32 bits floating point" );
    cl.AddSwitch( Key( 'd', "double" ), &DataType, cmtk::TYPE_DOUBLE, "64 bits floating point\n" );

    cl.AddOption( Key( 'P', "pad-out" ), &OutPaddingValue, "Padding value for output image", &OutPaddingValueFlag );

    cl.AddCallback( Key( "crop" ), CallbackCropImages, "Crop target image: x0,y0,z0,x1,y1,z2" );
    cl.AddCallback( Key( 'O', "target-offset" ), CallbackTargetImageOffset, "Set target image offset to x,y,z mm" );
    cl.AddCallback( Key( "target-offset-pixels" ), CallbackTargetImageOffsetPixels, "Set target image offset to dx,dy,dz pixels" );
    cl.AddOption( Key( 'F', "floating" ), &FloatingVolumeName, "Format and path of floating image." );
    cl.AddOption( Key( 'o', "outfile" ), &OutputImageName, "Format and path of output image." );

    cl.AddSwitch( Key( "linear" ), &Interpolation, cmtk::Interpolators::LINEAR, "Trilinear interpolation (default)" );
    cl.AddSwitch( Key( "nn" ), &Interpolation, cmtk::Interpolators::NEAREST_NEIGHBOR, "Nearest neighbor interpolation" );
    cl.AddSwitch( Key( "cubic" ), &Interpolation, cmtk::Interpolators::CUBIC, "Tricubic interpolation" );
    cl.AddSwitch( Key( "pv" ), &Interpolation, cmtk::Interpolators::PARTIALVOLUME, "Partial volume interpolation" );
    cl.AddSwitch( Key( "sinc-cosine" ), &Interpolation, cmtk::Interpolators::COSINE_SINC, "Sinc interpolation with cosine window" );
    cl.AddSwitch( Key( "sinc-hamming" ), &Interpolation, cmtk::Interpolators::HAMMING_SINC, "Sinc interpolation with Hamming window" );

    cl.AddOption( Key( "sinc-window-radius" ), &InterpolatorWindowRadius, "Window radius for Sinc interpolation [default: 3]" );

    cl.AddSwitch( Key( "push" ), &PushMode, true, "Push-forward reformating [default: pullback]" );
    
    cl.Parse();

    if ( ! UserDefinedTargetVolume )
      TargetVolumeName = cl.GetNext();
    
    const char* next = cl.GetNextOptional();
    while ( next ) 
      {
      if ( ! strcmp( next, "-j" ) || ! strcmp( next, "--jacobian" ) )
	break;
      
      bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
      if ( inverse ) next = cl.GetNext();
      
      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
      if ( ! xform ) 
	{
	cmtk::StdErr << "ERROR: could not read target-to-reference transformation from " << next << "\n";
	exit( 1 );
	}
      
      TargetToReference.Add( xform, inverse, xform->GetGlobalScaling() );
      next = cl.GetNextOptional();
      }
    
    if ( next ) 
      {
      if ( ! strcmp( next, "-j" ) || ! strcmp( next, "--jacobian" ) ) 
	{
	Mode = cmtk::ReformatVolume::REFORMAT_JACOBIAN;
	next = cl.GetNext();
	} 

      switch ( Mode )
	{
	case cmtk::ReformatVolume::REFORMAT_JACOBIAN:
	  while ( next ) 
	    {
	    bool inverse =
	      ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
	    if ( inverse ) next = cl.GetNext();
	    
	    cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
	    if ( ! xform ) 
	      {
	      cmtk::StdErr << "ERROR: could not read target-to-floating transformation from " << next << "\n";
	      exit( 1 );
	      }
	    
	    ReferenceToFloating.Add( xform, inverse, xform->GetGlobalScaling() );
	    next = cl.GetNextOptional();
	    }
	  break;
	case cmtk::ReformatVolume::REFORMAT_PLAIN:
	default:
	  break;
	}
    }
  }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  if ( PushMode )
    {
    ReformatPushforward();
    }
  else
    {
    ReformatPullback();
    }
  
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace reformatx
} // namespace apps
} // namespace cmtk
#endif

