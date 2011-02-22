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

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkProgressConsole.h>
#include <System/cmtkThreads.h>

#include <Base/cmtkXform.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypedArray.h>

#include <Base/cmtkUniformVolumeInterpolatorPartialVolume.h>
#include <Base/cmtkSincInterpolator.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkCubicInterpolator.h>
#include <Base/cmtkNearestNeighborInterpolator.h>
#include <Base/cmtkAnatomicalOrientation.h>

#include <IO/cmtkXformIO.h>
#include <IO/cmtkVolumeIO.h>

#include <Registration/cmtkReformatVolume.h>

#ifdef CMTK_USE_SQLITE
#  include <Registration/cmtkImageXformDB.h>
#endif

bool Verbose = false;

bool JacobianCorrectGlobal = true;
bool MassPreservingReformat = false;

bool TargetMask = false;

const char* TargetVolumeName = NULL;
const char* FloatingVolumeName = NULL;

cmtk::UniformVolume::SmartPtr UserDefinedTargetVolume;

#ifdef CMTK_USE_SQLITE
const char* updateDB = NULL;
#endif

void
CallbackTargetVolume( const char* arg )
{
  int gridDims[3] = { 0, 0, 0 };
  float gridDelta[3] = { 0, 0, 0 };
  float gridOffset[3] = { 0, 0, 0 };

  const size_t numArgs = sscanf( arg, "%d,%d,%d:%f,%f,%f:%f,%f,%f", gridDims, gridDims+1, gridDims+2, gridDelta, gridDelta+1, gridDelta+2, gridOffset, gridOffset+1, gridOffset+2 );
  if ( (numArgs != 6) && (numArgs != 9) )
    {
    cmtk::StdErr.printf( "ERROR: target volume definition must be int,int,int:float,float,float or int,int,int:float,float,float:float,float,float\n", arg );
    throw cmtk::ExitException( 1 );
    }
  
  UserDefinedTargetVolume = cmtk::UniformVolume::SmartPtr( new cmtk::UniformVolume( cmtk::UniformVolume::IndexType( gridDims ), gridDelta[0], gridDelta[1], gridDelta[2] ) );
  UserDefinedTargetVolume->m_MetaInformation[cmtk::META_SPACE] = 
    UserDefinedTargetVolume->m_MetaInformation[cmtk::META_SPACE_ORIGINAL] = cmtk::AnatomicalOrientation::ORIENTATION_STANDARD;

  if ( numArgs == 9 )
    {
    UserDefinedTargetVolume->SetOffset( cmtk::UniformVolume::CoordinateVectorType( gridOffset ) );
    }
}


cmtk::XformList TargetToReference;
cmtk::XformList ReferenceToFloating;

cmtk::ReformatVolume::Mode Mode = cmtk::ReformatVolume::REFORMAT_PLAIN;
cmtk::Interpolators::InterpolationEnum Interpolation = cmtk::Interpolators::LINEAR;
int InterpolatorWindowRadius = 3;

// default to auto-selection of data type
cmtk::ScalarDataType DataType = cmtk::TYPE_NONE;

cmtk::Types::DataItem OutPaddingValue = 0;
bool OutPaddingValueFlag = true;

const char* OutputImageName = "reformat.nii";

bool CropImages = false;
cmtk::DataGrid::RegionType CropImagesRegion;

void
CallbackCropImages( const char* arg )
{
  int cropFrom[3], cropTo[3];
  CropImages = (6 == sscanf( arg, "%d,%d,%d,%d,%d,%d", cropFrom, cropFrom+1, cropFrom+2, cropTo,cropTo+1,cropTo+2 ) );

  if ( CropImages )
    {
    CropImagesRegion = cmtk::DataGrid::RegionType( cmtk::DataGrid::IndexType( cropFrom ), cmtk::DataGrid::IndexType( cropTo ) );
    }
}

bool TargetImageOffsetReal = false;
bool TargetImageOffsetPixels = false;
cmtk::UniformVolume::CoordinateVectorType TargetImageOffset( cmtk::UniformVolume::CoordinateVectorType::Init( 0 ) );

void
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
}

void
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
}

template<class TInterpolator>
void InitializeReformatVolume( cmtk::TypedArray::SmartPtr& reformatData, cmtk::UniformVolume::SmartPtr& targetVolume, cmtk::UniformVolume::SmartPtr& floatingVolume )
{
  switch ( Mode ) 
    {
    default:
    case cmtk::ReformatVolume::REFORMAT_PLAIN: 
    {
    if ( ! floatingVolume )
      {
      cmtk::StdErr << "ERROR: floating volume must be defined in plain reformat mode; use '--floating' command line option\n";
      throw cmtk::ExitException( 1 );
      }

    cmtk::ReformatVolume::Plain plain( DataType );
    typename TInterpolator::SmartPtr interpolator ( new TInterpolator(*floatingVolume) );
    if ( OutPaddingValueFlag )
      plain.SetPaddingValue( OutPaddingValue );
    reformatData = cmtk::TypedArray::SmartPtr( cmtk::ReformatVolume::ReformatMasked( targetVolume, TargetToReference, ReferenceToFloating, plain, floatingVolume, interpolator ) );

    if ( MassPreservingReformat )
      {
      cmtk::ReformatVolume::Jacobian jacobian( cmtk::TYPE_DOUBLE, false /*correctGlobalScale*/ );
      cmtk::XformList emptyXformList;
      cmtk::TypedArray::SmartPtr jacobianData( cmtk::ReformatVolume::ReformatMasked( targetVolume, emptyXformList, TargetToReference, jacobian, NULL, TInterpolator::SmartConstPtr::Null ) );
      
      const size_t nPixels = reformatData->GetDataSize();
      for ( size_t i = 0; i < nPixels; ++i )
	{
	cmtk::Types::DataItem v, j;
	if ( reformatData->Get( v, i ) && jacobianData->Get( j, i ) )
	  reformatData->Set( v*j, i );
	else
	  reformatData->SetPaddingAt( i );
	}
      }
    break;
    }
    case cmtk::ReformatVolume::REFORMAT_JACOBIAN: 
    {
    cmtk::ReformatVolume::Jacobian jacobian( DataType, JacobianCorrectGlobal );
    if ( OutPaddingValueFlag )
      jacobian.SetPaddingValue( OutPaddingValue );
    reformatData = cmtk::TypedArray::SmartPtr( cmtk::ReformatVolume::ReformatMasked( targetVolume, TargetToReference, ReferenceToFloating, jacobian, NULL, TInterpolator::SmartConstPtr::Null ) );
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
    throw cmtk::ExitException( 1 );
    }

  if ( CropImages )
    {
    targetVolume->CropRegion() = CropImagesRegion;
    targetVolume = cmtk::UniformVolume::SmartPtr( targetVolume->GetCroppedVolume() );
    }
  
  TargetToReference.SetEpsilon( 0.1 * targetVolume->GetMinDelta() );
  ReferenceToFloating.SetEpsilon( 0.1 * targetVolume->GetMinDelta() );

  if ( TargetImageOffsetPixels )
    {
    for ( int dim = 0; dim < 3; ++dim )
      TargetImageOffset[dim] *= targetVolume->m_Delta[dim];
    targetVolume->SetOffset( TargetImageOffset );

    if ( Verbose )
      {
      cmtk::StdOut << "INFO: setting reference image offset to " << TargetImageOffset[0] << "/" << TargetImageOffset[1] << "/" << TargetImageOffset[2] << "\n";
      }
    }
  
  if ( TargetImageOffsetReal )
    {
    targetVolume->SetOffset( TargetImageOffset );
    
    if ( Verbose )
      {
      cmtk::StdOut << "INFO: setting reference image offset to " << TargetImageOffset[0] << "/" << TargetImageOffset[1] << "/" << TargetImageOffset[2] << "\n";
      }
    }
  
  cmtk::UniformVolume::SmartPtr floatingVolume;
  if ( FloatingVolumeName )
    {
    floatingVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( FloatingVolumeName, Verbose ) );
    if ( ! floatingVolume ) 
      {
      cmtk::StdErr << "ERROR: floating volume " << FloatingVolumeName << " could not be read\n";
      throw cmtk::ExitException( 1 );
      }
    }
  
  cmtk::TypedArray::SmartPtr floatingData( NULL );
  if ( floatingVolume ) 
    {
    floatingData = floatingVolume->GetData();
    if ( ! floatingData ) 
      {
      cmtk::StdErr << "ERROR: floating volume " << FloatingVolumeName << " seems to have no data\n";
      throw cmtk::ExitException( 1 );
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
      cmtk::StdOut << "INFO: Using target data as binary mask.\n";
    }
  
  cmtk::ProgressConsole progressIndicator;

  cmtk::TypedArray::SmartPtr reformatData;
  switch ( Interpolation ) 
    {
    default:
    case cmtk::Interpolators::LINEAR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> TInterpolator;
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
    break;
    }
    case cmtk::Interpolators::CUBIC:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Cubic> TInterpolator;
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
    break;
    }
    case cmtk::Interpolators::NEAREST_NEIGHBOR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor> TInterpolator;
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
    break;
    }
    case cmtk::Interpolators::COSINE_SINC:
    {
    switch ( InterpolatorWindowRadius )
      {
      case 2:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<2> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
      break;
      }
      case 3:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<3> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
      break;
      }
      case 4:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<4> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
      break;
      }
      case 5:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<5> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
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
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
      break;
      }
      case 3:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<3> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
      break;
      }
      case 4:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<4> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
      break;
      }
      case 5:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<5> > TInterpolator;
      InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
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
    InitializeReformatVolume<TInterpolator>( reformatData, targetVolume, floatingVolume );
    break;
    }
    }
  
  if ( reformatData ) 
    {
    targetVolume->SetData( reformatData );
    cmtk::VolumeIO::Write( *targetVolume, OutputImageName, Verbose );

#ifdef CMTK_USE_SQLITE
    if ( updateDB )
      {
      if ( TargetVolumeName )
	{
	cmtk::ImageXformDB db( updateDB );
	db.AddImage( OutputImageName, TargetVolumeName );
	}
      }
#endif
    } 
  else 
    {
    cmtk::StdErr << "ERROR: no reformatted data was generated.\n";
    }
}

int 
doMain( const int argc, const char* argv[] )
{
  cmtk::Threads::CheckEnvironment(); // need this to check for "CMTK_NUM_THREADS" and constrain OpenMP accordingly
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Volume reformatter" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Extended volume reformatter tool to compute reformatted images and Jacobian maps from arbitrary sequences of concatenated transformations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] --floating floatingImg target x0 [x1 ...]\n  OR\n"
		       "[options] target [x0 x1 ...] {-j,--jacobian} xx0 [xx1 ...]\n  OR\n"
		       "WHERE x0 ... xN and xx0 ... xxN is [{-i,--inverse}] transformation##" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.BeginGroup( "PlainOptions", "Options for Plain Reformatting" );
    cmtk::CommandLine::EnumGroup<cmtk::Interpolators::InterpolationEnum>::SmartPtr
      interpolationGroup = cl.AddEnum( "interpolation", &Interpolation, "Image interpolation method." );
    interpolationGroup->AddSwitch( Key( "linear" ), cmtk::Interpolators::LINEAR, "Trilinear interpolation" );
    interpolationGroup->AddSwitch( Key( "nn" ), cmtk::Interpolators::NEAREST_NEIGHBOR, "Nearest neighbor interpolation" );
    interpolationGroup->AddSwitch( Key( "cubic" ), cmtk::Interpolators::CUBIC, "Tricubic interpolation" );
    interpolationGroup->AddSwitch( Key( "pv" ), cmtk::Interpolators::PARTIALVOLUME, "Partial volume interpolation" );
    interpolationGroup->AddSwitch( Key( "sinc-cosine" ), cmtk::Interpolators::COSINE_SINC, "Sinc interpolation with cosine window" );
    interpolationGroup->AddSwitch( Key( "sinc-hamming" ), cmtk::Interpolators::HAMMING_SINC, "Sinc interpolation with Hamming window" );

    cl.AddOption( Key( "sinc-window-radius" ), &InterpolatorWindowRadius, "Window radius for Sinc interpolation" );

    cl.AddSwitch( Key( "preserve-mass" ), &MassPreservingReformat, true, "Mass-preserving reformatting: multiply every reformatted value with the Jacobian determinant of the applied transformation." );
    cl.EndGroup();

    cl.BeginGroup( "JacobianOptions", "Options for Jacobian Map Reformatting" );
    cl.AddSwitch( Key( "jacobian-correct-global" ), &JacobianCorrectGlobal, true, "Correct Jacobian maps for global scale." );
    cl.AddSwitch( Key( "no-jacobian-correct-global" ), &JacobianCorrectGlobal, false, "Do not correct Jacobian maps for global scale." );
    cl.EndGroup();

    cl.BeginGroup( "Input", "Input Options" );
    cl.AddCallback( Key( "target-grid" ), CallbackTargetVolume, "Define target grid for reformating as Nx,Ny,Nz:dX,dY,dZ[:Ox,Oy,Oz] (dims:pixel:offset)" );
    cl.AddSwitch( Key( 'm', "mask" ), &TargetMask, true, "Use target pixel data as binary mask." );
    cl.AddCallback( Key( "crop" ), CallbackCropImages, "Crop target image: x0,y0,z0,x1,y1,z2" );
    cl.AddCallback( Key( 'O', "target-offset" ), CallbackTargetImageOffset, "Override target image offset and set to x,y,z mm" );
    cl.AddCallback( Key( "target-offset-pixels" ), CallbackTargetImageOffsetPixels, "Override target image offset and set to dx,dy,dz pixels" );
    cl.AddOption( Key( 'F', "floating" ), &FloatingVolumeName, "Format and path of floating image." );

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'o', "outfile" ), &OutputImageName, "Format and path of output image." );
    cl.EndGroup();

    cmtk::CommandLine::EnumGroup<cmtk::ScalarDataType>::SmartPtr
      typeGroup = cl.AddEnum( "outputtype", &DataType, "Scalar data type for the output image." );
    typeGroup->AddSwitch( Key( "char" ), cmtk::TYPE_CHAR, "8 bits, signed" );
    typeGroup->AddSwitch( Key( "byte" ), cmtk::TYPE_BYTE, "8 bits, unsigned" );
    typeGroup->AddSwitch( Key( "short" ), cmtk::TYPE_SHORT, "16 bits, signed" );
    typeGroup->AddSwitch( Key( "ushort" ), cmtk::TYPE_USHORT, "16 bits, unsigned" );
    typeGroup->AddSwitch( Key( "int" ), cmtk::TYPE_INT, "32 bits signed" );
    typeGroup->AddSwitch( Key( "uint" ), cmtk::TYPE_UINT, "32 bits unsigned" );
    typeGroup->AddSwitch( Key( "float" ), cmtk::TYPE_FLOAT, "32 bits floating point" );
    typeGroup->AddSwitch( Key( "double" ), cmtk::TYPE_DOUBLE, "64 bits floating point\n" );

    cl.AddOption( Key( 'P', "pad-out" ), &OutPaddingValue, "Padding value for output image", &OutPaddingValueFlag );
    cl.EndGroup();

#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &updateDB, "Path to image/transformation database that should be updated with the newly created image." );
    cl.EndGroup();
#endif

    cl.Parse( argc, argv );

    if ( ! UserDefinedTargetVolume )
      TargetVolumeName = cl.GetNext();
    
    const char* next = cl.GetNextOptional();
    while ( next ) 
      {
      if ( ! strcmp( next, "-j" ) || ! strcmp( next, "--jacobian" ) )
	break;
      
      const bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
      if ( inverse ) 
	next = cl.GetNext();
      
      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
      if ( ! xform ) 
	{
	cmtk::StdErr << "ERROR: could not read target-to-reference transformation from " << next << "\n";
	throw cmtk::ExitException( 1 );
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
	    const bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
	    if ( inverse ) 
	      next = cl.GetNext();
	    
	    cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
	    if ( ! xform ) 
	      {
	      cmtk::StdErr << "ERROR: could not read target-to-floating transformation from " << next << "\n";
	      throw cmtk::ExitException( 1 );
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
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  ReformatPullback();
  
  return 0;
}

#include "cmtkSafeMain"
