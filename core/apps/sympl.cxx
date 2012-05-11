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
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkProgress.h>
#include <System/cmtkProgressConsole.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkUnits.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkVolume.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUniformVolumeInterpolator.h>
#include <Base/cmtkSincInterpolator.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkCubicInterpolator.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkClassStreamOutput.h>

#include <Registration/cmtkSymmetryPlaneFunctional.h>
#include <Registration/cmtkBestNeighbourOptimizer.h>
#include <Registration/cmtkReformatVolume.h>

#include <stdio.h>

float MinValue = -1e5;
bool MinValueSet = false;
float MaxValue = 1e5;
bool MaxValueSet = false;

cmtk::Types::Coordinate Sampling = 1.0;
cmtk::Types::Coordinate Accuracy = 0.1;

cmtk::Interpolators::InterpolationEnum Interpolation = cmtk::Interpolators::LINEAR;

int Levels = 4;

bool OutputOnly = false;
cmtk::Types::Coordinate Rho;
cmtk::Units::Degrees Theta;
cmtk::Units::Degrees Phi;

bool DoWriteMirror = false;
const char* MirrorOutFile = "symmetry_mirror.nii";

bool DoWriteAligned = false;
const char* AlignedOutFile = "symmetry_aligned.nii";
bool MarkPlaneAligned = false;

bool DoWriteMarked = false;
const char* MarkedOutFile = "symmetry_marked.nii";

bool DoWriteDifference = false;
const char* DifferenceOutFile = "symmetry_diff.nii";

const char* WriteXformPath = NULL;

cmtk::Types::DataItem MarkPlaneValue = 4095;

bool PadOutValueSet = false;
cmtk::Types::DataItem PadOutValue = 0;

const char* SymmetryOutFileName = NULL;

const char* SymmetryParameters = NULL;
const char* SymmetryParametersFile = NULL;

const char* InFileName = NULL;

/// Constants for initial plane orientation.
typedef enum
{
  /// XY plane (axial)
  SYMPL_INIT_XY,
  /// XZ plane (coronal)
  SYMPL_INIT_XZ,
  /// YZ plane (sagittal)
  SYMPL_INIT_YZ
} InitialPlaneEnum;

/// Initial plane orientation: default to sagittal for human images.
InitialPlaneEnum InitialPlane = SYMPL_INIT_YZ;

bool
ParseCommandLine ( const int argc, const char* argv[] )
{
  try
    {
    cmtk::CommandLine cl( cmtk::CommandLine::PROPS_XML  );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Symmetry plane computation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Compute the approximate symmetry plane of an image to determine, for example, the mid-sagittal plane in human brain images. "
		       "Various forms of output are supported, e.g., writing the input image with the symmetry plane drawn into it, or the input image realigned along the symmetry plane." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Registration" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Optimization", "Optimization" );
    cl.AddOption( Key( 'a', "accuracy" ), &Accuracy, "Accuracy (final optimization step size in [mm]." );
    cl.AddOption( Key( 's', "sampling" ), &Sampling, "Resampled image resolution. This is the resolution [in mm] of the first (finest) resampled image in the multi-scale pyramid, "
		  "which is derived directly from the original full-resolution images.");
    cl.AddOption( Key( 'l', "levels" ), &Levels, "Number of resolution levels. The algorithm will create (levels-1) resampled images with increasingly coarse resolution and use these "
		  "in successive order of increasing resolution before using the original images at the final level." );
    cl.EndGroup();
    
    cl.BeginGroup( "Initial", "Initial approximate symmetry plane orientation" );
    cmtk::CommandLine::EnumGroup<InitialPlaneEnum>::SmartPtr
      initialPlane = cl.AddEnum( "initial-plane", &InitialPlane, "Initial orientation of symmetry plane. This should be the closest orthogonal plane to the expected actual symmetry plane." );
    initialPlane->AddSwitch( Key( "initial-axial" ), SYMPL_INIT_XY, "Approximately axial symmetry" );
    initialPlane->AddSwitch( Key( "initial-coronal" ), SYMPL_INIT_XZ, "Approximately coronal symmetry" );
    initialPlane->AddSwitch( Key( "initial-sagittal" ), SYMPL_INIT_YZ, "Approximately sagittal symmetry" );
    initialPlane->AddSwitch( Key( "initial-xy" ), SYMPL_INIT_XY, "Approximately XY plane symmetry" );
    initialPlane->AddSwitch( Key( "initial-xz" ), SYMPL_INIT_XZ, "Approximately XZ plane symmetry" );
    initialPlane->AddSwitch( Key( "initial-yz" ), SYMPL_INIT_YZ, "Approximately YZ plane symmetry" );
    cl.EndGroup();
    
    cl.BeginGroup( "Pre-computed", "Pre-computed symmetry" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED | cmtk::CommandLine::PROPS_NOXML );
    cl.AddOption( Key( "output-only" ), &SymmetryParameters, "Give symmetry parameters [Rho Theta Phi] as option, skip search.", &OutputOnly );
    cl.AddOption( Key( "output-only-file" ), &SymmetryParametersFile, "Read symmetry parameters from file, skip search.", &OutputOnly );
    cl.EndGroup();
    
    cl.BeginGroup( "Preprocessing", "Data pre-processing" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( "min-value" ), &MinValue, "Force minumum data value.", &MinValueSet );
    cl.AddOption( Key( "max-value" ), &MaxValue, "Force maximum data value.", &MaxValueSet );
    cl.EndGroup();
    
    cl.BeginGroup( "OutputImages", "Output of Images" );
    cmtk::CommandLine::EnumGroup<cmtk::Interpolators::InterpolationEnum>::SmartPtr
      interpGroup = cl.AddEnum( "interpolation", &Interpolation, "Interpolation method used for reformatted output data" );
    interpGroup->AddSwitch( Key( 'L', "linear" ), cmtk::Interpolators::LINEAR, "Use linear image interpolation for output." );
    interpGroup->AddSwitch( Key( 'C', "cubic" ), cmtk::Interpolators::CUBIC, "Use cubic image interpolation for output." );
    interpGroup->AddSwitch( Key( 'S', "sinc" ), cmtk::Interpolators::COSINE_SINC, "Use cosine-windowed sinc image interpolation for output." );
    
    cl.AddOption( Key( 'P', "pad-out" ), &PadOutValue, "Padding value for output images.", &PadOutValueSet )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( "mark-value" ), &MarkPlaneValue, "Data value to mark (draw) symmetry plane.", &DoWriteMarked );
    cl.AddOption( Key( "write-marked" ), &MarkedOutFile, "File name for output image with marked symmetry plane.", &DoWriteMarked )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-aligned" ), &AlignedOutFile, "File name for symmetry plane-aligned output image.", &DoWriteAligned )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddSwitch( Key( "mark-aligned" ), &MarkPlaneAligned, true, "Mark symmetry plane in aligned output image." );
    cl.AddOption( Key( "write-subtract" ), &DifferenceOutFile, "File name for mirror subtraction image.", &DoWriteDifference )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-mirror" ), &MirrorOutFile, "File name for image mirrored w.r.t. symmetry plane.", &DoWriteMirror )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.EndGroup();

    cl.BeginGroup( "OutputParameters", "Output of Parameters" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( 'o', "outfile" ), &SymmetryOutFileName, "File name for symmetry plane parameter output." )->SetProperties( cmtk::CommandLine::PROPS_FILENAME | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-xform" ), &WriteXformPath, "Write affine alignment transformation to file" )
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT )
      ->SetAttribute( "reference", "InputImage" );
    cl.EndGroup();
    
    cl.AddParameter( &InFileName, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    
    if ( ! cl.Parse( argc, argv ) ) return false;
    
    if ( SymmetryParameters ) 
      {
      double rho, theta, phi;
      if ( 3 == sscanf( SymmetryParameters, "%lf %lf %lf", &rho, &theta, &phi ) ) 
	{
	Rho = rho; 
	Theta = cmtk::Units::Degrees( theta );
	Phi = cmtk::Units::Degrees( phi );
	}
      }
    
    if ( SymmetryParametersFile ) 
      {
      cmtk::ClassStreamInput inStream( SymmetryParametersFile );
      if ( inStream.IsValid() ) 
	{
	cmtk::ParametricPlane *plane = NULL;
	inStream >> plane;
	Rho = plane->GetRho(); 
	Theta = plane->GetTheta();
	Phi = plane->GetPhi();
	delete plane;
	} 
      else
	{
	cmtk::StdErr.printf( "ERROR: Could not open symmetry parameter file %s\n", SymmetryParametersFile );
	}
      }
    }
  catch ( const cmtk::CommandLine::Exception& ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  return true;
}

void
WriteDifference
( const cmtk::UniformVolume* originalVolume, const cmtk::UniformVolumeInterpolatorBase* interpolator, const cmtk::ParametricPlane& parametricPlane, const char* outFileName )
{
  cmtk::UniformVolume::SmartPtr diffVolume( originalVolume->CloneGrid() );
  const cmtk::TypedArray* originalData = originalVolume->GetData();
  cmtk::TypedArray::SmartPtr diffData = cmtk::TypedArray::SmartPtr( cmtk::TypedArray::Create( GetSignedDataType( originalData->GetType() ), originalData->GetDataSize() ) );
  diffVolume->SetData( diffData );

  cmtk::Types::DataItem dataV, dataW;

  int offset = 0;
  for ( int z = 0; z < originalVolume->GetDims()[2]; ++z )
    for ( int y = 0; y < originalVolume->GetDims()[1]; ++y )
      for ( int x = 0; x < originalVolume->GetDims()[0]; ++x, ++offset ) 
	{
	if ( ! originalData->Get( dataV, offset ) ) 
	  {
	  diffData->SetPaddingAt( offset );
	  continue;
	  }
	cmtk::UniformVolume::CoordinateVectorType w = originalVolume->GetGridLocation( x, y, z );
	parametricPlane.MirrorInPlace( w );

	if ( interpolator->GetDataAt( w, dataW ) )
	  {
	  diffData->Set( fabs( dataV - dataW ), offset );
	  }
	else
	  {
	  diffData->SetPaddingAt( offset );
	  }
	}
  
  cmtk::VolumeIO::Write( *diffVolume, outFileName );
}

void
WriteMirror
( const cmtk::UniformVolume* originalVolume, const cmtk::UniformVolumeInterpolatorBase* interpolator, const cmtk::ParametricPlane& parametricPlane, const char* outFileName )
{
  cmtk::TypedArray::SmartPtr mirrorData = cmtk::TypedArray::Create( originalVolume->GetData()->GetType(), originalVolume->GetData()->GetDataSize() );

  cmtk::Types::DataItem data;

  int offset = 0;
  for ( int z = 0; z < originalVolume->GetDims()[2]; ++z ) 
    {
    for ( int y = 0; y < originalVolume->GetDims()[1]; ++y )
      for ( int x = 0; x < originalVolume->GetDims()[0]; ++x, ++offset ) 
	{
	cmtk::UniformVolume::CoordinateVectorType v = originalVolume->GetGridLocation( x, y, z );
	parametricPlane.MirrorInPlace( v );

	if ( interpolator->GetDataAt( v, data ) )
	  {
	  mirrorData->Set( data, offset );
	  }
	else
	  {
	  mirrorData->SetPaddingAt( offset );
	  }
	}
    }

  cmtk::UniformVolume::SmartPtr mirrorVolume( originalVolume->CloneGrid() );
  mirrorVolume->SetData( mirrorData );
  cmtk::VolumeIO::Write( *mirrorVolume, outFileName );
}

void
WriteMarkPlane
( const cmtk::UniformVolume* originalVolume, const cmtk::ParametricPlane& parametricPlane, const cmtk::Types::DataItem markPlaneValue, const char* outFileName )
{
  cmtk::UniformVolume::SmartPtr markVolume( originalVolume->CloneGrid() );
  cmtk::TypedArray::SmartPtr markData( originalVolume->GetData()->Clone() );
  markVolume->SetData( markData );

  int offset = 0;
  for ( int z = 0; z < originalVolume->GetDims()[2]; ++z ) 
    {
    for ( int y = 0; y < originalVolume->GetDims()[1]; ++y ) 
      {
      int currentSideOfPlane = 0;
      for ( int x = 0; x < originalVolume->GetDims()[0]; ++x, ++offset ) 
	{
	int newSideOfPlane = parametricPlane.GetWhichSide( originalVolume->GetGridLocation( x, y, z ) );
	if ( ( newSideOfPlane != currentSideOfPlane ) && x )
	  markData->Set( markPlaneValue, offset );
	currentSideOfPlane = newSideOfPlane;
	}
      }    
    }

  cmtk::VolumeIO::Write( *markVolume, outFileName );
}

void
WriteAligned
( const cmtk::UniformVolume* originalVolume, const cmtk::UniformVolumeInterpolatorBase* interpolator, const cmtk::ParametricPlane& parametricPlane, const InitialPlaneEnum initialPlane, 
  const char* outFileName )
{
  const cmtk::TypedArray* originalData = originalVolume->GetData();

  cmtk::TypedArray::SmartPtr alignData = cmtk::TypedArray::Create( originalData->GetType(), originalData->GetDataSize() );
  if ( PadOutValueSet )
    {
    alignData->SetPaddingValue( PadOutValue );
    }

  cmtk::UniformVolume::SmartPtr alignVolume( originalVolume->CloneGrid() );
  alignVolume->SetData( alignData );

  const cmtk::Types::DataItem maxData = originalData->GetRange().m_UpperBound;

  cmtk::Types::DataItem data;

  int normalAxis = 0;
  switch ( initialPlane )
    {
    case SYMPL_INIT_XY: normalAxis = 2; break;
    case SYMPL_INIT_XZ: normalAxis = 1; break;
    case SYMPL_INIT_YZ: normalAxis = 0; break;
    }

  cmtk::AffineXform::SmartPtr alignment( parametricPlane.GetAlignmentXform( normalAxis ) );
  int offset = 0;
  for ( int z = 0; z < originalVolume->GetDims()[2]; ++z ) 
    {
    for ( int y = 0; y < originalVolume->GetDims()[1]; ++y ) 
      {
      for ( int x = 0; x < originalVolume->GetDims()[0]; ++x, ++offset ) 
	{
	cmtk::UniformVolume::CoordinateVectorType v = originalVolume->GetGridLocation( x, y, z );
	alignment->ApplyInPlace( v );

	if ( interpolator->GetDataAt( v, data ) )
	  {
	  if ( MarkPlaneAligned && (x == ( originalVolume->GetDims()[0] / 2 )) )
	    alignData->Set( 2 * maxData, offset );
	  else
	    alignData->Set( data, offset );
	  }
	else
	  {
	  alignData->SetPaddingAt( offset );
	  }
	}
      }
    }

  cmtk::VolumeIO::Write( *alignVolume, outFileName );
}

int
doMain ( const int argc, const char* argv[] ) 
{
  if ( ! ParseCommandLine( argc, argv ) ) return 1;

  cmtk::UniformVolume::SmartPtr originalVolume( cmtk::VolumeIO::ReadOriented( InFileName ) );
  if ( !originalVolume ) 
    {
    cmtk::StdErr.printf( "Could not read image file %s\n", InFileName );
    throw cmtk::ExitException(1);
    }

  cmtk::CoordinateVector v( 6 );
  // initialize plane as the mid-sagittal with respect to image orientation --
  // distance from coordinate origin (image center) is 0:
  v[0] = 0;
  // and angles are chosen so that the plane normal is (1,0,0)
  switch ( InitialPlane )
    {
    case SYMPL_INIT_XY:
      v[1] = 0;
      v[2] = 0;
      break;
    case SYMPL_INIT_XZ:
      v[1] = 90;
      v[2] = 90;
      break;
    default:
    case SYMPL_INIT_YZ:
      v[1] = 0;
      v[2] = 90;
      break;
    }
  
  // set center of volume (crop region) as coordinate origin.
  cmtk::Vector3D center = originalVolume->GetCenterCropRegion();
  v[3] = center[0]; v[4] = center[1]; v[5] = center[2];

  if ( OutputOnly ) 
    {
    v[0] = Rho;
    v[1] = Theta.Value();
    v[2] = Phi.Value();
    } 
  else
    {
    cmtk::BestNeighbourOptimizer optimizer;
    
    cmtk::ProgressConsole progressIndicator( "Symmetry Plane Computation" );
    cmtk::Progress::Begin( 0, Levels, 1, "Symmetry Plane Computation" );

    for ( int level = 0; level < Levels; ++level ) 
      {      
      cmtk::UniformVolume::SmartPtr volume;
      if ( level < Levels-1 ) 
	{
	cmtk::Types::Coordinate voxelSize = Sampling * pow( 2.0, (Levels-level-2) );
	volume = cmtk::UniformVolume::SmartPtr( new cmtk::UniformVolume( *originalVolume, voxelSize ) );
	cmtk::DebugOutput( 1 ).GetStream().printf( "Entering level %d out of %d (%.2f mm voxel size)\n", level+1, Levels, voxelSize );
	} 
      else
	{
	volume = originalVolume; 
	cmtk::DebugOutput( 1 ).GetStream().printf( "Entering level %d out of %d (original voxel size)\n", level+1, Levels );
	}
      
      cmtk::SmartPointer<cmtk::SymmetryPlaneFunctional> functional( NULL );
      if ( MinValueSet || MaxValueSet ) 
	{
	cmtk::Types::DataItemRange valueRange = volume->GetData()->GetRange();
	
	if ( MinValueSet ) 
	  valueRange.m_LowerBound = MinValue;
	if ( MaxValueSet ) 
	  valueRange.m_UpperBound = MaxValue;
	
	functional = cmtk::SmartPointer<cmtk::SymmetryPlaneFunctional>( new cmtk::SymmetryPlaneFunctional( volume, valueRange ) );
	} 
      else
	{
	functional = cmtk::SmartPointer<cmtk::SymmetryPlaneFunctional>( new cmtk::SymmetryPlaneFunctional( volume ) );
	}
      
      optimizer.SetFunctional( cmtk::Functional::SmartPtr::DynamicCastFrom( functional ) );
      optimizer.Optimize( v, pow( 2.0, Levels-level-1 ), Accuracy * pow( 2.0, Levels-level-1 ) );

      cmtk::Progress::SetProgress( level );      
      }

    cmtk::Progress::Done();

    cmtk::DebugOutput( 1 ).GetStream().printf( "rho=%f, theta=%f, phi=%f\n", v[0], v[1], v[2] );
    }
  
  cmtk::ParametricPlane parametricPlane;
  parametricPlane.SetParameters( v );
  
  if ( SymmetryOutFileName )
    {
    cmtk::ClassStreamOutput stream( SymmetryOutFileName, cmtk::ClassStreamOutput::MODE_WRITE );
    stream << parametricPlane;
    stream.Close();
    }

  const cmtk::UniformVolumeInterpolatorBase::SmartPtr interpolator( cmtk::ReformatVolume::CreateInterpolator( Interpolation, originalVolume ) );;
  
  if ( DoWriteAligned ) 
    WriteAligned( originalVolume, interpolator, parametricPlane, InitialPlane, AlignedOutFile );

  if ( DoWriteMarked ) 
    WriteMarkPlane( originalVolume, parametricPlane, MarkPlaneValue, MarkedOutFile );

  if ( DoWriteDifference )
    WriteDifference( originalVolume, interpolator, parametricPlane, DifferenceOutFile );

  if ( DoWriteMirror )
    WriteMirror( originalVolume, interpolator, parametricPlane, MirrorOutFile );

  if ( WriteXformPath )
    {
    cmtk::AffineXform::SmartPtr alignment( parametricPlane.GetAlignmentXform( 0 ) );
    cmtk::XformIO::Write( alignment, WriteXformPath );
    }

  return 0;
}

#include "cmtkSafeMain"
