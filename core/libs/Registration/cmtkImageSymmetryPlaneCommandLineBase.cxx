/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkImageSymmetryPlaneCommandLineBase.h"

#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>

#include <Registration/cmtkBestNeighbourOptimizer.h>
#include <Registration/cmtkReformatVolume.h>

#include <System/cmtkProgress.h>
#include <System/cmtkProgressConsole.h>
#include <System/cmtkDebugOutput.h>

cmtk::ImageSymmetryPlaneCommandLineBase
::ImageSymmetryPlaneCommandLineBase()
  : m_MinValue( 0.0 ), m_MinValueSet( false ),
    m_MaxValue( 0.0 ), m_MaxValueSet( false ),
    m_Sampling( 1.0 ),
    m_Accuracy( 0.1 ),
    m_Interpolation( Interpolators::LINEAR ),
    m_Levels( 4 ),
    m_DisableOptimization( false ),
    m_FixOffset( false ),
    m_MirrorOutFile( NULL ),
    m_AlignedOutFile( NULL ),
    m_MarkPlaneAligned( false ),
    m_MarkedOutFile( NULL ),
    m_DifferenceOutFile( NULL ),
    m_WriteXformPath( NULL ),
    m_MarkPlaneValue( 4095 ),
    m_PadOutValueSet( false ),
    m_SymmetryOutFileName( NULL ),
    m_SymmetryParameters( NULL ),
    m_SymmetryParametersFile( NULL ),
    m_InitialPlane( SYMPL_INIT_YZ ),
    m_CommandLine( CommandLine::PROPS_XML )
{
  this->m_CommandLine.SetProgramInfo( CommandLine::PRG_TITLE, "Symmetry plane computation" );
  this->m_CommandLine.SetProgramInfo( CommandLine::PRG_DESCR, "Compute the approximate symmetry plane of an image to determine, for example, the mid-sagittal plane in human brain images. "
				      "Various forms of output are supported, e.g., writing the input image with the symmetry plane drawn into it, or the input image realigned along the symmetry plane." );
  this->m_CommandLine.SetProgramInfo( CommandLine::PRG_CATEG, "CMTK.Registration" );
  
  typedef CommandLine::Key Key;

  this->m_CommandLine.BeginGroup( "Optimization", "Optimization" );
  this->m_CommandLine.AddOption( Key( 'a', "accuracy" ), &this->m_Accuracy, "Accuracy (final optimization step size in [mm]." );
  this->m_CommandLine.AddOption( Key( 's', "sampling" ), &this->m_Sampling, "Resampled image resolution. This is the resolution [in mm] of the first (finest) resampled image in the multi-scale pyramid, "
				 "which is derived directly from the original full-resolution images.");
  this->m_CommandLine.AddOption( Key( 'l', "levels" ), &this->m_Levels, "Number of resolution levels. The algorithm will create (levels-1) resampled images with increasingly coarse resolution and use these "
				 "in successive order of increasing resolution before using the original images at the final level." );
  this->m_CommandLine.AddSwitch( Key( "fix-offset" ), &this->m_FixOffset, true, "Fix symmetry plane offset. Reduces computation time and forces symmetry plane to cross center of mass, but may lead to less-than-accurate result." );
  this->m_CommandLine.EndGroup();
  
  this->m_CommandLine.BeginGroup( "Initial", "Initial approximate symmetry plane orientation" );
  CommandLine::EnumGroup<InitialPlaneEnum>::SmartPtr
    initialPlane = this->m_CommandLine.AddEnum( "initial-plane", &this->m_InitialPlane, "Initial orientation of symmetry plane. This should be the closest orthogonal plane to the expected actual symmetry plane." );
  initialPlane->AddSwitch( Key( "initial-axial" ), SYMPL_INIT_XY, "Approximately axial symmetry" );
  initialPlane->AddSwitch( Key( "initial-coronal" ), SYMPL_INIT_XZ, "Approximately coronal symmetry" );
  initialPlane->AddSwitch( Key( "initial-sagittal" ), SYMPL_INIT_YZ, "Approximately sagittal symmetry" );
  initialPlane->AddSwitch( Key( "initial-xy" ), SYMPL_INIT_XY, "Approximately XY plane symmetry" );
  initialPlane->AddSwitch( Key( "initial-xz" ), SYMPL_INIT_XZ, "Approximately XZ plane symmetry" );
  initialPlane->AddSwitch( Key( "initial-yz" ), SYMPL_INIT_YZ, "Approximately YZ plane symmetry" );
  this->m_CommandLine.EndGroup();
  
  this->m_CommandLine.BeginGroup( "Pre-computed", "Pre-computed symmetry" )->SetProperties( CommandLine::PROPS_ADVANCED | CommandLine::PROPS_NOXML );
  this->m_CommandLine.AddOption( Key( "output-only" ), &this->m_SymmetryParameters, "Give symmetry parameters [Rho Theta Phi] as option, skip search.", &this->m_DisableOptimization );
  this->m_CommandLine.AddOption( Key( "output-only-file" ), &this->m_SymmetryParametersFile, "Read symmetry parameters from file, skip search.", &this->m_DisableOptimization );
  this->m_CommandLine.EndGroup();
  
  this->m_CommandLine.BeginGroup( "Preprocessing", "Data pre-processing" )->SetProperties( CommandLine::PROPS_ADVANCED );
  this->m_CommandLine.AddOption( Key( "min-value" ), &this->m_MinValue, "Force minumum data value.", &this->m_MinValueSet );
  this->m_CommandLine.AddOption( Key( "max-value" ), &this->m_MaxValue, "Force maximum data value.", &this->m_MaxValueSet );
  this->m_CommandLine.EndGroup();
  
  this->m_CommandLine.BeginGroup( "OutputImages", "Output of Images" );
  CommandLine::EnumGroup<Interpolators::InterpolationEnum>::SmartPtr interpGroup = this->m_CommandLine.AddEnum( "interpolation", &this->m_Interpolation, "Interpolation method used for reformatted output data" );
  interpGroup->AddSwitch( Key( 'L', "linear" ), Interpolators::LINEAR, "Use linear image interpolation for output." );
  interpGroup->AddSwitch( Key( 'C', "cubic" ), Interpolators::CUBIC, "Use cubic image interpolation for output." );
  interpGroup->AddSwitch( Key( 'S', "sinc" ), Interpolators::COSINE_SINC, "Use cosine-windowed sinc image interpolation for output." );
  
  this->m_CommandLine.AddOption( Key( 'P', "pad-out" ), &this->m_PadOutValue, "Padding value for output images.", &this->m_PadOutValueSet )->SetProperties( CommandLine::PROPS_ADVANCED );
  this->m_CommandLine.AddOption( Key( "mark-value" ), &this->m_MarkPlaneValue, "Data value to mark (draw) symmetry plane." );
  this->m_CommandLine.AddOption( Key( "write-marked" ), &this->m_MarkedOutFile, "File name for output image with marked symmetry plane." )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OUTPUT );
  this->m_CommandLine.AddOption( Key( "write-aligned" ), &this->m_AlignedOutFile, "File name for symmetry plane-aligned output image." )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OUTPUT );
  this->m_CommandLine.AddSwitch( Key( "mark-aligned" ), &this->m_MarkPlaneAligned, true, "Mark symmetry plane in aligned output image." );
  this->m_CommandLine.AddOption( Key( "write-subtract" ), &this->m_DifferenceOutFile, "File name for mirror subtraction image." )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OUTPUT );
  this->m_CommandLine.AddOption( Key( "write-mirror" ), &this->m_MirrorOutFile, "File name for image mirrored w.r.t. symmetry plane." )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OUTPUT );
  this->m_CommandLine.EndGroup();
  
  this->m_CommandLine.BeginGroup( "OutputParameters", "Output of Parameters" )->SetProperties( CommandLine::PROPS_ADVANCED );
  this->m_CommandLine.AddOption( Key( 'o', "outfile" ), &this->m_SymmetryOutFileName, "File name for symmetry plane parameter output." )->SetProperties( CommandLine::PROPS_FILENAME | CommandLine::PROPS_OUTPUT );
  this->m_CommandLine.AddOption( Key( "write-xform" ), &this->m_WriteXformPath, "Write affine alignment transformation to file" )
    ->SetProperties( CommandLine::PROPS_XFORM | CommandLine::PROPS_OUTPUT )
    ->SetAttribute( "reference", "InputImage" );
  this->m_CommandLine.EndGroup();
  
  this->m_CommandLine.AddParameter( &this->m_InFileName, "InputImage", "Input image path" )->SetProperties( CommandLine::PROPS_IMAGE );
}

int
cmtk::ImageSymmetryPlaneCommandLineBase
::Run( const int argc, const char* argv[] )
{
  if ( ! this->ParseCommandLine( argc, argv ) )
    return 2;
  
  UniformVolume::SmartPtr originalVolume( VolumeIO::ReadOriented( this->m_InFileName ) );
  if ( !originalVolume ) 
    {
    StdErr.printf( "Could not read image file %s\n", this->m_InFileName.c_str() );
    return 1;
    }

  CoordinateVector v( 6 );
  // initialize plane as the mid-sagittal with respect to image orientation --
  // distance from coordinate origin (image center) is 0:
  v[0] = 0;
  // and angles are chosen so that the plane normal is (1,0,0)
  switch ( this->m_InitialPlane )
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
  Vector3D center = originalVolume->GetCenterOfMass();
  v[3] = center[0]; v[4] = center[1]; v[5] = center[2];

  if ( this->m_DisableOptimization ) 
    {
    v[0] = this->m_Rho;
    v[1] = this->m_Theta.Value();
    v[2] = this->m_Phi.Value();
    } 
  else
    {
    BestNeighbourOptimizer optimizer;

    // Instantiate programm progress indicator.
    ProgressConsole progressIndicator( "Symmetry Plane Computation" );
    Progress::Begin( 0, this->m_Levels, 1, "Symmetry Plane Computation" );
    
    for ( int level = 0; level < this->m_Levels; ++level ) 
      {      
      UniformVolume::SmartPtr volume;
      if ( level < this->m_Levels-1 ) 
	{
	Types::Coordinate voxelSize = this->m_Sampling * pow( 2.0, (this->m_Levels-level-2) );
	DebugOutput( 1 ).GetStream().printf( "Entering level %d out of %d (%.2f mm voxel size)\n", level+1, this->m_Levels, voxelSize );
	volume = UniformVolume::SmartPtr( originalVolume->GetResampled( voxelSize ) );
	} 
      else
	{
	DebugOutput( 1 ).GetStream().printf( "Entering level %d out of %d (original voxel size)\n", level+1, this->m_Levels );
	volume = originalVolume; 
	}
      
      ImageSymmetryPlaneFunctionalBase::SmartPtr functional;
      if ( this->m_MinValueSet || this->m_MaxValueSet ) 
	{
	Types::DataItemRange valueRange = volume->GetData()->GetRange();
	
	if ( this->m_MinValueSet ) 
	  valueRange.m_LowerBound = this->m_MinValue;
	if ( this->m_MaxValueSet ) 
	  valueRange.m_UpperBound = this->m_MaxValue;
	
	functional = this->CreateFunctional( volume, valueRange );
	} 
      else
	{
	functional = this->CreateFunctional( volume );
	}
      
      functional->SetFixOffset( this->m_FixOffset );
      
      optimizer.SetFunctional( functional );
      optimizer.Optimize( v, pow( 2.0, this->m_Levels-level-1 ), this->m_Accuracy * pow( 2.0, this->m_Levels-level-1 ) );

      Progress::SetProgress( level );
      }

    Progress::Done();

    DebugOutput( 1 ).GetStream().printf( "rho=%f, theta=%f, phi=%f\n", v[0], v[1], v[2] );
    }
  
  this->m_SymmetryPlane.SetParameters( v );
  
  if ( this->m_SymmetryOutFileName )
    {
    ClassStreamOutput stream( this->m_SymmetryOutFileName, ClassStreamOutput::MODE_WRITE );
    stream << this->m_SymmetryPlane;
    stream.Close();
    }

  if ( this->m_AlignedOutFile ) 
    this->WriteAligned( originalVolume );

  if ( this->m_MarkedOutFile ) 
    this->WriteMarkPlane( originalVolume );

  if ( this->m_DifferenceOutFile )
    this->WriteDifference( originalVolume );

  if ( this->m_MirrorOutFile )
    WriteMirror( originalVolume );

  if ( this->m_WriteXformPath )
    {
    AffineXform::SmartPtr alignment( this->m_SymmetryPlane.GetAlignmentXform( 0 ) );
    XformIO::Write( alignment, this->m_WriteXformPath );
    }
  
  return 0;
}

bool
cmtk::ImageSymmetryPlaneCommandLineBase
::ParseCommandLine( const int argc, const char* argv[] )
{
  try
    {
    if ( ! this->m_CommandLine.Parse( argc, argv ) ) return false;
    
    if ( this->m_SymmetryParameters ) 
      {
      double rho, theta, phi;
      if ( 3 == sscanf( this->m_SymmetryParameters, "%20lf %20lf %20lf", &rho, &theta, &phi ) ) 
	{
	this->m_Rho = rho; 
	this->m_Theta = Units::Degrees( theta );
	this->m_Phi = Units::Degrees( phi );
	}
      }
    
    if ( this->m_SymmetryParametersFile ) 
      {
      ClassStreamInput inStream( this->m_SymmetryParametersFile );
      if ( inStream.IsValid() ) 
	{
	ParametricPlane *plane = NULL;
	inStream >> plane;
	this->m_Rho = plane->GetRho(); 
	this->m_Theta = plane->GetTheta();
	this->m_Phi = plane->GetPhi();
	delete plane;
	} 
      else
	{
	StdErr.printf( "ERROR: Could not open symmetry parameter file %s\n", this->m_SymmetryParametersFile );
	}
      }
    }
  catch ( const CommandLine::Exception& ex ) 
    {
    StdErr << ex << "\n";
    return false;
    }
  
  return true;
}

void
cmtk::ImageSymmetryPlaneCommandLineBase
::WriteDifference
( UniformVolume::SmartConstPtr& originalVolume ) const
{
  UniformVolume::SmartPtr diffVolume( originalVolume->CloneGrid() );
  const TypedArray* originalData = originalVolume->GetData();
  TypedArray::SmartPtr diffData = TypedArray::SmartPtr( TypedArray::Create( GetSignedDataType( originalData->GetType() ), originalData->GetDataSize() ) );
  diffVolume->SetData( diffData );

  Types::DataItem dataV, dataW;

  const UniformVolumeInterpolatorBase::SmartPtr interpolator( ReformatVolume::CreateInterpolator( this->m_Interpolation, originalVolume ) );;
  
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
	UniformVolume::CoordinateVectorType w = originalVolume->GetGridLocation( x, y, z );
	this->m_SymmetryPlane.MirrorInPlace( w );

	if ( interpolator->GetDataAt( w, dataW ) )
	  {
	  diffData->Set( fabs( dataV - dataW ), offset );
	  }
	else
	  {
	  diffData->SetPaddingAt( offset );
	  }
	}
  
  VolumeIO::Write( *diffVolume, this->m_DifferenceOutFile );
}

void
cmtk::ImageSymmetryPlaneCommandLineBase
::WriteMirror
( UniformVolume::SmartConstPtr& originalVolume ) const
{
  TypedArray::SmartPtr mirrorData = TypedArray::Create( originalVolume->GetData()->GetType(), originalVolume->GetData()->GetDataSize() );

  Types::DataItem data;

  const UniformVolumeInterpolatorBase::SmartPtr interpolator( ReformatVolume::CreateInterpolator( this->m_Interpolation, originalVolume ) );;
  
  int offset = 0;
  for ( int z = 0; z < originalVolume->GetDims()[2]; ++z ) 
    {
    for ( int y = 0; y < originalVolume->GetDims()[1]; ++y )
      for ( int x = 0; x < originalVolume->GetDims()[0]; ++x, ++offset ) 
	{
	UniformVolume::CoordinateVectorType v = originalVolume->GetGridLocation( x, y, z );
	this->m_SymmetryPlane.MirrorInPlace( v );

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

  UniformVolume::SmartPtr mirrorVolume( originalVolume->CloneGrid() );
  mirrorVolume->SetData( mirrorData );
  VolumeIO::Write( *mirrorVolume, this->m_MirrorOutFile );
}

void
cmtk::ImageSymmetryPlaneCommandLineBase
::WriteMarkPlane
( UniformVolume::SmartConstPtr& originalVolume ) const
{
  UniformVolume::SmartPtr markVolume( originalVolume->CloneGrid() );
  TypedArray::SmartPtr markData( originalVolume->GetData()->Clone() );
  markVolume->SetData( markData );

  int offset = 0;
  for ( int z = 0; z < originalVolume->GetDims()[2]; ++z ) 
    {
    for ( int y = 0; y < originalVolume->GetDims()[1]; ++y ) 
      {
      int currentSideOfPlane = 0;
      for ( int x = 0; x < originalVolume->GetDims()[0]; ++x, ++offset ) 
	{
	int newSideOfPlane = this->m_SymmetryPlane.GetWhichSide( originalVolume->GetGridLocation( x, y, z ) );
	if ( ( newSideOfPlane != currentSideOfPlane ) && x )
	  markData->Set( this->m_MarkPlaneValue, offset );
	currentSideOfPlane = newSideOfPlane;
	}
      }    
    }

  VolumeIO::Write( *markVolume, this->m_MarkedOutFile );
}

void
cmtk::ImageSymmetryPlaneCommandLineBase
::WriteAligned
( UniformVolume::SmartConstPtr& originalVolume ) const
{
  const TypedArray* originalData = originalVolume->GetData();

  TypedArray::SmartPtr alignData = TypedArray::Create( originalData->GetType(), originalData->GetDataSize() );
  if ( this->m_PadOutValueSet )
    {
    alignData->SetPaddingValue( this->m_PadOutValue );
    }

  UniformVolume::SmartPtr alignVolume( originalVolume->CloneGrid() );
  alignVolume->SetData( alignData );

  const Types::DataItem maxData = originalData->GetRange().m_UpperBound;

  Types::DataItem data;

  int normalAxis = 0;
  switch ( this->m_InitialPlane )
    {
    case SYMPL_INIT_XY: normalAxis = 2; break;
    case SYMPL_INIT_XZ: normalAxis = 1; break;
    case SYMPL_INIT_YZ: normalAxis = 0; break;
    }

  const UniformVolumeInterpolatorBase::SmartPtr interpolator( ReformatVolume::CreateInterpolator( this->m_Interpolation, originalVolume ) );;
  
  AffineXform::SmartPtr alignment( this->m_SymmetryPlane.GetAlignmentXform( normalAxis ) );
  int offset = 0;
  for ( int z = 0; z < originalVolume->GetDims()[2]; ++z ) 
    {
    for ( int y = 0; y < originalVolume->GetDims()[1]; ++y ) 
      {
      for ( int x = 0; x < originalVolume->GetDims()[0]; ++x, ++offset ) 
	{
	const UniformVolume::CoordinateVectorType v = alignment->Apply( originalVolume->GetGridLocation( x, y, z ) );

	if ( interpolator->GetDataAt( v, data ) )
	  {
	  if ( this->m_MarkPlaneAligned && (x == ( originalVolume->GetDims()[0] / 2 )) )
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

  VolumeIO::Write( *alignVolume, this->m_AlignedOutFile );
}

