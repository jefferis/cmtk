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

#include "cmtkImagePairAffineRegistrationCommandLine.h"

#include <System/cmtkConsole.h>
#include <System/cmtkThreads.h>
#include <System/cmtkTimers.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkCompressedStream.h>
#include <System/cmtkMountPoints.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkAnatomicalOrientation.h>
#include <Base/cmtkTransformChangeToSpaceAffine.h>
#include <Base/cmtkTransformChangeFromSpaceAffine.h>

#include <Registration/cmtkRegistrationCallback.h>
#include <Registration/cmtkProtocolCallback.h>
#include <Registration/cmtkMakeInitialAffineTransformation.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkClassStream.h>
#include <IO/cmtkClassStreamAffineXform.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkAffineXformITKIO.h>

#ifdef CMTK_USE_SQLITE
#  include <Registration/cmtkImageXformDB.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef HAVE_SYS_UTSNAME_H
#  include <sys/utsname.h>
#endif

#ifdef _MSC_VER
#  include <direct.h>
#endif

#include <stdio.h>
#include <string.h>
#include <iostream>

namespace 
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairAffineRegistrationCommandLine
::ImagePairAffineRegistrationCommandLine 
( const int argc, const char* argv[] ) 
  : m_InitialXformPath( NULL ),
    m_ReformattedImagePath( NULL ),
    m_OutputPathITK( NULL ),
#ifdef CMTK_USE_SQLITE
    m_UpdateDB( NULL ),
#endif
    m_ProtocolFileName( NULL )
{
  this->m_Metric = 0;

  this->m_AutoMultiLevels = 0;
  this->m_CoarsestResolution = -1;
  this->m_MaxStepSize = 8;
  this->m_MinStepSize = 0.1;
  this->m_Sampling = 1.0;
  OutParametersName = OutMatrixName = Studylist = Time = NULL;

  Verbose = 0;

  bool forceOutsideFlag = false;
  Types::DataItem forceOutsideValue = 0;

  const char* inStudylist = NULL;
  const char *InitialStudylist = NULL;
  Study1 = Study2 = NULL;

  const char* clArg1 = NULL; // input studylist or reference image
  const char* clArg2 = NULL; // empty or floating image

  try 
    {
    CommandLine cl( CommandLine::PROPS_XML );
    cl.SetProgramInfo( CommandLine::PRG_TITLE, "Rigid and affine registration" );
    cl.SetProgramInfo( CommandLine::PRG_DESCR, "This program performs rigid and affine image registration using multi-resolution optimization of voxel-based image similarity measures." );
    cl.SetProgramInfo( CommandLine::PRG_CATEG, "CMTK.Registration.Experimental" );

    typedef CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" )->SetProperties( CommandLine::PROPS_NOXML );
    cl.AddSwitch( Key( 'q', "quiet" ), &Verbose, false, "Quiet mode" )->SetProperties( CommandLine::PROPS_NOXML );

    cl.BeginGroup( "Automation", "Automation Options" );
    cl.AddOption( Key( "auto-multi-levels" ), &this->m_AutoMultiLevels, "Automatic optimization and resolution parameter generation for <n> levels" );

    cl.BeginGroup( "Optimization", "Optimization settings" );
    cl.AddOption( Key( "max-stepsize" ), &this->m_MaxStepSize, "Maximum optimizer step size, which determines search space exploration." );
    cl.AddOption( Key( "min-stepsize" ), &this->m_MinStepSize, "Minimum optimizer step size, which determines precision." );
    cl.AddOption( Key( "stepfactor" ), &this->m_OptimizerStepFactor, "Factor for search step size reduction. Must be > 0.0 and < 1.0" );
    cl.AddOption( Key( "delta-f-threshold" ), &this->m_DeltaFThreshold, "Optional threshold to terminate optimization (level) if relative change of target function drops below this value." );
    cl.EndGroup();

    cl.BeginGroup( "Resolution", "Image resolution parameters" );
    cl.AddOption( Key( 's', "sampling" ), &this->m_Sampling, "Image sampling (finest resampled image resolution)" );
    cl.AddOption( Key( "coarsest" ), &this->m_CoarsestResolution, "Upper limit for image sampling in multiresolution hierarchy" );

    cl.AddSwitch( Key( "omit-original-data" ), &this->m_UseOriginalData, false, "Do not use original data in full resolution, omit final stage in multiresolution hierarchy, thus reducing computation time." );
    cl.EndGroup();

    cl.BeginGroup( "Transformation", "Transformation parameters" );
    cl.AddVector( Key( "dofs" ), this->NumberDOFs, "Add number of degrees of freedom. This can be 3 (translation), 6 (rigid: translation and rotation), "
		  "7 (rigid plus global scale), 9 (rigid plus anisotropic scales), 12 (rigid plus scales plus shears), or 603 (rigid plus shears, but no scale). "
		  "This option can be repeated, in which case DOFs are used for successive optimization runs in the order that they appear." );
    cl.AddVector( Key( "dofs-final" ), this->NumberDOFsFinal, "Add number of degrees of freedom for final level only [can be repeated]" );
    
    CommandLine::EnumGroup<MakeInitialAffineTransformation::Mode>::SmartPtr
      initGroup = cl.AddEnum( "init", &this->m_Initializer, "Select initializer for the affine trasnformation." );
    initGroup->AddSwitch( Key( "none" ), MakeInitialAffineTransformation::NONE, "Use input transformation, or identity transformation if none was provided ." );
    initGroup->AddSwitch( Key( "fov" ), MakeInitialAffineTransformation::FOV, "Align centers of field of view (or crop regions) using a translation." );
    initGroup->AddSwitch( Key( "com" ), MakeInitialAffineTransformation::COM, "Align centers of mass using a translation." );
    initGroup->AddSwitch( Key( "pax" ), MakeInitialAffineTransformation::PAX, "Align images by rotation using principal axes and translation using centers of mass." );
    initGroup->AddSwitch( Key( "physical" ), MakeInitialAffineTransformation::PHYS, "Align images by rotation using direction vectors stored in input images and translation using image origins." );
    
    cl.AddOption( Key( "initial" ), &InitialStudylist, "Initialize transformation from given path" )->SetProperties( CommandLine::PROPS_XFORM );
    cl.AddSwitch( Key( "initial-is-inverse" ), &this->m_InitialXformIsInverse, true, "Invert initial transformation before initializing registration" );
    cl.EndGroup();

    cl.BeginGroup( "Image data", "Image data" );
    CommandLine::EnumGroup<int>::SmartPtr
      metricGroup = cl.AddEnum( "registration-metric", &this->m_Metric, "Registration metric for motion estimation by image-to-image registration." );
    metricGroup->AddSwitch( Key( "nmi" ), 0, "Normalized Mutual Information metric" );
    metricGroup->AddSwitch( Key( "mi" ), 1, "Standard Mutual Information metric" );
    metricGroup->AddSwitch( Key( "cr" ), 2, "Correlation Ratio metric" );
    metricGroup->AddSwitch( Key( "rms" ), 3, "Root of Mean Squaresa metric (this is the square root of MSD)" );
    metricGroup->AddSwitch( Key( "msd" ), 4, "Mean Squared Difference metric" );
    metricGroup->AddSwitch( Key( "ncc" ), 5, "Normalized Cross Correlation metric" );

    cl.BeginGroup( "Interpolation", "Floating Image Interpolation Options" );
    cmtk::CommandLine::EnumGroup<Interpolators::InterpolationEnum>::SmartPtr kernelGroup = 
      cl.AddEnum( "interpolation", &this->m_FloatingImageInterpolation, "Interpolation method for floating image sampling:" );
    kernelGroup->AddSwitch( Key( "nearest-neighbor" ), Interpolators::NEAREST_NEIGHBOR, "Nearest neighbor interpolation (for intensity and label data)" );
    kernelGroup->AddSwitch( Key( "linear" ), Interpolators::LINEAR, "Trilinear interpolation" );
    kernelGroup->AddSwitch( Key( "cubic" ), Interpolators::CUBIC, "Tricubic interpolation" );
    kernelGroup->AddSwitch( Key( "cosine-sinc" ), Interpolators::COSINE_SINC, "Cosine-windowed sinc interpolation (most accurate but slowest)" );
    kernelGroup->AddSwitch( Key( "partial-volume" ), Interpolators::PARTIALVOLUME, "Partial volume interpolation (for label data)" );

    cl.AddSwitch( Key( "match-histograms" ), &this->m_MatchFltToRefHistogram, true, "Match floating image histogram to reference image histogram." );
    cl.AddOption( Key( "force-outside-value" ), &forceOutsideValue, "Force values outside field of view to this value rather than drop incomplete pixel pairs", &forceOutsideFlag );

    this->m_PreprocessorRef.AttachToCommandLine( cl );
    this->m_PreprocessorFlt.AttachToCommandLine( cl );

    cl.BeginGroup( "Output", "Output parameters" )->SetProperties( CommandLine::PROPS_NOXML );
    cl.AddOption( Key( 'o', "output" ), &this->Studylist, "Output path for final transformation" );
    cl.AddOption( Key( "write-matrix" ), &this->OutMatrixName, "Output path for final transformation in matrix format" );
    cl.AddOption( Key( "write-parameters" ), &this->OutParametersName, "Output path for final transformation in plain parameter list format" );
    cl.AddOption( Key( "write-protocol" ), &this->m_ProtocolFileName, "Optimization protocol output file name" );
    cl.AddOption( Key( "write-time" ), &this->Time, "Computation time statistics output file name" );
    cl.EndGroup();

    cl.BeginGroup( "SlicerImport", "Import Results into Slicer" );
    cl.AddOption( Key( "write-itk" ), &this->m_OutputPathITK, "Output path for final transformation in ITK format" )
      ->SetProperties( CommandLine::PROPS_XFORM | CommandLine::PROPS_OUTPUT )
      ->SetAttribute( "reference", "FloatingImage" );
    cl.AddOption( Key( "write-reformatted" ), &this->m_ReformattedImagePath, "Write reformatted floating image." )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OUTPUT );
    cl.EndGroup();
    
#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &this->m_UpdateDB, "Path to image/transformation database that should be updated with the new registration and/or reformatted image." );
    cl.EndGroup();
#endif

    cl.AddParameter( &clArg1, "ReferenceImage", "Reference (fixed) image path" )->SetProperties( CommandLine::PROPS_IMAGE );
    cl.AddParameter( &clArg2, "FloatingImage", "Floating (moving) image path" )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OPTIONAL );

    cl.Parse( argc, argv );
    }
  catch ( const CommandLine::Exception& ex )
    {
    StdErr << ex << "\n";
    exit( 1 );
    }
  
  if ( (this->m_OptimizerStepFactor <= 0) || (this->m_OptimizerStepFactor >= 1) ) 
    {
    StdErr << "ERROR: step factor value " << this->m_OptimizerStepFactor << " is invalid. Must be in range (0..1)\n";
    exit( 1 );
    }

  if ( clArg2 ) 
    {
    AffineXform::SmartPtr initialXform( new AffineXform() );
    this->SetInitialTransformation( initialXform );
    
    Study1 = const_cast<char*>( clArg1 );
    Study2 = const_cast<char*>( clArg2 );
    } 
  else
    {
    inStudylist = clArg1;
    AffineXform::SmartPtr initialXform( new AffineXform() );
    this->SetInitialTransformation( initialXform );

    if ( InitialStudylist ) 
      {
      StdErr << "Transformation will be overriden by '--initial' list.\n";
      }
    
    if ( Verbose )
      StdErr << "Reading input studylist " << inStudylist << ".\n";
    
    ClassStream typedStream( MountPoints::Translate(inStudylist), "registration", ClassStream::READ );
    if ( ! typedStream.IsValid() ) 
      {
      StdErr << "Could not open studylist archive " << inStudylist << ".\n";
      exit( 1 );
      }

    typedStream.Seek ( "registration" );
    Study1 = typedStream.ReadString( "reference_study" );
    Study2 = typedStream.ReadString( "floating_study" );
    if ( Study2 )
      {
      AffineXform::SmartPtr affineXform;
      typedStream >> affineXform;
      this->SetInitialTransformation( affineXform );
      }
    else
      {
      // legacy studylists have inverse transformation in them
      Study2 = typedStream.ReadString( "model_study" );
      AffineXform::SmartPtr affineXform;
      typedStream >> affineXform;
      this->SetInitialTransformation( affineXform->GetInverse() );
      }

    typedStream.Close();
    }

  if ( !Study1 )
    {
    StdErr << "ERROR: reference image path resolved to NULL.\n";
    exit( 1 );
    }
  
  if ( !Study2 )
    {
    StdErr << "ERROR: floating image path resolved to NULL.\n";
    exit( 1 );
    }
  
  UniformVolume::SmartPtr volume( VolumeIO::ReadOriented( Study1, Verbose ) );
  if ( !volume )
    {
    StdErr << "ERROR: volume " << this->Study1 << " could not be read\n";
    exit( 1 );
    }
  this->SetVolume_1( UniformVolume::SmartPtr( this->m_PreprocessorRef.GetProcessedImage( volume ) ) );

  volume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( Study2, Verbose ) );
  if ( !volume )
    {
    StdErr << "ERROR: volume " << this->Study2 << " could not be read\n";
    exit( 1 );
    }
  this->SetVolume_2(  UniformVolume::SmartPtr( this->m_PreprocessorFlt.GetProcessedImage( volume ) ) );

  if ( InitialStudylist ) 
    {
    Xform::SmartPtr xform( XformIO::Read( InitialStudylist, Verbose ) );
    if ( ! xform ) 
      {
      StdErr << "ERROR: could not read transformation from " << InitialStudylist << "\n";
      exit( 1 );
      }
    
    AffineXform::SmartPtr affine( AffineXform::SmartPtr::DynamicCastFrom( xform ) );
    if ( ! affine )
      {
      StdErr << "ERROR: transformation " << InitialStudylist << " is not affine.\n";
      exit( 1 );
      }

    if ( affine->m_MetaInformation[META_SPACE] != AnatomicalOrientation::ORIENTATION_STANDARD )
      {
      TransformChangeFromSpaceAffine toStandardSpace( *affine, *(this->m_Volume_1), *(this->m_Volume_2) );
      *affine = toStandardSpace.GetTransformation();
      affine->m_MetaInformation[META_SPACE] = AnatomicalOrientation::ORIENTATION_STANDARD;
      }
    
    this->SetInitialTransformation( affine );
    }
  
  if ( this->m_Initializer != MakeInitialAffineTransformation::NONE ) 
    {
    if ( inStudylist || InitialStudylist ) 
      {
      StdErr << "WARNING: Initial transformation was taken from studylist. Selected transformation initializer will be ignored.\n";
      } 
    }
  
  if ( this->m_ProtocolFileName ) 
    {
    RegistrationCallback::SmartPtr callback( new ProtocolCallback( this->m_ProtocolFileName ) );
    this->SetCallback( callback );
    }

  if ( forceOutsideFlag )
    {
    this->SetForceOutside( true, forceOutsideValue );
    }
}

CallbackResult
ImagePairAffineRegistrationCommandLine::InitRegistration ()
{
  CallbackResult Result = Superclass::InitRegistration();
  return Result;
}
	
void
ImagePairAffineRegistrationCommandLine::OutputResultMatrix( const char* matrixName ) const
{
  Types::Coordinate matrix[4][4];
  this->GetTransformation()->GetMatrix( matrix );

  FILE* mfile = fopen( matrixName, "w" );
  if ( mfile )
    {
    for ( int i = 0; i < 4; ++i )
      {
      fprintf( mfile, "%e\t%e\t%e\t%e\n", matrix[0][i], matrix[1][i], matrix[2][i], matrix[3][i] );
      }
    fclose( mfile );
    }
}

void
ImagePairAffineRegistrationCommandLine::OutputResultParameters
( const char* paramsName, const CoordinateVector& v ) const
{
  FILE* pfile = fopen( paramsName, "w" );
  if ( pfile )
    {
    for ( unsigned int idx=0; idx < v.Dim; ++idx )
      fprintf( pfile, "#%d: %f\n", idx, v.Elements[idx] );
    fclose( pfile );
    }
}

void
ImagePairAffineRegistrationCommandLine::OutputResultList( const char* studyList ) const
{
  ClassStream classStream( studyList, "studylist", ClassStream::WRITE );
  if ( !classStream.IsValid() ) return;
  
  classStream.Begin( "studylist" );
  classStream.WriteInt( "num_sources", 2 );
  classStream.End();
    
  classStream.Begin( "source" );
  classStream.WriteString( "studyname", CompressedStream::GetBaseName( Study1 ) );
  classStream.End();
    
  classStream.Begin( "source" );
  classStream.WriteString( "studyname", CompressedStream::GetBaseName( Study2 ) );
  classStream.End();
    
  classStream.Close();
    
  classStream.Open( studyList, "registration", ClassStream::WRITE );
    
  classStream.Begin( "registration" );
  classStream.WriteString( "reference_study", CompressedStream::GetBaseName( Study1 ) );
  classStream.WriteString( "floating_study", CompressedStream::GetBaseName( Study2 ) );
    
  classStream << *(this->GetTransformation());
    
  classStream.End();
  classStream.Close();
    
  classStream.Open( studyList, "settings", ClassStream::WRITE );
  classStream.WriteDouble( "exploration", this->m_MaxStepSize );
  classStream.WriteDouble( "accuracy", this->m_MinStepSize );
  classStream.WriteDouble( "min_sampling", this->m_Sampling );
  classStream.WriteDouble( "coarsest_resolution", this->m_CoarsestResolution );
  classStream.WriteInt( "metric", this->m_Metric );
  classStream.WriteDouble( "optimizer_step_factor", this->m_OptimizerStepFactor );
  classStream.WriteString( "initializer", MakeInitialAffineTransformation::GetModeName( this->m_Initializer ) );

  this->m_PreprocessorRef.WriteSettings( classStream );  
  this->m_PreprocessorFlt.WriteSettings( classStream );  

  classStream.Close();
    
  classStream.Open( studyList, "statistics", ClassStream::WRITE );
  classStream.WriteDouble( "time", this->GetTotalElapsedTime() );
  classStream.WriteDouble( "walltime", this->GetTotalElapsedWalltime() );
#ifdef CMTK_USE_THREADS
  classStream.WriteDouble( "thread_time", this->GetThreadTotalElapsedTime() );
#endif
    
#ifndef _MSC_VER
  struct utsname name;
  if ( uname( &name ) >= 0 ) 
    {
    classStream.WriteString( "host", name.nodename );
    classStream.WriteString( "system", name.sysname );
    }
#endif
  classStream.Close();
}

void
ImagePairAffineRegistrationCommandLine::OutputResult ( const CoordinateVector* v )
{
  if ( Verbose ) 
    {
    fprintf( stderr, "\rResulting transformation parameters: \n" );
    for ( unsigned int idx=0; idx<v->Dim; ++idx )
      fprintf( stderr, "#%d: %f\n", idx, v->Elements[idx] );
    }
  
  if ( this->OutMatrixName )
    {
    this->OutputResultMatrix( this->OutMatrixName );
    }

  if ( this->OutParametersName )
    {
    this->OutputResultParameters( this->OutParametersName, *v );
    }

  if ( this->Studylist ) 
    {
    this->OutputResultList( this->Studylist );
    }

  if ( this->m_OutputPathITK ) 
    {
    TransformChangeToSpaceAffine toNative( *(this->GetTransformation()), *(this->m_Volume_1), *(this->m_Volume_2) );
    AffineXformITKIO::Write( this->m_OutputPathITK, toNative.GetTransformation() );
    }

  if ( this->m_ReformattedImagePath )
    {
    VolumeIO::Write( *(this->GetReformattedFloatingImage()), this->m_ReformattedImagePath, this->Verbose );
    }

#ifdef CMTK_USE_SQLITE
  if ( this->m_UpdateDB )
    {
    try
      {
      cmtk::ImageXformDB db( this->m_UpdateDB );
      
      if ( this->m_ReformattedImagePath )
	{
	db.AddImage( this->m_ReformattedImagePath, this->m_ReferenceVolume->GetMetaInfo( META_FS_PATH ) );
	}
      
      if ( this->Studylist )
	{
	if ( this->m_InitialXformPath ) 
	  {
	  db.AddRefinedXform( this->Studylist, true /*invertible*/, this->m_InitialXformPath, this->m_InitialXformIsInverse );
	  }
	else
	  {
	  db.AddImagePairXform( this->Studylist, true /*invertible*/, this->m_ReferenceVolume->GetMetaInfo( META_FS_PATH ), this->m_FloatingVolume->GetMetaInfo( META_FS_PATH ) );
	  }
	}
      }
    catch ( const cmtk::ImageXformDB::Exception& ex )
      {
      StdErr << "DB ERROR: " << ex.what() << " on database " << this->m_UpdateDB << "\n";
      }
    }
#endif
}

void
ImagePairAffineRegistrationCommandLine::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f,
  const int index, const int total )
{
  if ( Verbose )
    fprintf( stderr, "\rEntering resolution level %d out of %d...\n", index, total );
  this->Superclass::EnterResolution( v, f, index, total );
}

CallbackResult
ImagePairAffineRegistrationCommandLine::Register ()
{
  const double baselineTime = Timers::GetTimeProcess();
  CallbackResult Result = Superclass::Register();
  const int elapsed = static_cast<int>( Timers::GetTimeProcess() - baselineTime );

  if ( Time ) 
    {
    FILE *tfp = fopen( Time, "w" );
    
    if ( tfp ) 
      {
      fprintf( tfp, "%d\n", elapsed );
      fclose( tfp );
      } 
    else 
      {
      std::cerr << "Could not open time file " << Time << "\n";
      }
    }
  return Result;
}

} // namespace cmtk

