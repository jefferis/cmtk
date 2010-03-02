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

#include <cmtkImagePairNonrigidRegistrationCommandLine.h>

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>
#include <cmtkMountPoints.h>
#include <cmtkTimers.h>
#include <cmtkThreads.h>

#include <cmtkXformIO.h>
#include <cmtkClassStream.h>
#include <cmtkClassStreamAffineXform.h>
#include <cmtkCompressedStream.h>
#include <cmtkTypedStreamStudylist.h>
#include <cmtkVolumeIO.h>
#include <cmtkAnatomicalOrientation.h>
#include <cmtkTransformChangeFromSpaceAffine.h>

#include <cmtkVoxelMatchingElasticFunctional.h>
#include <cmtkSymmetricElasticFunctional.h>

#include <cmtkSplineWarpXformITKIO.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef HAVE_SYS_UTSNAME_H
#  include <sys/utsname.h>
#endif

#ifdef HAVE_SIGNAL_H
#  include <signal.h>
#endif

#ifdef _MSC_VER
#  include <direct.h>
#endif // #ifdef _MSC_VER

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairNonrigidRegistrationCommandLine* ImagePairNonrigidRegistrationCommandLine::StaticThis = NULL;

ImagePairNonrigidRegistrationCommandLine
::ImagePairNonrigidRegistrationCommandLine
( int argc, char *argv[] ) :
  m_OutputPathITK( NULL ),
  m_ReformattedImagePath( NULL )
{
  Studylist = Time = NULL;

  this->m_OutputIntermediate = 0;
  Verbose = 0;

  InputStudylist = NULL;
  const char *initialTransformationFile = NULL;
  bool initialTransformationInverse = false;
  IntermediateResultIndex = 0;

  bool forceOutsideFlag = false;
  Types::DataItem forceOutsideValue = 0;

  const char* clArg1 = NULL; // input studylist or reference image
  const char* clArg2 = NULL; // empty or floating image
  const char* clArg3 = NULL; // empty or initial transformation

  try
    {
    CommandLine cl( argc, argv, CommandLine::PROPS_XML );
    cl.SetProgramInfo( CommandLine::PRG_TITLE, "B-spline nonrigid registration" );
    cl.SetProgramInfo( CommandLine::PRG_DESCR, "This program performs nonrigid image registration using multi-resolution optimization of voxel-based image similarity measures "
		       "and a multi-resolution B-spline transformation model." );
    cl.SetProgramInfo( CommandLine::PRG_CATEG, "CMTK.Registration" );

    typedef CommandLine::Key Key;
    cl.BeginGroup( "Console", "Console output control" )->SetProperties( CommandLine::PROPS_NOXML );
    cl.AddSwitch( Key( 'v', "verbose" ), &this->Verbose, true, "Verbose peration" );
    cl.AddSwitch( Key( 'q', "quiet" ), &Verbose, false, "Quiet mode" );
    cl.EndGroup();

    cl.BeginGroup( "TransformationIO", "Transformation import/export" );
    cl.AddOption( Key( "initial" ), &initialTransformationFile, "Initialize transformation from given path" )->SetProperties( CommandLine::PROPS_XFORM );
    cl.AddSwitch( Key( "invert-initial" ), &initialTransformationInverse, true, "Invert given (affine) initial transformation." );
    cl.AddOption( Key( "write-itk-xform" ), &this->m_OutputPathITK, "Output path for final transformation in ITK format" )
      ->SetProperties( CommandLine::PROPS_XFORM | CommandLine::PROPS_OUTPUT )
      ->SetAttribute( "type", "bspline" )->SetAttribute( "reference", "FloatingImage" );
    cl.AddOption( Key( "write-reformatted" ), &this->m_ReformattedImagePath, "Write reformatted floating image." )
      ->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OUTPUT );
    cl.EndGroup();

    cl.BeginGroup( "Transformation", "Transformation parameters" );
    cl.AddOption( Key( 'g', "grid-spacing" ), &this->m_GridSpacing, "Control point grid spacing" );
    cl.AddSwitch( Key( "exact-spacing" ), &this->m_ExactGridSpacing, true, "Use exact control point spacing; do not modify spacing to fit reference image bounding box" );
    cl.AddOption( Key( "refine" ), &this->m_RefineGrid, "Number of refinements (control point grid resolution levels)" );
    cl.AddSwitch( Key( "delay-refine" ), &this->m_DelayRefineGrid, true, "Delay control point grid refinement; first switch to next higher image resolution" );

    cl.AddOption( Key( "ignore-edge" ), &this->IgnoreEdge, "Ignore n control point layers along each image face" );
    cl.AddOption( Key( "restrict" ), &this->RestrictToAxes, "Restrict deformation to coordinate dimension(s) [one or more of 'x','y','z']" );

    cl.AddSwitch( Key( "no-adaptive-fix" ), &this->m_AdaptiveFixParameters, false, "Disable adaptive fixing of control points; optimize all deformation parameters" );
    cl.AddOption( Key( "adaptive-fix-thresh" ), &this->m_AdaptiveFixThreshFactor, "Threshold factor for entropy criterion to fix local control points" );
    cl.AddSwitch( Key( "accurate" ), &this->m_FastMode, false, "Accurate computation mode: may give slightly better results after substantially longer computation" );
    cl.AddSwitch( Key( "fast" ), &this->m_FastMode, true, "Fast computation mode: may give slightly worse results than accurate mode, but saves substantial CPU time" );
    cl.EndGroup();

    cl.BeginGroup( "Optimization", "Optimization parameters" );
    cl.AddOption( Key( "max-stepsize" ), &this->m_MaxStepSize, "Maximum optimizer step size, which determines search space exploration." );
    cl.AddOption( Key( "min-stepsize" ), &this->m_MinStepSize, "Minimum optimizer step size, which determines precision." );
    cl.AddOption( Key( "stepfactor" ), &this->m_OptimizerStepFactor, "Factor for search step size reduction. Must be > 0.0 and < 1.0 [default: 0.5]" );
    cl.AddOption( Key( "delta-f-threshold" ), &this->m_DeltaFThreshold, "Optional threshold to terminate optimization (level) if relative change of target function drops below this value." );

    cl.AddSwitch( Key( "no-maxnorm" ), &this->m_UseMaxNorm, false, "Use Euclid norm for gradient normalication in optimization, rather than maximum norm" );

    cl.AddOption( Key( "jacobian-constraint-weight" ), &this->m_JacobianConstraintWeight, "Weight for Jacobian-based local volume preservation constraint" );
    cl.AddOption( Key( "smoothness-constraint-weight" ), &this->m_GridEnergyWeight, "Weight for smoothness constraint based on second-order grid bending energy." );
    cl.AddOption( Key( "landmark-constraint-weight" ), &this->m_LandmarkErrorWeight, "Weight for landmark misregistration registration" );
    cl.AddOption( Key( "inverse-consistency-weight" ), &this->m_InverseConsistencyWeight, "Weight for inverse consistency constraint" );
    cl.AddOption( Key( "constraint-relaxation-factor" ), &this->m_RelaxWeight, "Weight relaxation factor for alternating under-constrained iterations" );
    cl.EndGroup();

    cl.BeginGroup( "Resolution", "Image resolution parameters" );
    cl.AddOption( Key( 's', "sampling" ), &this->m_Sampling, "Image sampling (finest resampled image resolution)" );
    cl.AddOption( Key( "coarsest" ), &this->m_CoarsestResolution, "Upper limit for image sampling in multiresolution hierarchy" );

    cl.AddSwitch( Key( "omit-original-data" ), &this->m_UseOriginalData, false, "Do not use original data in full resolution for final registration stage." );
    cl.EndGroup();

    cl.BeginGroup( "Images", "Image data" );
    CommandLine::EnumGroup<int>::SmartPtr
      metricGroup = cl.AddEnum( "registration-metric", &this->m_Metric, "Registration metric for motion estimation by image-to-image registration." );
    metricGroup->AddSwitch( Key( "nmi" ), 0, "Normalized Mutual Information metric" );
    metricGroup->AddSwitch( Key( "mi" ), 1, "Standard Mutual Information metric" );
    metricGroup->AddSwitch( Key( "cr" ), 2, "Correlation Ratio metric" );
    metricGroup->AddSwitch( Key( "msd" ), 4, "Mean Squared Difference metric" );
    metricGroup->AddSwitch( Key( "ncc" ), 5, "Normalized Cross Correlation metric" );

    cl.AddSwitch( Key( "match-histograms" ), &this->m_MatchFltToRefHistogram, true, "Match floating image histogram to reference image histogram." );
    cl.AddSwitch( Key( "repeat-match-histograms" ), &this->m_RepeatMatchFltToRefHistogram, true, "Repeat histogram matching after every level of the registration to account for volume changes." );
    cl.AddOption( Key( "force-outside-value" ), &forceOutsideValue, "Force values outside field of view to this value rather than drop incomplete pixel pairs", &forceOutsideFlag );

    this->m_PreprocessorRef.AttachToCommandLine( cl );
    this->m_PreprocessorFlt.AttachToCommandLine( cl );

    cl.BeginGroup( "Output", "Output parameters" )->SetProperties( CommandLine::PROPS_NOXML );
    cl.AddOption( Key( 'o', "outlist" ), &this->Studylist, "Output path for final transformation" );
    cl.AddOption( Key( 't', "time" ), &this->Time, "Computation time statistics output file name" );
    cl.AddSwitch( Key( "output-intermediate" ), &this->m_OutputIntermediate, true, "Write transformation for each level [default: only write final transformation]" );
    cl.EndGroup();

    cl.AddParameter( &clArg1, "ReferenceImage", "Reference (fixed) image path" )
      ->SetProperties( CommandLine::PROPS_IMAGE );
    cl.AddParameter( &clArg2, "FloatingImage", "Floating (moving) image path" )
      ->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OPTIONAL);
    cl.AddParameter( &clArg3, "InitialXform", "Initial affine transformation from reference to floating image" )
      ->SetProperties( CommandLine::PROPS_NOXML | CommandLine::PROPS_XFORM | CommandLine::PROPS_OPTIONAL );

    cl.Parse();
    }
  catch ( CommandLine::Exception ex )
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
    this->SetInitialTransformation( AffineXform::SmartPtr( new AffineXform() ) );
    
    Study1 = const_cast<char*>( clArg1 );
    Study2 = const_cast<char*>( clArg2 );

    if ( clArg3 )
      {
      initialTransformationFile = clArg3;
      }
    }
  else 
    {
    InputStudylist = clArg1;
    
    if ( Verbose )
      fprintf( stderr, "Reading input studylist %s.\n", InputStudylist );
    
    ClassStream classStream( MountPoints::Translate(InputStudylist),"registration", ClassStream::READ );
    if ( ! classStream.IsValid() ) 
      {
      std::cerr << "Could not open studylist archive " << InputStudylist << ".\n";
      exit( 1 );
      }
    
    classStream.Seek ( "registration" );
    Study1 = classStream.ReadString( "reference_study" );
    Study2 = classStream.ReadString( "floating_study" );
    if ( Study2 )
      {
      AffineXform::SmartPtr affineXform;
      classStream >> affineXform;
      this->SetInitialTransformation( affineXform );
      }
    else
      {
      // legacy studylists have inverse transformation stored in them
      Study2 = classStream.ReadString( "model_study" );
      AffineXform::SmartPtr affineXform;
      classStream >> affineXform;
      this->SetInitialTransformation( affineXform->GetInverse() );
      }
    
    classStream.Close();
    }

  /// Was an initial studylist given? If so, get warp 
  if ( initialTransformationFile ) 
    {
    Xform::SmartPtr initialXform( XformIO::Read( initialTransformationFile ) );
    AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( initialXform );
    if ( affineXform )
      {
      if ( initialTransformationInverse )
	this->SetInitialTransformation( affineXform->GetInverse() );
      else
	this->SetInitialTransformation( affineXform );
      }
    else
      {
      InitialWarpXform = SplineWarpXform::SmartPtr::DynamicCastFrom( initialXform );
      }
    }

  UniformVolume::SmartPtr volume( VolumeIO::ReadOriented( Study1, Verbose ) );
  if ( !volume ) throw ConstructorFailed();
  this->SetVolume_1( UniformVolume::SmartPtr( this->m_PreprocessorRef.GetProcessedImage( volume ) ) );

  volume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( Study2, Verbose ) );
  if ( !volume ) throw ConstructorFailed();
  this->SetVolume_2( UniformVolume::SmartPtr( this->m_PreprocessorFlt.GetProcessedImage( volume ) ) );

  AffineXform::SmartPtr affineXform( AffineXform::SmartPtr::DynamicCastFrom( this->m_InitialTransformation ) );
  if ( affineXform )
    {
    if ( affineXform->m_MetaInformation[CMTK_META_SPACE] != AnatomicalOrientation::ORIENTATION_STANDARD )
      {
      TransformChangeFromSpaceAffine toStandardSpace( *affineXform, *(this->m_Volume_1), *(this->m_Volume_2) );
      *affineXform = toStandardSpace.GetTransformation();
      affineXform->m_MetaInformation[CMTK_META_SPACE] = AnatomicalOrientation::ORIENTATION_STANDARD;
      this->SetInitialTransformation( affineXform );
      }
    }

  if ( forceOutsideFlag )
    {
    this->SetForceOutside( true, forceOutsideValue );
    }
}

ImagePairNonrigidRegistrationCommandLine::~ImagePairNonrigidRegistrationCommandLine ()
{
#ifdef HAVE_SIGRELSE
  // release signal handler
  sigrelse( SIGUSR1 );
#endif
}

CallbackResult
ImagePairNonrigidRegistrationCommandLine
::InitRegistration ()
{
  CallbackResult result = this->Superclass::InitRegistration();
  if ( result != CALLBACK_OK )
    return result;

  if ( this->m_OutputIntermediate )
    this->OutputIntermediate();
  
  // register signal handler for intermediate result output.
  Self::StaticThis = this;
#ifndef _MSC_VER
  signal( SIGUSR1, cmtkImagePairNonrigidRegistrationCommandLineDispatchSIGUSR1 );
#endif

  return CALLBACK_OK;
}
	
void
ImagePairNonrigidRegistrationCommandLine
::OutputResult
( const CoordinateVector* )
{
  if ( Studylist ) 
    this->OutputWarp( Studylist );

  if ( this->m_OutputPathITK ) 
    {
    SplineWarpXformITKIO::Write( this->m_OutputPathITK, *(this->GetTransformation()), *(this->m_ReferenceVolume), *(this->m_FloatingVolume) );
    }

  if ( this->m_ReformattedImagePath )
    {
    VolumeIO::Write( UniformVolume::SmartPtr( this->GetReformattedFloatingImage() ), this->m_ReformattedImagePath, this->Verbose );
    }
}

void
ImagePairNonrigidRegistrationCommandLine
::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, 
  const int index, const int total )
{
  if ( Verbose )
    fprintf( stderr, "\rEntering resolution level %d out of %d...\n", index, total );
  
  this->Superclass::EnterResolution( v, f, index, total );
}

int
ImagePairNonrigidRegistrationCommandLine::DoneResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, const int index, const int total )
{
  if ( this->m_OutputIntermediate )
    this->OutputIntermediate();
  return this->Superclass::DoneResolution( v, f, index, total );
}

void
ImagePairNonrigidRegistrationCommandLine::OutputWarp ( const char* path ) const
{
  ClassStream classStream( path, "studylist", ClassStream::WRITE );
  if ( ! classStream.IsValid() ) return;

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

  classStream.Open( path, "settings", ClassStream::WRITE );
  classStream.WriteInt( "algorithm", this->m_Algorithm );
  classStream.WriteBool( "use_maxnorm", this->m_UseMaxNorm );
  classStream.WriteDouble( "exploration", this->m_MaxStepSize );
  classStream.WriteDouble( "accuracy", this->m_MinStepSize );
  classStream.WriteDouble( "min_sampling", this->m_Sampling );
  classStream.WriteDouble( "coarsest_resolution", this->m_CoarsestResolution );
  classStream.WriteBool( "use_original_data", this->m_UseOriginalData );
  classStream.WriteBool( "fast_mode", this->m_FastMode );
  classStream.WriteInt( "metric", this->m_Metric );
  classStream.WriteDouble( "optimizer_step_factor", this->m_OptimizerStepFactor );
  classStream.WriteDouble( "grid_spacing", this->m_GridSpacing );
  classStream.WriteInt( "ignore_edge", IgnoreEdge );
  classStream.WriteDouble( "jacobian_constraint_weight", this->m_JacobianConstraintWeight );
  classStream.WriteDouble( "energy_constraint_weight", this->m_GridEnergyWeight );
  classStream.WriteDouble( "inverse_consistency_weight", this->m_InverseConsistencyWeight );
  classStream.WriteDouble( "weight_relaxation", this->m_RelaxWeight );
  classStream.WriteDouble( "landmark_error_weight", this->m_LandmarkErrorWeight );
  classStream.WriteInt( "refine_grid", this->m_RefineGrid );
  classStream.WriteBool( "delay_refine_grid", this->m_DelayRefineGrid );
  classStream.WriteBool( "adaptive_fix_parameters", this->m_AdaptiveFixParameters );
  classStream.WriteDouble( "adaptive_fix_parameters_thresh", this->m_AdaptiveFixThreshFactor );

  this->m_PreprocessorRef.WriteSettings( classStream );  
  this->m_PreprocessorFlt.WriteSettings( classStream );  

  classStream.Close();
      
  classStream.Open( path, "statistics", ClassStream::WRITE );
  classStream.WriteDouble( "time_level", this->GetLevelElapsedTime() );
  classStream.WriteDouble( "time_total", this->GetTotalElapsedTime() );
  classStream.WriteDouble( "walltime_level", this->GetLevelElapsedWalltime() );
  classStream.WriteDouble( "walltime_total", this->GetTotalElapsedWalltime() );
  classStream.WriteDouble( "thread_time_level", this->GetThreadLevelElapsedTime() );
  classStream.WriteDouble( "thread_time_total", this->GetThreadTotalElapsedTime() );
  classStream.WriteInt( "number_of_threads", Threads::NumberOfThreads );
  classStream.WriteInt( "number_of_cpus", Threads::GetNumberOfProcessors() );

#ifndef _MSC_VER
  struct utsname name;
  if ( uname( &name ) >= 0 ) 
    {
    classStream.WriteString( "host", name.nodename );
    classStream.WriteString( "system", name.sysname );
    }
#endif
  classStream.Close();

  const WarpXform::SmartPtr warp = WarpXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  if ( warp ) 
    {
    classStream.Open( path, "registration", ClassStream::WRITE );
    if ( classStream.IsValid() ) 
      {
      classStream.Begin( "registration" );
      classStream.WriteString( "reference_study", CompressedStream::GetBaseName( this->Study1 ) );
      classStream.WriteString( "floating_study", CompressedStream::GetBaseName( this->Study2 ) );
    
      if ( warp->GetInitialAffineXform() ) 
	{
	classStream << (*warp->GetInitialAffineXform());
	} 
      else 
	{
	classStream << *this->m_InitialTransformation;
	}
      classStream << warp;
      classStream.End();
      }
    classStream.Close();
    }
}

void
ImagePairNonrigidRegistrationCommandLine::OutputIntermediate( const bool incrementCount )
{
  char path[PATH_MAX];
  if ( Studylist ) 
    {
    snprintf( path, sizeof( path ), "%s/level-%02d.list", Studylist, IntermediateResultIndex );
    } 
  else
    {
    snprintf( path, sizeof( path ), "level-%02d.list", IntermediateResultIndex );
    }
  this->OutputWarp( path );
  
  if ( incrementCount )
    ++IntermediateResultIndex;
}

CallbackResult 
ImagePairNonrigidRegistrationCommandLine::Register ()
{
  const double baselineTime = Timers::GetTimeProcess();
  CallbackResult Result = this->Superclass::Register();
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

void
cmtkImagePairNonrigidRegistrationCommandLineDispatchSIGUSR1( int sig )
{
  fprintf( stderr, "Received USR1 (%d) signal. Writing intermediate result #%d.\nNote that this result is not final.\n", sig, cmtk::ImagePairNonrigidRegistrationCommandLine::StaticThis->IntermediateResultIndex );

#ifndef _MSC_VER
  // set signal handler again.
  signal( sig, cmtkImagePairNonrigidRegistrationCommandLineDispatchSIGUSR1 );
#endif
  
  // write intermediate result. give "false" flag for index increment to 
  // preserve to final numbering of levels.
  cmtk::ImagePairNonrigidRegistrationCommandLine::StaticThis->OutputIntermediate( false );
}

