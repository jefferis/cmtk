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

#include <cmtkElasticRegistrationCommandLine.h>

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

#include <cmtkVoxelMatchingElasticFunctional.h>
#include <cmtkSymmetricElasticFunctional.h>
#include <cmtkProtocolCallback.h>

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

ElasticRegistrationCommandLine* ElasticRegistrationCommandLine::StaticThis = NULL;

ElasticRegistrationCommandLine
::ElasticRegistrationCommandLine
( int argc, char *argv[] ) 
{
  Metric = 0;
  Algorithm = 3;

  CoarsestResolution = -1;
  Exploration = 4.0;
  GridSpacing = 15;
  ExactGridSpacing = 0;
  Accuracy = 0.1;
  Sampling = 1.0;
  Studylist = Protocol = Time = NULL;

  this->m_OutputIntermediate = 0;
  Verbose = 0;

  InputStudylist = NULL;
  const char *InitialStudylist = NULL;
  IntermediateResultIndex = 0;

  bool Switch = false;

  this->RigidityConstraintMapFilename = NULL;

  bool forceOutsideFlag = false;
  Types::DataItem forceOutsideValue = 0;

  const char* clArg1; // input studylist or reference image
  const char* clArg2; // empty or floating image

  try
    {
    CommandLine cl( argc, argv );    
    cl.SetProgramInfo( CommandLine::PRG_TITLE, "B-spline nonrigid registration" );
    cl.SetProgramInfo( CommandLine::PRG_DESCR, "This program performs nonrigid image registration using multi-resolution optimization of voxel-based image similarity measures "
		       "and a multi-resolution B-spline transformation model." );
    cl.SetProgramInfo( CommandLine::PRG_SYNTX, "[options] [refImage fltImage | initialStudylist]" );

    typedef CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &this->Verbose, true, "Verbose peration" );
    cl.AddSwitch( Key( 'q', "quiet" ), &Verbose, false, "Quiet mode" );

    cl.BeginGroup( "Transformation", "Transformation parameters" );
    cl.AddOption( Key( 'g', "grid-spacing" ), &this->GridSpacing, "Control point grid spacing" );
    cl.AddSwitch( Key( "exact-spacing" ), &this->ExactGridSpacing, true, "Use exact control point spacing; do not modify spacing to fit reference image bounding box" );
    cl.AddOption( Key( "refine" ), &this->RefineGrid, "Number of refinements (control point grid resolution levels)" );
    cl.AddSwitch( Key( "no-delay-refine" ), &this->DelayRefineGrid, false, "Always refine control point grid (up to maximum number of levels) when switching next higher image resolution [default]" );
    cl.AddSwitch( Key( "delay-refine" ), &this->DelayRefineGrid, true, "Delay control point grid refinement; first switch to next higher image resolution" );

    cl.AddOption( Key( "ignore-edge" ), &this->IgnoreEdge, "Ignore n control point layers along each image face" );
    cl.AddOption( Key( "restrict" ), &this->RestrictToAxes, "Restrict deformation to coordinate dimension(s) [one or more of 'x','y','z']" );

    cl.AddSwitch( Key( "adaptive-fix" ), &this->AdaptiveFixParameters, true, "Adaptive fixing of control points in image background [default]" );
    cl.AddSwitch( Key( "no-adaptive-fix" ), &this->AdaptiveFixParameters, false, "Disable adaptive fixing of control points; optimize all deformation parameters" );
    cl.AddOption( Key( "adaptive-fix-thresh" ), &this->AdaptiveFixThreshFactor, "Threshold factor for entropy criterion to fix local control points" );
    cl.AddSwitch( Key( "fast" ), &this->FastMode, true, "Fast computation mode: take some numerical short cuts to save about 80% computation time [default]" );
    cl.AddSwitch( Key( "accurate" ), &this->FastMode, false, "Accurate computation mode: may give slightly better results after substantially longer computation" );

    cl.AddSwitch( Key( 'S', "switch" ), &Switch, true, "Switch reference and floating image" );
    cl.AddSwitch( Key( 'x', "exchange" ), &this->ForceSwitchVolumes, true, "Exchange reference and floating image");
    cl.AddOption( Key( "initial" ), &InitialStudylist, "Initialize transformation from given path" );
    cl.EndGroup();

    cl.BeginGroup( "Optimization", "Optimization parameters" );
    cl.AddOption( Key( 'e', "exploration" ), &this->Exploration, "Search space exploration (initial step size)" );
    cl.AddOption( Key( 'a', "accuracy" ), &this->Accuracy, "Search accuracy (initial step size)" );
    cl.AddOption( Key( 'f', "stepfactor" ), &this->OptimizerStepFactor, "Factor for search step size reduction. Must be > 0.0 and < 1.0 [default: 0.5]" );

    cl.AddSwitch( Key( "maxnorm" ), &this->UseMaxNorm, true, "Normalized optimization gradient using maximum norm [default]" );
    cl.AddSwitch( Key( "no-maxnorm" ), &this->UseMaxNorm, false, "Disable optimization gradient maximum normalication; use Euclid norm instead" );

    cl.AddOption( Key( "jacobian-weight" ), &this->JacobianConstraintWeight, "Weight for Jacobian-based local volume preservation constraint" );
    cl.AddOption( Key( "energy-weight" ), &this->GridEnergyWeight, "Weight for grid bending energy constraint" );
    cl.AddOption( Key( "rigidity-weight" ), &this->RigidityConstraintWeight, "Weight for local rigidity constraint" );
    cl.AddOption( Key( "landmark-weight" ), &this->LandmarkErrorWeight, "Weight for landmark misregistration registration" );
    cl.AddOption( Key( "ic-weight" ), &this->InverseConsistencyWeight, "Weight for inverse consistency constraint" );
    cl.AddOption( Key( "relax" ), &this->RelaxWeight, "Weight relaxation factor for alternating under-constrained iterations" );

    cl.AddOption( Key( "rigidity-weight-map" ), &this->RigidityConstraintMapFilename, "Filename for rigidity map weight image" );
    cl.EndGroup();

    cl.BeginGroup( "Resolution", "Image resolution parameters" );
    cl.AddOption( Key( 's', "sampling" ), &this->Sampling, "Image sampling (finest resampled image resolution)" );
    cl.AddOption( Key( "coarsest" ), &this->CoarsestResolution, "Upper limit for image sampling in multiresolution hierarchy" );

    cl.AddSwitch( Key( "use-original-data" ), &this->UseOriginalData, true, "Use original data in full resolution as final level [default]" );
    cl.AddSwitch( Key( "omit-original-data" ), &this->UseOriginalData, false, "Do NOT use original data in full resolution" );
    cl.EndGroup();

    cl.BeginGroup( "Images", "Image data" );
    cl.AddSwitch( Key( "nmi" ), &this->Metric, 0, "Normalized Mutual Information metric" );
    cl.AddSwitch( Key( "mi" ), &this->Metric, 1, "Standard Mutual Information metric" );
    cl.AddSwitch( Key( "cr" ), &this->Metric, 2, "Correlation Ratio metric" );
    cl.AddSwitch( Key( "msd" ), &this->Metric, 4, "Mean Squared Difference metric" );
    cl.AddSwitch( Key( "ncc" ), &this->Metric, 5, "Normalized Cross Correlation metric" );

    cl.AddOption( Key( "force-outside-value" ), &forceOutsideValue, "Force values outside field of view to this value rather than drop incomplete pixel pairs", &forceOutsideFlag );

    this->m_PreprocessorRef.AttachToCommandLine( cl );
    this->m_PreprocessorFlt.AttachToCommandLine( cl );

    cl.BeginGroup( "Output", "Output parameters" );
    cl.AddOption( Key( 'o', "outlist" ), &this->Studylist, "Output path for final transformation" );
    cl.AddOption( Key( 'p', "protocol" ), &this->Protocol, "Optimization protocol output file name" );
    cl.AddOption( Key( 't', "time" ), &this->Time, "Computation time statistics output file name" );
    cl.AddSwitch( Key( "output-intermediate" ), &this->m_OutputIntermediate, true, "Write transformation for each level [default: only write final transformation]" );
    cl.EndGroup();

    cl.Parse();

    clArg1 = cl.GetNext();
    clArg2 = cl.GetNextOptional();

//      cl.PrintHelp( "B-spline nonrigid registration", "[options] [studylist | refImage fltImage]");
    }
  catch ( CommandLine::Exception ex )
    {
    StdErr << ex << "\n";
    exit( 1 );
    }

  if ( (OptimizerStepFactor <= 0) || (OptimizerStepFactor >= 1) ) 
    {
    StdErr << "ERROR: step factor value " << OptimizerStepFactor << " is invalid. Must be in range (0..1)\n";
    exit( 1 );
    }

  if ( clArg2 ) 
    {
    this->SetInitialXform( AffineXform::SmartPtr( new AffineXform() ) );
    
    Study1 = const_cast<char*>( clArg1 );
    Study2 = const_cast<char*>( clArg2 );
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
      this->SetInitialXform( affineXform->GetInverse() );
      }
    else
      {
      Study2 = classStream.ReadString( "model_study" );
      AffineXform::SmartPtr affineXform;
      classStream >> affineXform;
      this->SetInitialXform( affineXform );
      }
    
    classStream.Close();
    }

  /// Was an initial studylist given? If so, get warp 
  if ( InitialStudylist ) 
    {
    Xform::SmartPtr initialXform( XformIO::Read( InitialStudylist ) );
    AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( initialXform );
    if ( affineXform )
      {
      this->SetInitialXform( affineXform->GetInverse() );
      }
    else
      {
      InitialWarpXform = SplineWarpXform::SmartPtr::DynamicCastFrom( initialXform );
      }
    }

  // Did user ask for exchanged reference/floating images?
  if ( Switch ) 
    {
    AffineXform::SmartPtr affineXform( dynamic_cast<AffineXform*>( InitialXform->MakeInverse() ) );
    this->SetInitialXform( affineXform );
    
    char *swap = Study1;
    Study1 = Study2;
    Study2 = swap;
    }
  
  UniformVolume::SmartPtr volume( VolumeIO::ReadOriented( Study1, Verbose ) );
  if ( !volume ) throw ConstructorFailed();
  this->SetVolume_1( UniformVolume::SmartPtr( this->m_PreprocessorRef.GetProcessedImage( volume ) ) );

  volume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( Study2, Verbose ) );
  if ( !volume ) throw ConstructorFailed();
  this->SetVolume_2( UniformVolume::SmartPtr( this->m_PreprocessorFlt.GetProcessedImage( volume ) ) );
  
  if ( this->RigidityConstraintMapFilename )
    {
    UniformVolume::SmartPtr rigidityWeightMap( VolumeIO::ReadOriented( this->RigidityConstraintMapFilename, Verbose ) );
    if ( rigidityWeightMap )
      {
      this->SetRigidityConstraintMap( rigidityWeightMap );
      }
    else
      {
      StdErr << "ERROR: rigidity constraint mapcould not be read from " << this->RigidityConstraintMapFilename << "\n";
      exit( 1 );
      }
    }

  if ( forceOutsideFlag )
    {
    this->SetForceOutside( true, forceOutsideValue );
    }

  if ( Protocol )
    {
    RegistrationCallback::SmartPtr callback( new ProtocolCallback( Protocol ) );
    this->SetCallback( callback );
    }
}

ElasticRegistrationCommandLine::~ElasticRegistrationCommandLine ()
{
#ifdef HAVE_SIGRELSE
  // release signal handler
  sigrelse( SIGUSR1 );
#endif
}

CallbackResult
ElasticRegistrationCommandLine
::InitRegistration ()
{
  CallbackResult result = this->Superclass::InitRegistration();
  if ( result != CALLBACK_OK )
    return result;

  if ( this->m_OutputIntermediate )
    this->OutputIntermediate();
  
  // register signal handler for intermediate result output.
  ElasticRegistrationCommandLine::StaticThis = this;
#ifndef _MSC_VER
  signal( SIGUSR1, ElasticRegistrationCommandLine::DispatchSIGUSR1 );
#endif

  return CALLBACK_OK;
}
	
void
ElasticRegistrationCommandLine
::OutputResult
( const CoordinateVector* )
  const
{
  if ( Studylist ) 
    this->OutputWarp( Studylist );
}

void
ElasticRegistrationCommandLine
::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, 
  const int index, const int total )
{
  if ( Verbose )
    fprintf( stderr, "\rEntering resolution level %d out of %d...\n", index, total );
  
  this->Superclass::EnterResolution( v, f, index, total );
}

int
ElasticRegistrationCommandLine::DoneResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, const int index, const int total )
{
  if ( this->m_OutputIntermediate )
    this->OutputIntermediate();
  return this->Superclass::DoneResolution( v, f, index, total );
}

void
ElasticRegistrationCommandLine::OutputWarp ( const char* path ) const
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
  classStream.WriteInt( "algorithm", Algorithm );
  classStream.WriteBool( "use_maxnorm", UseMaxNorm );
  classStream.WriteDouble( "exploration", Exploration );
  classStream.WriteDouble( "accuracy", Accuracy );
  classStream.WriteDouble( "min_sampling", Sampling );
  classStream.WriteDouble( "coarsest_resolution", CoarsestResolution );
  classStream.WriteBool( "use_original_data", UseOriginalData );
  classStream.WriteInt( "metric", Metric );
  classStream.WriteDouble( "optimizer_step_factor", OptimizerStepFactor );
  classStream.WriteDouble( "grid_spacing", GridSpacing );
  classStream.WriteInt( "ignore_edge", IgnoreEdge );
  classStream.WriteDouble( "jacobian_constraint_weight", JacobianConstraintWeight );
  classStream.WriteDouble( "rigidity_constraint_weight", RigidityConstraintWeight );
  if ( this->RigidityConstraintMapFilename )
    {
    classStream.WriteString( "rigidity_constraint_map_filename", RigidityConstraintMapFilename );
    }
  classStream.WriteDouble( "energy_constraint_weight", GridEnergyWeight );
  classStream.WriteDouble( "inverse_consistency_weight", InverseConsistencyWeight );
  classStream.WriteDouble( "weight_relaxation", RelaxWeight );
  classStream.WriteDouble( "landmark_error_weight", LandmarkErrorWeight );
  classStream.WriteBool( "force_switch", ForceSwitchVolumes );
  classStream.WriteInt( "refine_grid", RefineGrid );
  classStream.WriteBool( "delay_refine_grid", DelayRefineGrid );
  classStream.WriteBool( "adaptive_fix_parameters", AdaptiveFixParameters );
  classStream.WriteDouble( "adaptive_fix_parameters_thresh", AdaptiveFixThreshFactor );

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
      if ( SwitchVolumes ) 
	{
	classStream.WriteString( "reference_study", CompressedStream::GetBaseName( Study2 ) );
	classStream.WriteString( "floating_study", CompressedStream::GetBaseName( Study1 ) );
	} 
      else
	{
	classStream.WriteString( "reference_study", CompressedStream::GetBaseName( Study1 ) );
	classStream.WriteString( "floating_study", CompressedStream::GetBaseName( Study2 ) );
	}
      
      if ( warp->GetInitialAffineXform() ) 
	{
	classStream << (*warp->GetInitialAffineXform());
	} 
      else 
	{
	classStream << (*InitialXform->GetInverse());
	}
      classStream << warp;
      classStream.End();
      }
    classStream.Close();
    }
}

void
ElasticRegistrationCommandLine::OutputIntermediate( const bool incrementCount )
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

CallbackResult ElasticRegistrationCommandLine::Register ()
{
  const double baselineTime = cmtk::Timers::GetTimeProcess();
  CallbackResult Result = this->Superclass::Register();
  const int elapsed = static_cast<int>( cmtk::Timers::GetTimeProcess() - baselineTime );

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

void
ElasticRegistrationCommandLine::DispatchSIGUSR1( int sig )
{
  fprintf( stderr, "Received USR1 (%d) signal. Writing intermediate result #%d.\n"
	   "Note that this result is not final.\n", sig, StaticThis->IntermediateResultIndex );

#ifndef _MSC_VER
  // set signal handler again.
  signal( sig, ElasticRegistrationCommandLine::DispatchSIGUSR1 );
#endif

  // write intermediate result. give "false" flag for index increment to 
  // preserve to final numbering of levels.
  StaticThis->OutputIntermediate( false );
}

} // namespace cmtk
