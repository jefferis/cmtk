/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkTimers.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkFilterVolume.h>
#include <Base/cmtkWarpXform.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkClassStream.h>

#include <Registration/cmtkSplineWarpGroupwiseRegistrationRMIFunctional.h>
#include <Registration/cmtkBestDirectionOptimizer.h>
#include <Registration/cmtkGroupwiseRegistrationOutput.h>

#include <vector>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

bool Verbose = false;
bool Help = false;

int DownsampleFrom = 4;
int DownsampleTo = 1;

int RefineTransformationGrid = 0;

cmtk::Types::Coordinate GridSpacing = 40.0;
bool GridSpacingExact = true;

bool ForceZeroSum = false;
bool ForceZeroSumNoAffine = false;
size_t ForceZeroSumFirstN = 0;
size_t NormalGroupFirstN = 0;

float SamplingDensity = -1.0;

bool DeactivateUninformative = true;
cmtk::Types::Coordinate PartialGradientThreshold = 0.01;

cmtk::Types::Coordinate SmoothSigmaFactor = 0.0;

cmtk::Types::Coordinate Accuracy = 0.01;
cmtk::Types::Coordinate Exploration = 0.25;
cmtk::Types::Coordinate OptimizerStepFactor = 0.5;
cmtk::Optimizer::ReturnType OptimizerDeltaFThreshold = 0;
bool OptimizerAggressive = true;
int OptimizerRepeatLevel = 2;

bool DisableOptimization = false;

const char* AffineGroupRegistration = NULL;

const char* OutputRootDirectory = NULL;
const char* OutputArchive = "groupwise_rmi_warp.xforms";
const char* OutputStudyListGroup = "groupwise_rmi_warp.list";
const char* OutputStudyListIndividual = "groupwise_rmi_warp_pairs";
const char* AverageImagePath = "average_groupwise_rmi_warp.nii";
cmtk::Interpolators::InterpolationEnum AverageImageInterpolation = cmtk::Interpolators::LINEAR;

float UserBackgroundValue = 0;
bool UserBackgroundFlag = false;

int
doMain( int argc, char ** argv )
{
#ifdef CMTK_BUILD_MPI
#  ifdef CMTK_BUILD_SMP
  const int threadLevelSupportedMPI =
    MPI::Init_thread( argc, argv, MPI::THREAD_FUNNELED );
  if ( threadLevelSupportedMPI < MPI::THREAD_FUNNELED )
    {
    cmtk::StdErr << "WARNING: your MPI implementation does not seem to support THREAD_FUNNELED.\n";
    }
#  else
  MPI::Init( argc, argv );
#  endif
#endif

  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Nonrigid RMI population registration" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] affineGroupwiseAlignment" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Image Registration" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose operation." );

    cl.AddOption( Key( 'd', "downsample-from" ), &DownsampleFrom, "Initial downsampling factor [4]." );
    cl.AddOption( Key( 'D', "downsample-to" ), &DownsampleTo, "Final downsampling factor [1]." );
    cl.AddOption( Key( 'r', "refine-grid" ), &RefineTransformationGrid, "Number of times to refine transformation grid [default: 0]." );

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'O', "output-root" ), &OutputRootDirectory, "Root directory for all output files." );
    cl.AddOption( Key( 'o', "output" ), &OutputArchive, "Output filename for groupwise registration archive." );
    cl.AddOption( Key( "output-average" ), &AverageImagePath, "Output filename for registered average image." );
    cl.AddSwitch( Key( "no-output-average" ), &AverageImagePath, (const char*)NULL, "Do not write average image." );
    cl.AddSwitch( Key( "average-cubic" ), &AverageImageInterpolation, cmtk::Interpolators::CUBIC, "Use cubic interpolation for average image (default: linear)" );
    cl.EndGroup();
    
    cl.AddOption( Key( 's', "sampling-density" ), &SamplingDensity, "Probabilistic sampling density [default: off]." );
    cl.AddOption( Key( 'B', "force-background" ), &UserBackgroundValue, "Force background pixels (outside FOV) to given value.", &UserBackgroundFlag );

    cl.AddOption( Key( 'p', "partial-gradient-thresh" ), &PartialGradientThreshold, "Threshold factor for partial gradient zeroing [default: 0.01; <0 turn off]" );
    cl.AddSwitch( Key( "deactivate-uninformative" ), &DeactivateUninformative, true, "Deactivate uninformative control points [default: on]" );
    cl.AddSwitch( Key( "activate-uninformative" ), &DeactivateUninformative, false, "Activate potentially uninformative control points [default: off]" );

    cl.AddOption( Key( "smooth" ), &SmoothSigmaFactor, "Sigma of Gaussian smoothing kernel in multiples of template image pixel size [default: off] )" );

    cl.AddOption( Key( "grid-spacing" ), &GridSpacing, "Control point grid spacing." );
    cl.AddSwitch( Key( "grid-spacing-fit" ), &GridSpacingExact, false, "Use grid spacing that fits volume FOV" );
    cl.AddSwitch( Key( 'z', "zero-sum" ), &ForceZeroSum, true, "Enforce zero-sum computation." );
    cl.AddSwitch( Key( "zero-sum-no-affine" ), &ForceZeroSumNoAffine, true, "Enforce zero-sum computation EXCLUDING affine components." );
    cl.AddOption( Key( 'N', "normal-group-first-n" ), &NormalGroupFirstN, "First N images are from the normal group and should be registered unbiased." );
    cl.AddOption( Key( 'Z', "zero-sum-first-n" ), &ForceZeroSumFirstN, "Enforce zero-sum computation for first N images.", &ForceZeroSum );

    cl.AddOption( Key( 'e', "exploration" ), &Exploration, "Exploration of optimization in pixels [0.25]" );
    cl.AddOption( Key( 'a', "accuracy" ), &Accuracy, "Accuracy of optimization in pixels [0.01]" );
    cl.AddOption( Key( 'S', "step-factor" ), &OptimizerStepFactor, "Step factor for successive optimization passes [0.5]" );
    cl.AddOption( Key( "delta-f-threshold" ), &OptimizerDeltaFThreshold, "Optional threshold to terminate optimization (level) if relative change of target function drops below this value." );

    cl.AddSwitch( Key( "disable-optimization" ), &DisableOptimization, true, "Disable optimization and output initial configuration." );
      
    cl.Parse( argc, const_cast<const char**>( argv ) );

    AffineGroupRegistration = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

#ifdef CMTK_BUILD_MPI
  const int mpiRank = MPI::COMM_WORLD.Get_rank();
  const int mpiSize = MPI::COMM_WORLD.Get_size();
  
  if ( Verbose )
    {
    std::cerr << "MPI started process with rank " << mpiRank << " out of " << mpiSize << "\n";
    }
#endif

  typedef cmtk::SplineWarpGroupwiseRegistrationRMIFunctional FunctionalType;
  FunctionalType::SmartPtr functional( new FunctionalType );
  functional->SetForceZeroSum( ForceZeroSum );
  functional->SetForceZeroSumFirstN( ForceZeroSumFirstN );
  functional->SetForceZeroSumNoAffine( ForceZeroSumNoAffine );
  if ( UserBackgroundFlag )
    functional->SetUserBackgroundValue( UserBackgroundValue );
  
  // disable parameters with less than 1% of maximum contribution
  functional->SetPartialGradientMode( (PartialGradientThreshold > 0) , PartialGradientThreshold );
  functional->SetDeactivateUninformativeMode( DeactivateUninformative );

  cmtk::ClassStream stream( AffineGroupRegistration, cmtk::ClassStream::READ );
  stream >> *functional;
  stream.Close();

  // this vector holds the original (not downsampled) images.
  std::vector<cmtk::UniformVolume::SmartPtr> imageListOriginal;
  functional->GetOriginalTargetImages( imageListOriginal );

  cmtk::UniformVolume::SmartPtr originalTemplateGrid = functional->GetTemplateGrid();
  functional->InitializeXforms( GridSpacing, GridSpacingExact ); // must do this before downsampling template grid

  const double timeBaselineProcess = cmtk::Timers::GetTimeProcess();

  if ( ! DisableOptimization )
    {
    cmtk::CoordinateVector v;
    for ( int downsample = DownsampleFrom; downsample >= DownsampleTo; 
	  downsample /= 2 )
      {
      if ( (RefineTransformationGrid > 0) && (downsample != DownsampleFrom) )
	{
	functional->RefineTransformationGrids();
	--RefineTransformationGrid;
	}
      functional->GetParamVector( v );
      
      functional->SetTemplateGrid( originalTemplateGrid, downsample ); 
      cmtk::UniformVolume::SmartPtr templateGrid = functional->GetTemplateGrid();

      if ( SmoothSigmaFactor > 0.0 )
	{
	functional->SetGaussianSmoothImagesSigma( SmoothSigmaFactor * templateGrid->GetMinDelta() );
	functional->SetTargetImages( imageListOriginal );
	}

      if ( Verbose )
	{
	cmtk::StdErr.printf( "Template grid is %d x %d x %d pixels of size %f x %f x %f\n",
			     templateGrid->m_Dims[0], templateGrid->m_Dims[1], templateGrid->m_Dims[2], templateGrid->m_Delta[0], templateGrid->m_Delta[1], templateGrid->m_Delta[2] );
	}
      
      if ( SamplingDensity > 0 )
	{
	functional->SetProbabilisticSampleDensity( SamplingDensity );
	functional->SetProbabilisticSampleUpdatesAfter( 10 );
	}
      
      functional->AllocateStorage();

      cmtk::BestDirectionOptimizer optimizer( OptimizerStepFactor );
      optimizer.SetAggressiveMode( OptimizerAggressive );
      optimizer.SetDeltaFThreshold( OptimizerDeltaFThreshold );
      optimizer.SetRepeatLevelCount( OptimizerRepeatLevel );
      optimizer.SetFunctional( functional );

      try
	{
	// do we have a normal subgroup?
	if ( NormalGroupFirstN )
	  {
	  // yes: first run normal group by itself
	  cmtk::StdErr << "Running normal subgroup...\n";
	  functional->SetForceZeroSum( ForceZeroSum );
	  functional->SetActiveImagesFromTo( 0, NormalGroupFirstN );
	  functional->SetActiveXformsFromTo( 0, NormalGroupFirstN );
	  optimizer.Optimize( v, Exploration * templateGrid->GetMinDelta(), Accuracy * templateGrid->GetMinDelta() );
	  
	  // second: run abnormal group, but keep using normal group's data for reference
	  cmtk::StdErr << "Running diseased subgroup...\n";
	  functional->SetForceZeroSum( false ); // no point here
	  functional->SetActiveImagesFromTo( 0, imageListOriginal.size() );
	  functional->SetActiveXformsFromTo( NormalGroupFirstN, imageListOriginal.size() );
	  optimizer.Optimize( v, Exploration * templateGrid->GetMinDelta(), Accuracy * templateGrid->GetMinDelta() );
	  }
	else
	  {
	  optimizer.Optimize( v, Exploration * templateGrid->GetMinDelta(), Accuracy * templateGrid->GetMinDelta() );
	  }
	}
      catch ( cmtk::GroupwiseRegistrationFunctionalBase::BadXform )
	{
	cmtk::StdErr << "FAILED: at least one image has too few pixels in the template area.\n";
	return 1;
	}	      
      }
    }

  // determine and print CPU time (by node, if using MPI)
  const double timeElapsedProcess = cmtk::Timers::GetTimeProcess() - timeBaselineProcess;
#ifdef CMTK_BUILD_MPI    
  std::vector<double> timeElapsedByNodeProcess( mpiSize );
  MPI::COMM_WORLD.Gather( &timeElapsedProcess, 1, MPI::DOUBLE, &timeElapsedByNodeProcess[0], 1, MPI::DOUBLE, 0 /*root*/ );
  
  if ( mpiRank == 0 )
    {
    cmtk::StdErr << "Process CPU time [s] by node:\n";
    for ( int node = 0; node < mpiSize; ++node )
      {
      cmtk::StdErr.printf( "%d\t%lf\n", node, timeElapsedByNodeProcess[node] );
      }
    }
#else
  cmtk::StdErr.printf( "Process CPU time [s]: %lf\n", timeElapsedProcess );  
#endif  

  cmtk::GroupwiseRegistrationOutput output;
  output.SetFunctional( functional );
  output.SetOutputRootDirectory( OutputRootDirectory );

  output.WriteGroupwiseArchive( OutputArchive );
  output.WriteXformsSeparateArchives( OutputStudyListIndividual, AverageImagePath );
  output.WriteAverageImage( AverageImagePath, AverageImageInterpolation );

#ifdef CMTK_BUILD_MPI    
  MPI::Finalize();
#endif

  return 0;
}

#include "cmtkSafeMainMPI"
