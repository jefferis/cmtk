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

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>
#include <cmtkTimers.h>

#include <cmtkUniformVolume.h>
#include <cmtkFilterVolume.h>

#include <cmtkVolumeIO.h>
#include <cmtkUniformVolume.h>
#include <cmtkClassStream.h>

#include <cmtkAffineCongealingFunctional.h>
#include <cmtkBestDirectionOptimizer.h>

#include <cmtkGroupwiseRegistrationOutput.h>

#include <vector>

bool Verbose = false;

int DownsampleFrom = 4;
int DownsampleTo = 1;

byte UserBackgroundValue = 0;
bool UserBackgroundFlag = false;

float SamplingDensity = -1.0;
bool UseSamplingDensity = false;

bool ForceZeroSum = false;
unsigned int ForceZeroSumFirstN = 0;
unsigned int  NormalGroupFirstN = 0;

bool UseSmoothSigmaFactor = false;
cmtk::Types::Coordinate SmoothSigmaFactor = -1.0;

bool UseNumberOfHistogramBins = false;
unsigned int NumberOfHistogramBins = 0;
bool CropImageHistograms = false;

const char* PreDefinedTemplatePath = NULL;
cmtk::UniformVolume::SmartPtr PreDefinedTemplate;
bool UseTemplateData = false;
bool TransformationsFromArchive = false;

const char* OutputRootDirectory = NULL;
const char* OutputArchive = "congeal.xforms";
const char* OutputStudyListGroup = "congeal.list";
const char* OutputStudyListIndividual = "congeal_pairs";
const char* AverageImagePath = "average_congeal.hdr";
cmtk::Interpolators::InterpolationEnum AverageImageInterpolation = cmtk::Interpolators::LINEAR;

std::vector<int> NumberDOFs;

bool AlignCentersOfMass = false;
bool InitScales = false;
bool HistogramMatching = false;

cmtk::Types::Coordinate Accuracy = 0.01;
cmtk::Types::Coordinate Exploration = 0.25;
cmtk::Types::Coordinate OptimizerStepFactor = 0.5;
cmtk::Optimizer::ReturnType OptimizerDeltaFThreshold = 0;
bool DisableOptimization = false;
bool OptimizerAggressive = false;
int OptimizerRepeatLevel = 5;

// this vector holds all target image filenames
std::vector<const char*> fileNameList;

// this vector holds the original (not downsampled) images.
std::vector<cmtk::UniformVolume::SmartPtr> imageListOriginal;

int
main( int argc, char* argv[] )
{
#ifdef CMTK_BUILD_MPI
#  ifdef CMTK_BUILD_SMP
  const int threadLevelSupportedMPI = MPI::Init_thread( argc, argv, MPI::THREAD_FUNNELED );
  if ( threadLevelSupportedMPI < MPI::THREAD_FUNNELED )
    {
    cmtk::StdErr << "WARNING: your MPI implementation does not seem to support THREAD_FUNNELED.\n";
    }
#  else
  MPI::Init( argc, argv );
#  endif
#endif

#ifdef CMTK_BUILD_MPI
  const int mpiRank = MPI::COMM_WORLD.Get_rank();
  const int mpiSize = MPI::COMM_WORLD.Get_size();
#endif

  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Affine population registration" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool registers a population of input images simultaneously, without a template, using the 'congealing' algorithm" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] image0 [image1 ...]" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Image Registration" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose operation." )->SetProperties( cmtk::CommandLine::PROPS_NOXML );

    cl.BeginGroup( "Template", "Template Image Options" );
    cl.AddOption( Key( 't', "template" ), &PreDefinedTemplatePath, "Input filename for pre-defined template image." );
    cl.AddOption( Key( 'T', "template-with-data" ), &PreDefinedTemplatePath, "Use user-supplied template images's pixel data in registration", &UseTemplateData );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'O', "output-root" ), &OutputRootDirectory, "Root directory for all output files." );
    cl.AddOption( Key( 'o', "output" ), &OutputArchive, "Output filename for groupwise registration archive." );
    cl.AddOption( Key( "output-average" ), &AverageImagePath, "Output filename for registered average image." );
    cl.AddSwitch( Key( "average-linear" ), &AverageImageInterpolation, cmtk::Interpolators::LINEAR, "Use linear interpolation for average image" );
    cl.AddSwitch( Key( "average-cubic" ), &AverageImageInterpolation, cmtk::Interpolators::CUBIC, "Use cubic interpolation for average image" );
    cl.AddSwitch( Key( "no-output-average" ), &AverageImagePath, (const char*)NULL, "Do not write average image." );
    cl.EndGroup();

    cl.BeginGroup( "Multuresolution", "Multuresolution Parameters" );
    cl.AddOption( Key( 'd', "downsample-from" ), &DownsampleFrom, "Initial downsampling factor" );
    cl.AddOption( Key( 'D', "downsample-to" ), &DownsampleTo, "Final downsampling factor." );
    cl.AddOption( Key( 's', "sampling-density" ), &SamplingDensity, "Probabilistic sampling density. Legal values between 0 and 1.", &UseSamplingDensity );
    cl.EndGroup();

    cl.BeginGroup( "Image", "Image Options and Operations" );
    cl.AddOption( Key( 'B', "force-background" ), &UserBackgroundValue, "Force background pixels (outside FOV) to given (bin) value.", &UserBackgroundFlag );
    cl.AddOption( Key( 'H', "histogram-bins" ), &NumberOfHistogramBins, "Manually set number of histogram bins for entropy evaluation", &UseNumberOfHistogramBins );
    cl.AddSwitch( Key( "crop-histograms" ), &CropImageHistograms, true, "Crop image histograms to make better use of histogram bins." );
    cl.AddOption( Key( "smooth" ), &SmoothSigmaFactor, "Sigma of Gaussian smoothing kernel in multiples of template image pixel size.", &UseSmoothSigmaFactor );
    cl.AddSwitch( Key( "match-histograms" ), &HistogramMatching, true, "Match all image histograms to template data (or first image, if no template image is given)" );
    cl.EndGroup();

    cl.BeginGroup( "Transformation", "Transformation Parameters" );
    cl.AddVector( Key( "dofs" ), NumberDOFs, "Add DOFs to list [default: one pass, 6 DOF]" );
    cl.AddSwitch( Key( 'z', "zero-sum" ), &ForceZeroSum, true, "Enforce zero-sum computation." );
    cl.AddOption( Key( 'N', "normal-group-first-n" ), &NormalGroupFirstN, "First N images are from the normal group and should be registered unbiased." );
    cl.AddOption( Key( 'Z', "zero-sum-first-n" ), &ForceZeroSumFirstN, "Enforce zero-sum computation for first N images.", &ForceZeroSum );
    
    cl.BeginGroup( "Initialization", "Transformation Initialization" );
    cl.AddSwitch( Key( "align-bounding-boxes" ), &AlignCentersOfMass, false, "Initially align centers of bounding boxes of all images by translations" );
    cl.AddSwitch( Key( "align-centers-of-mass" ), &AlignCentersOfMass, true, "Initially align centers of mass by translations" );
    cl.AddSwitch( Key( "init-scales" ), &InitScales, true, "Initialize scale factors using first-order moments" );
    cl.EndGroup();

    cl.BeginGroup( "Optimization", "Optimization Parameters" );
    cl.AddOption( Key( 'e', "exploration" ), &Exploration, "Exploration of optimization in pixels" );
    cl.AddOption( Key( 'a', "accuracy" ), &Accuracy, "Accuracy of optimization in pixels" );
    cl.AddOption( Key( 'r', "repeat-level" ), &OptimizerRepeatLevel, "Number of repetitions per optimization level" );
    cl.AddOption( Key( 'S', "step-factor" ), &OptimizerStepFactor, "Step factor for successive optimization passes" );
    cl.AddOption( Key( "delta-f-threshold" ), &OptimizerDeltaFThreshold, "Optional threshold to terminate optimization (level) if relative change of target function drops below this value." );
    cl.AddSwitch( Key( "disable-optimization" ), &DisableOptimization, true, "Disable optimization and output initial configuration." );
    cl.EndGroup();
      
    cl.Parse();

    const char* next = cl.GetNext();
    while ( next )
      {      
      fileNameList.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  if ( !NumberDOFs.size() )
    NumberDOFs.push_back( 6 );

  cmtk::AffineCongealingFunctional::SmartPtr functional( new cmtk::AffineCongealingFunctional );
  functional->SetForceZeroSum( ForceZeroSum );
  functional->SetForceZeroSumFirstN( ForceZeroSumFirstN );
  functional->SetFreeAndRereadImages( ! (AlignCentersOfMass || InitScales) );
  functional->SetCropImageHistograms( CropImageHistograms );

  if ( UserBackgroundFlag )
    functional->SetUserBackgroundValue( UserBackgroundValue );

  if ( UseNumberOfHistogramBins )
    {
    // must be done IMMEDIATELY after creating the functional!
    // otherwise, scaling and conversion of input images goes
    // WRONG!
    functional->SetNumberOfHistogramBins( NumberOfHistogramBins );
    }

  if ( cmtk::FileFormat::Identify( fileNameList[0] ) == cmtk::FILEFORMAT_TYPEDSTREAM )
    {
    if ( fileNameList.size() > 1 )
      {
      cmtk::StdErr << "First input file is an archive, but additional arguments are given.\n"
		<< "I am terminating just to make sure not to do something stupid.\n";
      exit( 1 );
      }
    
    cmtk::ClassStream inStream( fileNameList[0], cmtk::ClassStream::READ );
    if ( inStream.IsValid() )
      {
      inStream >> *functional;
      PreDefinedTemplate = functional->GetTemplateGrid();
      imageListOriginal = functional->GetOriginalTargetImages();
      TransformationsFromArchive = true;
      }
    else
      {
      cmtk::StdErr << "Could not open input groupwise archive " << fileNameList[0] << "\n";
      exit( 1 );
      }
    }
  else
    {
    int idx = 0;
    for ( std::vector<const char*>::const_iterator fnIt = fileNameList.begin(); fnIt != fileNameList.end(); ++fnIt, ++idx )
      {
      cmtk::UniformVolume::SmartPtr nextImage;
#ifdef CMTK_BUILD_MPI
      if ( (idx % mpiSize) == mpiRank )
#endif
	{
	cmtk::UniformVolume::SmartPtr image( cmtk::VolumeIO::ReadOriented( *fnIt, Verbose ) );
	if ( ! image || ! image->GetData() )
	  {
	  cmtk::StdErr << "ERROR: Could not read image " << *fnIt << "\n";
	  exit( 1 );
	  }
	nextImage = image;
	}
      
#ifdef CMTK_BUILD_MPI
      else
	{
	cmtk::UniformVolume::SmartPtr image( cmtk::VolumeIO::ReadGridOriented( *fnIt, false ) );
	if ( ! image )
	  {
	  cmtk::StdErr << "ERROR: Could not read image " << *fnIt << "\n";
	  exit( 1 );
	  }
	nextImage = image;
	}
#endif
      imageListOriginal.push_back( nextImage );
      }
    
    if ( PreDefinedTemplatePath )
      {
      if ( UseTemplateData )
	{
	PreDefinedTemplate = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( PreDefinedTemplatePath, Verbose ) );
	}
      else
	{
	PreDefinedTemplate = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( PreDefinedTemplatePath, Verbose ) );
	}

      if ( ! PreDefinedTemplate )
	{
	cmtk::StdErr << "ERROR: could not read template grid/image " << PreDefinedTemplatePath << "\n";
	exit( 2 );
	}
      }
    }

  if ( HistogramMatching )
    {
    const cmtk::TypedArray* referenceDataForHistogramMatching = NULL; 
    if ( PreDefinedTemplate && UseTemplateData )
      {
      referenceDataForHistogramMatching = PreDefinedTemplate->GetData();
      }
    if ( !referenceDataForHistogramMatching )
      {
      referenceDataForHistogramMatching = imageListOriginal[0]->GetData();
      }
    
    for ( size_t idx = 1; idx < imageListOriginal.size(); ++idx )
      {
      imageListOriginal[idx]->GetData()->MatchHistogramToReference( referenceDataForHistogramMatching );
      }
    }

  if ( UseSamplingDensity )
    {
    functional->SetProbabilisticSampleDensity( SamplingDensity );
    }
  
  const double timeBaselineProcess = cmtk::Timers::GetTimeProcess();

  const int downsampleFrom = std::max( DownsampleFrom, DownsampleTo );
  const int downsampleTo = std::min( DownsampleFrom, DownsampleTo );

  cmtk::CoordinateVector v;
  for ( int downsample = downsampleFrom; downsample >= downsampleTo; downsample = downsample?downsample/2:-1 )
    {
    functional->SetTargetImages( imageListOriginal );
    if ( PreDefinedTemplate )
      functional->SetTemplateGrid( PreDefinedTemplate, std::max( 1, downsample ), UseTemplateData );
    else
      functional->CreateTemplateGridFromTargets( imageListOriginal, std::max( 1, downsample ) );
    
    cmtk::UniformVolume::SmartPtr templateGrid = functional->GetTemplateGrid();
    if ( UseSmoothSigmaFactor && downsample )
      {
      functional->SetGaussianSmoothImagesSigma( SmoothSigmaFactor * templateGrid->GetMinDelta() );
      }
    functional->AllocateStorage();    
    
    if ( Verbose )
      {
      cmtk::StdErr.printf( "Template grid is %d x %d x %d pixels of size %f x %f x %f\n",
			   templateGrid->m_Dims[0], templateGrid->m_Dims[1], templateGrid->m_Dims[2], 
			   templateGrid->m_Delta[0], templateGrid->m_Delta[1], templateGrid->m_Delta[2] );
      }
    
    if ( downsampleFrom == downsample )
      {
      if ( !TransformationsFromArchive )
	functional->InitializeXforms( true /*alignCenters*/, AlignCentersOfMass, InitScales );
      functional->SetFreeAndRereadImages( true ); // can now get rid of unused original images
      functional->GetParamVector( v );
      }
    else
      {
      functional->SetParamVector( v );
      }
    
    if ( ! DisableOptimization )
      {
      cmtk::BestDirectionOptimizer optimizer( OptimizerStepFactor );
      optimizer.SetAggressiveMode( OptimizerAggressive );
      optimizer.SetRepeatLevelCount( OptimizerRepeatLevel );
      optimizer.SetDeltaFThreshold( OptimizerDeltaFThreshold );
      optimizer.SetFunctional( functional );
      
      for ( std::vector<int>::const_iterator itDOF = NumberDOFs.begin(); itDOF != NumberDOFs.end(); ++itDOF )
	{
	functional->SetXformNumberDOFs( *itDOF );
	
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
      }
    }
  
  // determine and print CPU time (by node, if using MPI)
  const double timeElapsedProcess = cmtk::Timers::GetTimeProcess() - timeBaselineProcess;
#ifdef CMTK_BUILD_MPI    
  std::vector<float> timeElapsedByNodeProcess( mpiSize );
  MPI::COMM_WORLD.Gather( &timeElapsedProcess, 1, MPI::FLOAT, &timeElapsedByNodeProcess[0], 1, MPI::FLOAT, 0 /*root*/ );

  if ( mpiRank == 0 )
    {
    cmtk::StdErr << "Process CPU time [s] by node:\n";
    for ( int node = 0; node < mpiSize; ++node )
      {
      cmtk::StdErr.printf( "%d\t%f\n", node, timeElapsedByNodeProcess[node] );
      }
    }
#else
  cmtk::StdErr.printf( "Process CPU time [s]: %f\n", timeElapsedProcess );  
#endif  

  functional->SetTargetImages( imageListOriginal );
  if ( PreDefinedTemplate )
    functional->SetTemplateGrid( PreDefinedTemplate );
  else
    functional->CreateTemplateGridFromTargets( imageListOriginal );

  cmtk::GroupwiseRegistrationOutput output;
  output.SetFunctional( functional );
  output.SetOutputRootDirectory( OutputRootDirectory );

  if ( ! UseTemplateData )
    {
    PreDefinedTemplatePath = AverageImagePath;
    output.SetExistingTemplatePath( false );
    }
  else
    {
    output.SetExistingTemplatePath( true );
    }

  output.WriteGroupwiseArchive( OutputArchive );
  output.WriteXformsSeparateArchives( OutputStudyListIndividual, PreDefinedTemplatePath );
  output.WriteAverageImage( AverageImagePath, AverageImageInterpolation, UseTemplateData );

#ifdef CMTK_BUILD_MPI    
  MPI::Finalize();
#endif

  return 0;
}
