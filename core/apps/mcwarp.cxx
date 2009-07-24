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
#include <cmtkVolumeIO.h>

#include <cmtkUniformVolume.h>
#include <cmtkUniformVolumeInterpolator.h>
#include <cmtkLinearInterpolator.h>
#include <cmtkCubicInterpolator.h>

#include <cmtkMultiChannelRMIRegistrationFunctional.h>
#include <cmtkMultiChannelHistogramRegistrationFunctional.h>

#include <cmtkAffineMultiChannelRegistrationFunctional.h>
#include <cmtkSplineWarpMultiChannelRegistrationFunctional.h>
#include <cmtkSplineWarpMultiChannelIntensityCorrectionRegistrationFunctional.h>

#include <cmtkBestDirectionOptimizer.h>
#include <cmtkRegistrationCallback.h>

#include <cmtkClassStreamMultiChannelRegistration.h>

#include <list>
#include <queue>
#include <algorithm>

#ifdef HAVE_STDINT_H
#  include <stdint.h>
#else
typedef long int uint64_t;
#endif

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace mcwarp
{
#endif
bool verbose = false;

std::list<const char*> fileListRef;
std::list<const char*> fileListFlt;

const char* mcaffineOutput = NULL;
const char* outArchive = NULL;

std::list<cmtk::UniformVolume::SmartPtr> refChannelList;
std::list<cmtk::UniformVolume::SmartPtr> fltChannelList;

cmtk::Types::Coordinate gridSpacing = 40;
bool exactGridSpacing = false;
int gridRefinements = 0;
bool delayGridRefine = false;

bool adaptiveFixEntropy = false;
float adaptiveFixThreshFactor = 0.0;

bool fixWarpX = false;
bool fixWarpY = false;
bool fixWarpZ = false;

float jacobianConstraintWeight = 0.0;

cmtk::Types::Coordinate initialStepSize = 1.0;
cmtk::Types::Coordinate finalStepSize = 0.125;

bool metricNMI = false;
bool intensityCorrection = false;
bool useHistograms = true;
bool useCubicInterpolation = false;

int downsampleFrom = 1;
int downsampleTo = 1;
bool downsampleWithAverage = false;

float smoothSigmaFactor = 0.0;

cmtk::UniformVolume::SmartPtr
MakeDownsampled( const cmtk::UniformVolume::SmartPtr& image, const int downsample, const cmtk::Types::Coordinate smoothSigmaFactor )
{
  if ( downsampleWithAverage )
    return cmtk::UniformVolume::SmartPtr( new cmtk::UniformVolume( *image, downsample * image->GetMinDelta() ) );

  cmtk::UniformVolume::SmartPtr result( image->CloneGrid() );

  if ( (smoothSigmaFactor > 0) && downsample )
    {
    const cmtk::Types::Coordinate sigma = smoothSigmaFactor * downsample * image->GetMinDelta();
    cmtk::TypedArray::SmartPtr data( image->GetDataGaussFiltered( sigma ) );
    result->SetData( data );
    }
  else
    {
    result->SetData( image->GetData() );
    }

  if ( downsample > 1 )
    result = cmtk::UniformVolume::SmartPtr( result->GetDownsampled( downsample, true /*approxIsotropic*/ ) );
  return result;
}

template<class TMetricFunctional>
void
DoRegistration( int argc, char* argv[] )
{
  typedef cmtk::SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional> FunctionalType;
  typedef cmtk::SplineWarpMultiChannelIntensityCorrectionRegistrationFunctional<TMetricFunctional> ICFunctionalType;
  typename FunctionalType::SmartPtr functional;
  if ( intensityCorrection )
    {
    functional = typename FunctionalType::SmartPtr( new ICFunctionalType );
    }
  else
    {
    functional = typename FunctionalType::SmartPtr( new FunctionalType );
    }
  functional->SetNormalizedMI( metricNMI );  
  functional->SetAdaptiveFixEntropyThreshold( adaptiveFixEntropy );
  functional->SetAdaptiveFixThreshFactor( adaptiveFixThreshFactor );
  functional->SetJacobianConstraintWeight( jacobianConstraintWeight );

  if ( fixWarpX ) functional->AddFixedCoordinateDimension( 0 );
  if ( fixWarpY ) functional->AddFixedCoordinateDimension( 1 );
  if ( fixWarpZ ) functional->AddFixedCoordinateDimension( 2 );

  if ( mcaffineOutput )
    {
    cmtk::AffineMultiChannelRegistrationFunctional<TMetricFunctional> affineFunctional;
    cmtk::ClassStream inStream( mcaffineOutput, cmtk::ClassStream::READ );
    if ( !inStream.IsValid() )
      {
      cmtk::StdErr << "ERROR: could not open '" << mcaffineOutput << "' for reading.\n";
      exit( 1 );
      }
    inStream >> affineFunctional;
    inStream.Close();

    functional->SetInitialAffineTransformation( affineFunctional.GetTransformation() );
    for ( size_t idx = 0; idx < affineFunctional.GetNumberOfReferenceChannels(); ++idx )
      {
      cmtk::UniformVolume::SmartPtr channel = affineFunctional.GetReferenceChannel(idx);
      refChannelList.push_back( channel );
      }
    for ( size_t idx = 0; idx < affineFunctional.GetNumberOfFloatingChannels(); ++idx )
      {
      cmtk::UniformVolume::SmartPtr channel = affineFunctional.GetFloatingChannel(idx);
      fltChannelList.push_back( channel );
      }
    }

  const cmtk::Types::Coordinate* referenceImageDomain = NULL;
  cmtk::Types::Coordinate minPixelSize = FLT_MAX;
  for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it = refChannelList.begin(); it != refChannelList.end(); ++it )
    {
    if ( !referenceImageDomain )
      referenceImageDomain = (*it)->Size;

    minPixelSize = std::min( (*it)->GetMinDelta(), minPixelSize );
    }
  for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it = fltChannelList.begin(); it != fltChannelList.end(); ++it )
    {
    minPixelSize = std::min( (*it)->GetMinDelta(), minPixelSize );
    }
  
  cmtk::BestDirectionOptimizer optimizer;
  optimizer.SetCallback( cmtk::RegistrationCallback::SmartPtr( new cmtk::RegistrationCallback ) );
  optimizer.SetFunctional( functional );

  int gridRefinementsLeft = gridRefinements;
  std::queue<int> refinementSchedule;
  refinementSchedule.push( -1 ); // init marker
  for ( int downsample = std::max(downsampleFrom, downsampleTo); (downsample >= std::min(downsampleFrom,downsampleTo)) || gridRefinementsLeft; )
    {
    int refinementFlags = 0;
    if ( gridRefinementsLeft )
      {
      refinementFlags |= 2;
      --gridRefinementsLeft;
      }

    if ( downsample > downsampleTo )
      {
      refinementFlags |= 1;
      --downsample;
      }

    if ( ! refinementFlags )
      break;

    if ( delayGridRefine && (refinementFlags == 3) )
      {
      refinementSchedule.push( 1 );
      refinementSchedule.push( 2 );
      }
    else
      {
      refinementSchedule.push( refinementFlags );
      }
    }
  refinementSchedule.push( 0 ); // last iteration marker

  cmtk::CoordinateVector params;
  for ( int downsample = downsampleFrom; !refinementSchedule.empty(); refinementSchedule.pop() )
    {
    if ( verbose )
      cmtk::StdErr.printf( "Downsampling stage 1:%d\n", downsample );

    functional->ClearAllChannels();
    if ( (downsample == 0) || ( (downsample==1) && (smoothSigmaFactor==0) ) )
      {
      functional->AddReferenceChannels( refChannelList.begin(), refChannelList.end() );
      functional->AddFloatingChannels( fltChannelList.begin(), fltChannelList.end() );
      }
    else
      {
      for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it = refChannelList.begin(); it != refChannelList.end(); ++it )
	{
	cmtk::UniformVolume::SmartPtr image = MakeDownsampled( (*it), downsample, smoothSigmaFactor );
	image->m_MetaInformation[CMTK_META_FS_PATH] = (*it)->m_MetaInformation[CMTK_META_FS_PATH];
	functional->AddReferenceChannel( image );
	}

      for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it = fltChannelList.begin(); it != fltChannelList.end(); ++it )
	{
	cmtk::UniformVolume::SmartPtr image = MakeDownsampled( (*it), downsample, smoothSigmaFactor );
	image->m_MetaInformation[CMTK_META_FS_PATH] = (*it)->m_MetaInformation[CMTK_META_FS_PATH];
	functional->AddFloatingChannel( image );
	}
      }

    if ( refinementSchedule.front() == -1 )
      {
      functional->InitTransformation( referenceImageDomain, gridSpacing, exactGridSpacing );
      functional->GetParamVector( params );
      refinementSchedule.pop();
      }
    
    if ( verbose )
      cmtk::StdErr.printf( "Number of parameters is %d\n", functional->VariableParamVectorDim() );
    
    if ( optimizer.Optimize( params, initialStepSize * downsample * minPixelSize, finalStepSize * downsample * minPixelSize )
	 != cmtk::CALLBACK_OK )
      break;

    if ( refinementSchedule.front() & 1 )
      {
      --downsample;
      }

    if ( refinementSchedule.front() & 2 )
      {
      functional->RefineTransformation();
      functional->GetParamVector( params );
      }
    }
  
  if ( outArchive )
    {
    cmtk::ClassStream stream( outArchive, cmtk::ClassStream::WRITE );
    stream << *functional;
    stream.Close();
    }
}

int
main( int argc, char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Multi-channel nonrigid registration" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Multi-channel nonrigid B-spline image registration using histogram-based or covariance-based joint entropy measures" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] mcaffineOutput" );    
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Image Registration" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'o', "out-archive" ), &outArchive, "Output archive path." );

    cl.AddOption( Key( 'd', "downsample-from" ), &downsampleFrom, "Initial downsampling factor [1]." );
    cl.AddOption( Key( 'D', "downsample-to" ), &downsampleTo, "Final downsampling factor [1]." );
    cl.AddOption( Key( "downsample-average" ), &downsampleWithAverage, "Downsample using sliding-window averaging [default: off] )" );
    cl.AddOption( Key( "smooth" ), &smoothSigmaFactor, "Sigma of Gaussian smoothing kernel in multiples of template image pixel size [default: off] )" );

    cl.AddOption( Key( "grid-spacing" ), &gridSpacing, "Initial control point grid spacing in mm." );
    cl.AddOption( Key( "grid-spacing-exact" ), &gridSpacing, "Exact initial grid spacing, even if it doesn't match image FOV.", &exactGridSpacing );
    cl.AddOption( Key( "refine-grid" ), &gridRefinements, "Number of control point grid refinements." );
    cl.AddSwitch( Key( "delay-refine-grid" ), &delayGridRefine, true, "Delay control point grid refinement until after pixel refinement." );
    cl.AddOption( Key( "adaptive-fix-thresh-factor" ), &adaptiveFixThreshFactor, "Intensity threshold factor [0..1] for adaptive parameter fixing. [default: 0 -- no fixing]" );
    cl.AddOption( Key( "adaptive-fix-thresh-factor-entropy" ), &adaptiveFixThreshFactor, "Entropy threshold factor [0..1] for adaptive parameter fixing. [default: 0 -- no fixing]", &adaptiveFixEntropy );

    cl.AddSwitch( Key( "fix-warp-x" ), &fixWarpX, true, "Fix transformation (do not warp) in x-direction." );
    cl.AddSwitch( Key( "fix-warp-y" ), &fixWarpY, true, "Fix transformation (do not warp) in y-direction." );
    cl.AddSwitch( Key( "fix-warp-z" ), &fixWarpZ, true, "Fix transformation (do not warp) in z-direction." );

    cl.AddOption( Key( "jacobian-constraint-weight" ), &jacobianConstraintWeight, "Weight for Jacobian volume preservation constraint" );

    cl.AddSwitch( Key( 'H', "histograms" ), &useHistograms, true, "Use multi-dimensional histograms to compute entropies [default." );
    cl.AddSwitch( Key( 'C', "covariance" ), &useHistograms, false, "Use covariance matrix determinants to compute entropies." );
    cl.AddSwitch( Key( 'c', "cubic" ), &useCubicInterpolation, true, "Use cubic interpolation [default: linear]" );

    cl.AddSwitch( Key( "mi" ), &metricNMI, false, "Use standard mutual information metric [default]" );
    cl.AddSwitch( Key( "nmi" ), &metricNMI, true, "Use normalized mutual information metric" );
    cl.AddSwitch( Key( 'I', "intensity-correction" ), &intensityCorrection, true, "Correct image intensities using local transformation Jacobian to preserve total signal" );
    
    cl.AddOption( Key( "initial-step-size" ), &initialStepSize, "Initial optimizer step size in pixels." );
    cl.AddOption( Key( "final-step-size" ), &finalStepSize, "Initial optimizer step size in pixels." );

    cl.Parse();

    mcaffineOutput = cl.GetNext();
    }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  if ( useCubicInterpolation )
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Cubic> InterpolatorType;
    if ( useHistograms )
      {
      typedef cmtk::MultiChannelHistogramRegistrationFunctional<float,InterpolatorType,uint64_t,6> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>( argc, argv );
      }
    else
      {
      typedef cmtk::MultiChannelRMIRegistrationFunctional<float> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>( argc, argv );
      }
    }
  else
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> InterpolatorType;
    if ( useHistograms )
      {
      typedef cmtk::MultiChannelHistogramRegistrationFunctional<float,InterpolatorType,uint64_t,6> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>( argc, argv );
      }
    else
      {
      typedef cmtk::MultiChannelRMIRegistrationFunctional<float> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>( argc, argv );
      }
    }

  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace mcwarp
} // namespace apps
} // namespace cmtk
#endif

