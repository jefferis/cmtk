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

#include <cmtkAffineMultiChannelRegistrationFunctional.h>
#include <cmtkMultiChannelRMIRegistrationFunctional.h>
#include <cmtkMultiChannelHistogramRegistrationFunctional.h>

#include <cmtkBestNeighbourOptimizer.h>
#include <cmtkRegistrationCallback.h>

#include <cmtkClassStreamMultiChannelRegistration.h>

#include <list>
#include <algorithm>

#if defined(HAVE_STDINT_H)
#  include <stdint.h>
#else
typedef long int uint64_t;
#endif

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace mcaffine
{
#endif
bool verbose = false;

std::list<const char*> fileListRef;
std::list<const char*> fileListFlt;

const char* outArchive = NULL;

std::list<cmtk::UniformVolume::SmartPtr> refChannelList;
std::list<cmtk::UniformVolume::SmartPtr> fltChannelList;

std::list<int> numberDOFs;
const char*
CallbackNumberDOFs( const char* argv )
{
  const int dofs = atoi( argv );
  if ( (dofs == 3) || (dofs == 6) || (dofs == 7) || (dofs == 9) || (dofs == 12) )
    numberDOFs.push_back( dofs );
  else
    cmtk::StdErr.printf( "WARNING: number of DOFs cannot be %d (from '%s'). Ignoring this.\n", dofs, argv );

  return NULL;
}

cmtk::Types::Coordinate initialStepSize = 1.0;
cmtk::Types::Coordinate finalStepSize = 0.125;
bool alignCenters = true;
bool metricNMI = true;
bool useHistograms = true;
bool useCubicInterpolation = false;

int downsampleFrom = 1;
int downsampleTo = 1;
bool downsampleWithAverage = false;

float smoothSigmaFactor = 0.0;
cmtk::Types::Coordinate minPixelSize = FLT_MAX;

const char* cropReferenceFromIndex = NULL;
const char* cropReferenceToIndex = NULL;

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
DoRegistration() 
{
  typedef cmtk::AffineMultiChannelRegistrationFunctional<TMetricFunctional> FunctionalType;
  typename FunctionalType::SmartPtr functional( new FunctionalType );  
  functional->SetNormalizedMI( metricNMI );
  
  cmtk::BestNeighbourOptimizer optimizer;
  optimizer.SetCallback( cmtk::RegistrationCallback::SmartPtr( new cmtk::RegistrationCallback ) );
  optimizer.SetFunctional( functional );

  cmtk::CoordinateVector params;
  for ( int downsample = std::max(downsampleFrom, downsampleTo); downsample >= std::min(downsampleFrom, downsampleTo); --downsample )
    {
    if ( verbose )
      cmtk::StdErr.printf( "Downsampling stage 1:%d\n", downsample );

    functional->ClearAllChannels();
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

    if ( downsample == downsampleFrom )
      {
      functional->InitTransformation( alignCenters );
      functional->GetParamVector( params );
      }
    
    for ( std::list<int>::const_iterator itDOF = numberDOFs.begin(); itDOF != numberDOFs.end(); ++itDOF )
      {
      if ( verbose )
	cmtk::StdErr.printf( "Setting number of DOFs to %d\n", *itDOF );
      
      functional->SetNumberDOFs( *itDOF );
      const cmtk::Types::Coordinate effectiveMinPixelSize = std::max( 1, downsample ) * minPixelSize;
      optimizer.Optimize( params, initialStepSize * effectiveMinPixelSize, finalStepSize * effectiveMinPixelSize );
      }
    }

  if ( outArchive )
    {
    cmtk::ClassStream stream( outArchive, cmtk::ClassStream::WRITE );
    if ( stream.IsValid() )
      {
      stream << *functional;
      stream.Close();
      }
    else
      {
      cmtk::StdErr << "ERROR: could not open archive " << outArchive << " for writing.\n";
      }
    }
}

int
main( int argc, char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Multi-channel affine registration" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Multi-channel affine image registration using histogram-based or covariance-based joint entropy measures" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] refChannel0 [refChannel1 ...] -- fltChannel0 [fltChannel1 ...]" );    
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Image Registration" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'o', "out-archive" ), &outArchive, "Output archive path." );

    cl.AddOption( Key( 'd', "downsample-from" ), &downsampleFrom, "Initial downsampling factor [1]." );
    cl.AddOption( Key( 'D', "downsample-to" ), &downsampleTo, "Final downsampling factor [1]. Factor 0 is full resolution with smoothing turned off" );
    cl.AddOption( Key( "smooth" ), &smoothSigmaFactor, "Sigma of Gaussian smoothing kernel in multiples of template image pixel size [default: off] )" );
    cl.AddOption( Key( "downsample-average" ), &downsampleWithAverage, "Downsample using sliding-window averaging [default: off] )" );

    cl.AddCallback( Key( "dofs" ), CallbackNumberDOFs, "Add number of DOFs to optimization schedule [can be repeated]." );

    cl.AddSwitch( Key( "nmi" ), &metricNMI, true, "Use normalized mutual information metric [default]" );
    cl.AddSwitch( Key( "mi" ), &metricNMI, false, "Use standard mutual information metric" );

    cl.AddSwitch( Key( 'H', "histograms" ), &useHistograms, true, "Use multi-dimensional histograms to compute entropies [default]" );
    cl.AddSwitch( Key( 'C', "covariance" ), &useHistograms, false, "Use covariance matrix determinants to compute entropies" );
    cl.AddSwitch( Key( 'c', "cubic" ), &useCubicInterpolation, true, "Use cubic interpolation [default: linear]" );

    cl.AddOption( Key( "initial-step-size" ), &initialStepSize, "Initial optimizer step size in pixels." );
    cl.AddOption( Key( "final-step-size" ), &finalStepSize, "Initial optimizer step size in pixels." );

    cl.AddOption( Key( "crop-reference-from-index" ), &cropReferenceFromIndex, "Crop reference image from index x,y,z." );
    cl.AddOption( Key( "crop-reference-to-index" ), &cropReferenceToIndex, "Crop reference image to index x,y,z." );

    cl.Parse();

    const char* next = cl.GetNext();
    while ( next && strcmp( next, "--" ) ) 
      {
      fileListRef.push_back( next );
      next = cl.GetNextOptional();
      }
    
    next = cl.GetNext();
    while ( next ) 
      {
      fileListFlt.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  for ( std::list<const char*>::const_iterator refIt = fileListRef.begin(); refIt != fileListRef.end(); ++refIt )
    {
    cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( *refIt, verbose ) );
    if ( !volume || !volume->GetData() )
      {
      cmtk::StdErr << "ERROR: Cannot read image " << *refIt << "\n";
      exit( 1 );
      }
    minPixelSize = std::min( volume->GetMinDelta(), minPixelSize );
    refChannelList.push_back( volume );
    }

  if ( cropReferenceFromIndex )
    {
    int xyz[3];
    if ( 3 != sscanf( cropReferenceFromIndex, "%d,%d,%d", &xyz[0], &xyz[1], &xyz[2] ) )
      {
      cmtk::StdErr << "ERROR: reference crop from index could not parse index '" << cropReferenceFromIndex << "' as valid x,y,z index.\n";
      exit( 1 );
      }
    for ( std::list<cmtk::UniformVolume::SmartPtr>::iterator refIt = refChannelList.begin(); refIt != refChannelList.end(); ++refIt )
      {
      (*refIt)->SetCropRegionFrom( xyz );
      }
    }

  if ( cropReferenceToIndex )
    {
    int xyz[3];
    if ( 3 != sscanf( cropReferenceToIndex, "%d,%d,%d", &xyz[0], &xyz[1], &xyz[2] ) )
      {
      cmtk::StdErr << "ERROR: reference crop to index could not parse index '" << cropReferenceFromIndex << "' as valid x,y,z index.\n";
      exit( 1 );
      }
    for ( std::list<cmtk::UniformVolume::SmartPtr>::iterator refIt = refChannelList.begin(); refIt != refChannelList.end(); ++refIt )
      {
      (*refIt)->SetCropRegionTo( xyz );
      }
    }

  for ( std::list<const char*>::const_iterator fltIt = fileListFlt.begin(); fltIt != fileListFlt.end(); ++fltIt )
    {
    cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( *fltIt, verbose ) );
    if ( !volume || !volume->GetData() )
      {
      cmtk::StdErr << "ERROR: Cannot read image " << *fltIt << "\n";
      exit( 1 );
      }
    minPixelSize = std::min( volume->GetMinDelta(), minPixelSize );
    fltChannelList.push_back( volume );
    }

  if ( useCubicInterpolation )
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Cubic> InterpolatorType;
    if ( useHistograms )
      {
      typedef cmtk::MultiChannelHistogramRegistrationFunctional<float,InterpolatorType,uint64_t,6> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>();
      }
    else
      {
      typedef cmtk::MultiChannelRMIRegistrationFunctional<float> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>();
      }
    }
  else
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> InterpolatorType;
    if ( useHistograms )
      {
      typedef cmtk::MultiChannelHistogramRegistrationFunctional<float,InterpolatorType,uint64_t,6> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>();
      }
    else
      {
      typedef cmtk::MultiChannelRMIRegistrationFunctional<float> MetricFunctionalType;
      DoRegistration<MetricFunctionalType>();
      }
    }
  
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace mcaffine
} // namespace apps
} // namespace cmtk
#endif

