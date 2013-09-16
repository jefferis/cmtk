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
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>

#include <Base/cmtkValueSequence.h>

#include <Registration/cmtkTypedArraySimilarity.h>

#include <IO/cmtkVolumeIO.h>

#include <vector>
#include <string>

int
doMain
( const int argc, const char *argv[] )
{
  std::vector<std::string> imagePaths;

  int sliceAxis = cmtk::AXIS_Z;

  cmtk::Types::DataItem standardDeviationsThreshold = 3.0;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Find bad slices in a time series of interleaved images (e.g., a resting-state fMRI series)." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool reads a time series of 3D images and detects outliers." );
    
    cl.AddParameterVector( &imagePaths, "ImagePaths", "List of file system paths for all images in the time series.");

    typedef cmtk::CommandLine::Key Key;

    cl.BeginGroup( "Input", "Input Options" );
    cmtk::CommandLine::EnumGroup<int>::SmartPtr sliceGroup = cl.AddEnum( "slice-orientation", &sliceAxis, "Define slice orientation of the diffusion images." );
    sliceGroup->AddSwitch( Key( "axial" ), (int)cmtk::AXIS_Z, "Axial slices" );
    sliceGroup->AddSwitch( Key( "sagittal" ),(int)cmtk::AXIS_X, "Sagittal slices" );
    sliceGroup->AddSwitch( Key( "coronal" ), (int)cmtk::AXIS_Y, "Coronal slices" );
    sliceGroup->AddSwitch( Key( "slice-x" ), (int)cmtk::AXIS_X, "X coordinate axis is slice direction" );
    sliceGroup->AddSwitch( Key( "slice-y" ), (int)cmtk::AXIS_Y, "Y coordinate axis is slice direction" );
    sliceGroup->AddSwitch( Key( "slice-z" ), (int)cmtk::AXIS_Z, "Z coordinate axis is slice direction" );
    cl.EndGroup();

    cl.BeginGroup( "detection", "Bad Slice Detection" );
    cl.AddOption( Key( "stdev-thresh" ), &standardDeviationsThreshold, "Threshold for bad slice identification in units of intensity standard deviations over all corresponding slices from the remaining diffusion images." );
    cl.EndGroup();

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  // read all diffusion images and make sure their grids match
  std::vector<cmtk::UniformVolume::SmartConstPtr> images( imagePaths.size() );
  for ( size_t i = 0; i < imagePaths.size(); ++i )
    {
    images[i] = cmtk::UniformVolume::SmartConstPtr( cmtk::VolumeIO::Read( imagePaths[i] ) );

    if ( i && ! images[0]->GridMatches( *images[i] ) )
      {
      cmtk::StdErr << "ERROR: geometry of image '" << imagePaths[i] << "' does not match that of image '" << imagePaths[0] << "'\n";
      throw cmtk::ExitException( 1 );
      }
    }
  
  // Build slice pair difference tables and statistics
  std::vector< cmtk::ValueSequence<double> > slicePairStatistics( images[0]->m_Dims[sliceAxis]-1 );
  std::vector< std::vector< double> > slicePairSamples( images[0]->m_Dims[sliceAxis]-1 );
  for ( int slice = 0; slice < images[0]->m_Dims[sliceAxis]-1; ++slice )
    {
    slicePairSamples[slice].resize( images.size() );
    for ( size_t i = 0; i < images.size(); ++i )
      {    
      const double difference = cmtk::TypedArraySimilarity::GetCrossCorrelation( images[i]->ExtractSlice( sliceAxis, slice )->GetData(), images[i]->ExtractSlice( sliceAxis, slice+1 )->GetData() );

      slicePairSamples[slice][i] = difference;
      slicePairStatistics[slice].Proceed( difference );
      }
    }
  
  // Search for outliers
  std::vector<double> means( images[0]->m_Dims[sliceAxis]-1 );
  std::vector<double> sdevs( images[0]->m_Dims[sliceAxis]-1 );
  for ( int slice = 0; slice < images[0]->m_Dims[sliceAxis]-1; ++slice )
    {
    means[slice] = slicePairStatistics[slice].GetAverage();
    sdevs[slice] = sqrt( slicePairStatistics[slice].GetVariance() );
    }

  for ( size_t i = 0; i < images.size(); ++i )
    {    
    for ( int slice = 0; slice < images[0]->m_Dims[sliceAxis]-1; ++slice )
      {
      if ( fabs( slicePairSamples[slice][i]-means[slice] ) / sdevs[slice] > standardDeviationsThreshold )
	cmtk::StdOut << imagePaths[i] << "\t" << slice << "\n";
      }
    }
    
  return 0;
}

#include "cmtkSafeMain"
