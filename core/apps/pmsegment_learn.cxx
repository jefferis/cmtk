/*
//
//  Copyright 2011 SRI International
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

#include <System/cmtkExitException.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkConsole.h>
#include <System/cmtkFileUtils.h>
#include <System/cmtkProgressConsole.h>

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUniformDistanceMap.h>

#include <Registration/cmtkTypedArraySimilarity.h>

#include <IO/cmtkVolumeIO.h>

#include <math.h>

#include <fstream>

const std::pair<size_t,float> 
GetLeastUnique( const std::vector< std::pair<float,cmtk::TypedArray::SmartPtr> >& patterns )
{
  const size_t nPatterns = patterns.size();

  std::vector<float> ncc( (nPatterns * (nPatterns+1)) / 2 );

  size_t ofs = 0;
  for ( size_t i = 0; i < nPatterns; ++i )
    {
    for ( size_t j = 0; j < i; ++j, ++ofs )
      {
      ncc[ofs] = cmtk::TypedArraySimilarity::GetCrossCorrelation( patterns[i].second, patterns[j].second );
      }
    }
  
  float maxValue = 0;
  size_t maxIndex = 0;
  
  for ( size_t i = 0; i < nPatterns; ++i )
    {
    float average = 0;
    for ( size_t j = 0; j < nPatterns; ++j )
      {
      if ( j < i )
	average += ncc[j + i*nPatterns];
      else if ( i < j )
	average += ncc[i + j*nPatterns];      
      }
    
    if ( average >= maxValue )
      {
      maxValue = average;
      maxIndex = i;
      }
    }
  
  return std::pair<size_t,float>( maxIndex, maxValue );
}

int regionRadius[3] = { 8, 8, 1 };
double pThreshold = 0.2;
int patternSetSize = 100;

std::vector< std::pair<float,cmtk::TypedArray::SmartPtr> > patternsPos;
std::vector< std::pair<float,cmtk::TypedArray::SmartPtr> > patternsNeg;

int
doMain( const int argc, const char* argv[] )
{
  const char* exampleImagePath = NULL;
  const char* exampleMaskPath = NULL;
  const char* modelDirectory = NULL;

  try
    {
    cmtk::CommandLine  cl( cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Learn pattern matching-based segmentation model" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Segmentation" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddParameter( &exampleImagePath, "ExampleImage", "Example intensity image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &exampleMaskPath, "ExampleMask", "Example binary mask image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_LABELS );
    cl.AddParameter( &modelDirectory, "OutputDirectory", "Segmentation model output directory path." )->SetProperties( cmtk::CommandLine::PROPS_DIRNAME | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }  

  cmtk::UniformVolume::SmartPtr exampleImage( cmtk::VolumeIO::Read( exampleImagePath ) );
  if ( ! exampleImage )
    {
    cmtk::StdErr << "ERROR: could not read example image " << exampleImagePath << "\n";
    return 1;
    }

  cmtk::UniformVolume::SmartConstPtr exampleMask( cmtk::VolumeIO::Read( exampleMaskPath ) );
  if ( ! exampleMask )
    {
    cmtk::StdErr << "ERROR: could not read example mask " << exampleMaskPath << "\n";
    return 1;
    }
  
  cmtk::FileUtils::RecursiveMkDir( modelDirectory );

  cmtk::UniformVolume::SmartConstPtr distanceMap( cmtk::UniformDistanceMap<float>( *exampleMask, cmtk::UniformDistanceMap<float>::SIGNED, 1 ).Get() );
  cmtk::DataGrid::RegionType insideRegion( exampleImage->GetWholeImageRegion() );

  const cmtk::DataGrid::IndexType radius = cmtk::DataGrid::IndexType( regionRadius );
  const cmtk::DataGrid::IndexType radiusPlus = cmtk::DataGrid::IndexType( regionRadius ).AddScalar(1);

  insideRegion.From() += radius;
  insideRegion.To() -= radius;

  cmtk::ProgressConsole progressIndicator;
  cmtk::Progress::Begin( insideRegion.From()[2], insideRegion.To()[2], 1, "Extracting example patterns" );
  
  cmtk::DataGrid::IndexType center;
  for ( center[2] = insideRegion.From()[2]; center[2] < insideRegion.To()[2]; ++center[2] )
    {
    cmtk::Progress::SetProgress( center[2] );      
    for ( center[1] = insideRegion.From()[1]; center[1] < insideRegion.To()[1]; ++center[1] )
      {
      for ( center[0] = insideRegion.From()[0]; center[0] < insideRegion.To()[0]; ++center[0] )
	{
	const size_t offset = distanceMap->GetOffsetFromIndex( center );

	const float distance = distanceMap->GetDataAt( offset );
	if ( (1.0-pThreshold) < exp( - distance*distance ) * cmtk::MathUtil::UniformRandom() )
	  {
	  exampleImage->SetCropRegion( cmtk::DataGrid::RegionType( center-radius, center+radiusPlus ) );
	  cmtk::TypedArray::SmartPtr cropData( exampleImage->GetCroppedData() );

	  const bool inside = exampleMask->GetDataAt( offset ) > 0;
	  std::pair<float,cmtk::TypedArray::SmartPtr> distanceAndData( distance, cropData );
	  if ( inside )
	    {
	    patternsPos.push_back( distanceAndData );
	    }
	  else
	    {
	    patternsNeg.push_back( distanceAndData );
	    }

	  if ( (patternsPos.size() + patternsNeg.size()) > patternSetSize )
	    {
	    const std::pair<size_t,float> leastUniquePos = GetLeastUnique( patternsPos );
	    const std::pair<size_t,float> leastUniqueNeg = GetLeastUnique( patternsNeg );

	    if ( leastUniquePos.second > leastUniqueNeg.second )
	      {
	      patternsPos[leastUniquePos.first] = *(patternsPos.rbegin());
	      patternsPos.pop_back();
	      }
	    else
	      {
	      patternsNeg[leastUniqueNeg.first] = *(patternsNeg.rbegin());
	      patternsNeg.pop_back();
	      }
	    }
	  }
	}
      }
    }
  cmtk::StdErr << "\n";

  cmtk::DebugOutput( 1 ) << "Collected " << patternsPos.size() << " positive and " << patternsNeg.size() << " negative patterns.\n";

  cmtk::Progress::Done();

  // write model
  char fname[32];
  char path[PATH_MAX];

  exampleImage->SetCropRegion( cmtk::DataGrid::RegionType( cmtk::DataGrid::IndexType( cmtk::DataGrid::IndexType::Init( 0 ) ), radius+radiusPlus ) );
  cmtk::UniformVolume::SmartPtr region( exampleImage->GetCroppedVolume() );

  snprintf( path, PATH_MAX, "%s%c%s", modelDirectory, CMTK_PATH_SEPARATOR, "info.txt" );
  std::ofstream infoStream( path );

  if ( infoStream.good() )
    {
    infoStream << regionRadius[0] << " " << regionRadius[1] << " " << regionRadius[2] << "\n";
    
    for ( size_t i = 0; i < patternsPos.size(); ++i )
      {
      snprintf( fname, 32, "pos%c%05zd.nii", CMTK_PATH_SEPARATOR, i );
      snprintf( path, PATH_MAX, "%s%c%s", modelDirectory, CMTK_PATH_SEPARATOR, fname );
      
      region->SetData( patternsPos[i].second );
      cmtk::VolumeIO::Write( *region, path );
      
      infoStream << fname << "\t" << patternsPos[i].first << "\n";
      }
    
    for ( size_t i = 0; i < patternsNeg.size(); ++i )
      {
      snprintf( fname, 32, "neg%c%05zd.nii", CMTK_PATH_SEPARATOR, i );
      snprintf( path, PATH_MAX, "%s%c%s", modelDirectory, CMTK_PATH_SEPARATOR, fname );
      
      region->SetData( patternsNeg[i].second );
      cmtk::VolumeIO::Write( *region, path );
      
      infoStream << fname << "\t" << patternsNeg[i].first << "\n";
      }
    }
  else
    {
    cmtk::StdErr << "ERROR: could not open model info file " << path << "\n";
    return 1;
    }
  
  return 0;
}

#include "cmtkSafeMain"
