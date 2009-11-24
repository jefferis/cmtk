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
#include <cmtkTypedArray.h>
#include <cmtkSmartPtr.h>

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>
#include <cmtkMathFunctionWrappers.h>

#include <cmtkHypothesisTests.h>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace ttest
{
#endif
typedef enum { TTEST, TTEST_PAIRED, CORRELATION_PAIRED, ZSCORES } ModeEnum;

int
main ( int argc, char* argv[] ) 
{
  bool Verbose = false;

  bool Symmetric = false;

  ModeEnum Mode = TTEST;

  bool UseLogData = false;
  bool UseAbsData = false;

  bool ConvertToShort = false;
  bool AbsoluteOutput = false;
  bool Invert = false;

  bool TextFileMode = false;
  
  const char* OutFileName = "ttest.hdr";
  const char* TStatFileName = NULL;
  const char* MaskFileName = NULL;
  
  std::list<const char*> FileListX;
  std::list<const char*> FileListY;

  
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "T-tests" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Pixelwise tests of statistical significance. Also compute correlations and z-scores" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] imageX0 [imageX1 ...] -- imageY0 [imageY1] ...\n"
		       "[options] --symmetric imageX0 [imageX1 ...]");
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Statistics and Modeling" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddSwitch( Key( 'S', "symmetric" ), &Symmetric, true, "Use mirrored X data as Y data" );
    cl.AddSwitch( Key( 'l', "log" ), &UseLogData, true, "Use log data for testing" );
    cl.AddSwitch( Key( 'a', "abs" ), &UseAbsData, true, "Use absolute data for testing" );

    cl.AddSwitch( Key( 'A', "abs-out" ), &AbsoluteOutput, true, "Generate absolute (unsigned) data" );
    cl.AddSwitch( Key( 'i', "invert" ), &Invert, true, "Invert data, i.e., subtract from 1.0" );
    cl.AddSwitch( Key( 's', "short" ), &ConvertToShort, true, "Convert probabilities to short (and scale by 1000)" );

    cl.AddSwitch( Key( 'p', "paired" ), &Mode, TTEST_PAIRED, "Compute paired t-test [default: groupwise]" );
    cl.AddSwitch( Key( 'c', "cross-correlation" ), &Mode, CORRELATION_PAIRED, "Compute paired cross-correlation" );
    cl.AddSwitch( Key( 'Z', "zscores" ), &Mode, ZSCORES, "Compute z-scores of test (y) distribution averages" );

    cl.AddSwitch( Key( 't', "text" ), &TextFileMode, true, "Text file input and output rather than image files." );
    cl.AddOption( Key( 'm', "mask" ), &MaskFileName, "Mask file name" );
    cl.AddOption( Key( 'o', "outfile" ), &OutFileName, "Output file name" );
    cl.AddOption( Key( "tstats-file" ), &TStatFileName, "T-statistics output file name (for correlation: p-value output file)" );

    cl.Parse();

    const char* next = cl.GetNext();
    while ( next && strcmp( next, "--" ) ) 
      {
      FileListX.push_back( next );
      next = cl.GetNextOptional();
      }
    
    if ( ! Symmetric ) 
      {
      next = cl.GetNextOptional();
      while ( next ) 
	{
	FileListY.push_back( next );
	next = cl.GetNextOptional();
	}
      }    
  }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr refVolume;

  std::vector<cmtk::TypedArray::SmartPtr> dataX;
  std::vector<cmtk::TypedArray::SmartPtr> dataY;
  std::list<const char*>::const_iterator fnameIt;

  if ( TextFileMode )
    {
    fnameIt = FileListX.begin();
    for ( ; fnameIt != FileListX.end(); ++fnameIt ) 
      {
      if ( Verbose )
	{
	std::cerr << "Reading X data file " << *fnameIt << "...\n";
	}
      std::ifstream stream( *fnameIt );
      std::vector<cmtk::Types::DataItem> data;
      while ( ! stream.eof() )
	{
	std::string line;
	getline( stream, line );
	if ( line[0] != '#' )
	  {
	  float value;
	  if ( 1 == sscanf( line.c_str(), "%f", &value ) )
	    {
	    data.push_back( value );
	    }
	  }
	}
      
      cmtk::TypedArray::SmartPtr array
	( cmtk::TypedArray::Create( cmtk::TYPE_ITEM, data.size() ) );
      for ( size_t i = 0; i < data.size(); ++i )
	array->Set( data[i], i );

      dataX.push_back( array );
      }

    fnameIt = FileListY.begin();
    for ( ; fnameIt != FileListY.end(); ++fnameIt ) 
      {
      if ( Verbose )
	{
	std::cerr << "Reading Y data file " << *fnameIt << "...\n";
	}
      std::ifstream stream( *fnameIt );
      std::vector<cmtk::Types::DataItem> data;
      
      while ( ! stream.eof() )
	{
	std::string line;
	getline( stream, line );
	if ( line[0] != '#' )
	  {
	  float value;
	  if ( 1 == sscanf( line.c_str(), "%f", &value ) )
	    {
	    data.push_back( value );
	    }
	  }
	}

      cmtk::TypedArray::SmartPtr array( cmtk::TypedArray::Create( cmtk::TYPE_ITEM, data.size() ) );
      for ( size_t i = 0; i < data.size(); ++i )
	array->Set( data[i], i );

      dataY.push_back( array );
      }
    } 
  else
    {
    fnameIt = FileListX.begin();
    for ( ; fnameIt != FileListX.end(); ++fnameIt ) 
      {
      std::cerr << "Reading X volume " << *fnameIt << "...\n";
      cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( *fnameIt, Verbose ) );
      if ( volume ) 
	{
	if ( ! refVolume ) refVolume = volume;
	
	if ( volume->GetData() ) 
	  {
	  dataX.push_back( volume->GetData() );

	  if ( Symmetric )
	    {
	    cmtk::TypedArray::SmartPtr mirroredData
	      ( volume->GetDataMirrored( cmtk::AXIS_X ) );
	    dataY.push_back( mirroredData );
	    }
	  }
	}
      }
    
    if ( ! Symmetric )
      {
      fnameIt = FileListY.begin();
      for ( ; fnameIt != FileListY.end(); ++fnameIt ) 
	{
	std::cerr << "Reading Y volume " << *fnameIt << "...\n";
	cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( *fnameIt, Verbose ) );
	if ( volume ) 
	  {
	  if ( ! refVolume ) refVolume = volume;
	  
	  if ( volume->GetData() ) 
	    {
	    dataY.push_back( volume->GetData() );
	    }
	  }
	}
      }
    }
  
  std::vector<cmtk::TypedArray::SmartPtr>::iterator dataIt;
  if ( UseLogData )
    {
    for ( dataIt = dataX.begin(); dataIt != dataX.end(); ++dataIt )
      (*dataIt)->ApplyFunctionDouble( cmtk::Wrappers::Log );
    for ( dataIt = dataY.begin(); dataIt != dataY.end(); ++dataIt )
      (*dataIt)->ApplyFunctionDouble( cmtk::Wrappers::Log );
    }
  
  if ( UseAbsData )
    {
    for ( dataIt = dataX.begin(); dataIt != dataX.end(); ++dataIt )
      (*dataIt)->ApplyFunctionDouble( cmtk::Wrappers::Abs );
    for ( dataIt = dataY.begin(); dataIt != dataY.end(); ++dataIt )
      (*dataIt)->ApplyFunctionDouble( cmtk::Wrappers::Abs );
    }

  cmtk::TypedArray::SmartPtr maskData( NULL );
  if ( MaskFileName ) 
    {
    if ( TextFileMode )
      {
      cmtk::StdErr << "WARNING: mask image not supported in text file mode\n";
      }
    else
      {
      cmtk::UniformVolume::SmartPtr volume
	( cmtk::VolumeIO::ReadOriented( MaskFileName, Verbose ) );
      if ( volume ) 
	maskData = volume->GetData();
      if ( ! maskData )
	{
	cmtk::StdErr << "WARNING: could not read mask image " << MaskFileName << "\n";
	}
      }
    }
  
  if ( TextFileMode || refVolume ) 
    {
    cmtk::TypedArray::SmartPtr probData;
    
    switch ( Mode )
      {
      case TTEST:
      {
      // allocated by GetUnpairedTTest; freed by SP:
      cmtk::TypedArray *tstatsData = NULL, *avgXData = NULL, *avgYData = NULL;
      
      if ( TextFileMode )
	{
	if ( dataY.size() )
	  {
	  probData = cmtk::TypedArray::SmartPtr( cmtk::HypothesisTests::GetUnpairedTwoTailedTTest( dataX, dataY, &tstatsData, &avgXData, &avgYData, maskData ) );
	  }
	else
	  {
	  probData = cmtk::TypedArray::SmartPtr( cmtk::HypothesisTests::GetOneSampleTTest( dataX, &tstatsData, &avgXData, maskData ) );
	  }
	}
      else
	{
	if ( dataY.size() )
	  {
	  probData = cmtk::TypedArray::SmartPtr( cmtk::HypothesisTests::GetUnpairedTwoTailedTTest( dataX, dataY, &tstatsData, NULL /*avgXData*/, NULL /*avgYData*/, maskData ) );
	  }
	else
	  {
	  probData = cmtk::TypedArray::SmartPtr( cmtk::HypothesisTests::GetOneSampleTTest( dataX, &tstatsData, NULL /*avgXData*/, maskData ) );
	  }
	}
      
      if ( refVolume && TStatFileName )
	{	
	if ( Verbose )
	  {
	  cmtk::StdErr << "Writing T-statistics to file " << TStatFileName << "\n";
	  }
	
	refVolume->SetData( cmtk::TypedArray::SmartPtr( tstatsData ) );
	cmtk::VolumeIO::Write( refVolume, TStatFileName, Verbose );
	}

      if ( AbsoluteOutput ) probData->ApplyFunctionDouble( cmtk::Wrappers::Abs );
      if ( Invert ) probData->Rescale( -1.0, 1.0 );
      
      if ( ConvertToShort ) 
	{
	probData->Rescale( 1000, 0 );
	probData = cmtk::TypedArray::SmartPtr( probData->Convert( cmtk::TYPE_SHORT ) );
	}
      
      if ( TextFileMode )
	{
	fprintf( stdout, "#M\tp\t\tT\t\tavgX\t\tavgY\n" );
	for ( size_t i = 0; i < probData->GetDataSize(); ++i )
	  {
	  cmtk::Types::DataItem p = 0, t = 0, avgX = 0, avgY = 0;
	  probData->Get( p, i );
	  tstatsData->Get( t, i );
	  avgXData->Get( avgX, i );
	  if ( avgYData ) avgYData->Get( avgY, i );
	  fprintf( stdout, "%d\t%+-15g\t%+-15g\t%+-15g\t%+-15g\n", (int)i, p, t, avgX, avgY );
	  }
	}
      else
	{
	if ( Verbose )
	  {
	  cmtk::StdErr << "Writing T-probablities to file " << OutFileName << "\n";
	  }
	
	refVolume->SetData( probData );
	cmtk::VolumeIO::Write( refVolume, OutFileName );
	}
      break;
      }
      case TTEST_PAIRED:
      {
      cmtk::TypedArray *tstatsData = NULL;
      
      if ( dataY.size() )
	{
	probData = cmtk::TypedArray::SmartPtr( cmtk::HypothesisTests::GetPairedTwoTailedTTest( dataX, dataY, &tstatsData, NULL /*avgXData*/, NULL /*avgYData*/, maskData ) );
	}
      
      if ( refVolume && TStatFileName )
	{	
	if ( Verbose )
	  {
	  cmtk::StdErr << "Writing T-statistics to file " << TStatFileName << "\n";
	  }
	
	refVolume->SetData( cmtk::TypedArray::SmartPtr( tstatsData ) );
	cmtk::VolumeIO::Write( refVolume, TStatFileName, Verbose );
	}
      
      if ( AbsoluteOutput ) probData->ApplyFunctionDouble( cmtk::Wrappers::Abs );
      if ( Invert ) probData->Rescale( -1.0, 1.0 );
      
      if ( ConvertToShort ) 
	{
	probData->Rescale( 1000, 0 );
	probData = cmtk::TypedArray::SmartPtr( probData->Convert( cmtk::TYPE_SHORT ) );
	}
      
      if ( Verbose )
	{
	cmtk::StdErr << "Writing T-probablities to file " << OutFileName << "\n";
	}
      
      refVolume->SetData( probData );
      cmtk::VolumeIO::Write( refVolume, OutFileName );
      break;
      }
      case CORRELATION_PAIRED:
      {
      cmtk::TypedArray *pData = NULL;
      if ( dataY.size() )
	{
	probData = cmtk::TypedArray::SmartPtr( cmtk::HypothesisTests::GetPairedCorrelation( dataX, dataY, &pData, maskData ) );
	}
      
      if ( AbsoluteOutput ) probData->ApplyFunctionDouble( cmtk::Wrappers::Abs );
      if ( Invert ) probData->Rescale( -1.0, 1.0 );
      
      if ( ConvertToShort ) 
	{
	probData->Rescale( 1000, 0 );
	probData = cmtk::TypedArray::SmartPtr( probData->Convert( cmtk::TYPE_SHORT ) );
	}
      
      if ( refVolume && TStatFileName )
	{	
	if ( Verbose )
	  {
	  cmtk::StdErr << "Writing probabilities to file " << TStatFileName << "\n";
	  }
	
	refVolume->SetData( cmtk::TypedArray::SmartPtr( pData ) );
	cmtk::VolumeIO::Write( refVolume, TStatFileName, Verbose );
	}
      
      if ( Verbose )
	{
	cmtk::StdErr << "Writing correlations to file " << OutFileName << "\n";
	}
      
      refVolume->SetData( probData );
      cmtk::VolumeIO::Write( refVolume, OutFileName );
      break;
      }
      case ZSCORES:
      {
      cmtk::TypedArray::SmartPtr zscoreData = cmtk::TypedArray::SmartPtr( cmtk::HypothesisTests::GetZScores( dataX, dataY, maskData ) );
      
      if ( AbsoluteOutput ) zscoreData->ApplyFunctionDouble( cmtk::Wrappers::Abs );
      if ( Invert ) zscoreData->Rescale( -1.0, 1.0 );
      
      if ( ConvertToShort ) 
	{
	zscoreData->Rescale( 1000, 0 );
	zscoreData = cmtk::TypedArray::SmartPtr( zscoreData->Convert( cmtk::TYPE_SHORT ) );
	}

      if ( Verbose )
	{
	cmtk::StdErr << "Writing z-scores map to file " << OutFileName << "\n";
	}
      
      refVolume->SetData( zscoreData );
      cmtk::VolumeIO::Write( refVolume, OutFileName );
      }
      break;
      }
    }
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace model
} // namespace apps
} // namespace cmtk
#endif

