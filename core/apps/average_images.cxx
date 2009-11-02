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

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>
#include <cmtkProgressConsole.h>

#include <list>
#include <algorithm>

#include <cmtkVolumeIO.h>
#include <cmtkUniformVolume.h>
#include <cmtkTypedArray.h>
#include <cmtkHistogram.h>

#include <math.h>

bool Verbose = false;
const char* OutputFileName = "average.hdr";

bool ApplyLog = false;
bool ApplyAbs = false;
bool Normalize = false;

bool Padding = false;
cmtk::Types::DataItem PaddingValue = 0;

size_t NumberHistogramBins = 64;

typedef enum
{
  MODE_AVG,
  MODE_STDEV,
  MODE_VAR,
  MODE_ZSCORE,
  MODE_ENTROPY
} ModeEnum;

ModeEnum Mode = MODE_AVG;

std::list<const char*> imagePathList;

void
GetNormalizationCoefficients
( const cmtk::TypedArray* floatingData, const cmtk::Types::DataItem refMean, const cmtk::Types::DataItem refVariance, cmtk::Types::DataItem& scale, cmtk::Types::DataItem& offset )
{
  cmtk::Types::DataItem fltMean, fltVariance;
  floatingData->GetStatistics( fltMean, fltVariance );
  scale = sqrt( refVariance ) / sqrt( fltVariance );
  offset = refMean - scale * fltMean;

  if ( Verbose ) 
    {
    fprintf( stderr, "Converted grey values: %f +- %f -> %f +- %f\n", fltMean, sqrt( fltVariance ), refMean, sqrt( refVariance ) );
    }
}

int main( const int argc, const char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Average images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] image0 ..." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );
    cl.AddOption( Key( 'o', "outfile-name" ), &OutputFileName, "Output file name" );

    cl.AddOption( Key( "pad" ), &PaddingValue, "Define padding value in input images", &Padding );

    cl.AddSwitch( Key( 'l', "log" ), &ApplyLog, true, "Apply log to input data" );
    cl.AddSwitch( Key( 'a', "abs" ), &ApplyAbs, true, "Use absolute input values" );
    cl.AddSwitch( Key( 'n', "normalize" ), &Normalize, true, "Normalize image intensities" );

    cl.AddSwitch( Key( 'A', "avg" ), &Mode, MODE_AVG, "Output average image" );
    cl.AddSwitch( Key( 'V', "var" ), &Mode, MODE_VAR, "Output variance image" );
    cl.AddSwitch( Key( 'S', "stdev" ), &Mode, MODE_STDEV, "Output standard deviation image" );
    cl.AddSwitch( Key( 'Z', "zscore" ), &Mode, MODE_ZSCORE, "Output zscore image" );
    cl.AddSwitch( Key( 'E', "entropy" ), &Mode, MODE_ENTROPY, "Output pixel-by-pixel population entropy image" );

    cl.Parse();

    const char* next = cl.GetNext();
    while ( next )
      {
      imagePathList.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e << "\n";
    exit(1);
    }

  cmtk::UniformVolume::SmartPtr volume( NULL );
  cmtk::TypedArray::SmartPtr outputData( NULL );

  bool firstImage = true;
  cmtk::Types::DataItem refMean, refVariance;

  cmtk::Types::DataItem imagesMinValue = 0, imagesMaxValue = 0;

  std::list<cmtk::TypedArray::SmartPtr> dataList;

  std::list<const char*>::const_iterator it;
  for ( it = imagePathList.begin(); it != imagePathList.end(); ++it ) 
    {
    cmtk::UniformVolume::SmartPtr nextVolume( dynamic_cast<cmtk::UniformVolume*>( cmtk::VolumeIO::ReadOriented( *it, Verbose ) ) );
    
    if ( ! nextVolume ) 
      {
      cmtk::StdErr << "ERROR: Could not open image " << *it << "\n";
      exit( 1 );
      }
    
    if ( ! volume )
      {
      volume = nextVolume;
      }
    
    cmtk::TypedArray::SmartPtr data = nextVolume->GetData();

    if ( Padding ) 
      {
      data->SetPaddingValue( PaddingValue );
      }

    if ( ApplyLog )
      {
      data->ApplyFunctionDouble( log );
      }

    if ( ApplyAbs )
      {
      data->ApplyFunctionDouble( fabs );
      }

    if ( Normalize ) 
      {
      if ( firstImage )
	{
	firstImage = false;
	data->GetStatistics( refMean, refVariance );
	}
      else
	{
	cmtk::Types::DataItem normFactor, normOffset;
	GetNormalizationCoefficients( data, refMean, refVariance, normFactor, normOffset );
	data->Rescale( normFactor, normOffset );
	}
      }
    
    cmtk::Types::DataItem dataMin, dataMax;
    data->GetRange( dataMin, dataMax );
    if ( firstImage )
      {
      imagesMinValue = dataMin;
      imagesMaxValue = dataMax;
      }
    else
      {
      imagesMinValue = std::min( imagesMinValue, dataMin );
      imagesMaxValue = std::max( imagesMaxValue, dataMax );
      }
    
    dataList.push_back( data );
    }

  // this is only used in "Entropy" mode, but we'll instantiate it anyway to save time
  cmtk::Histogram<float> histogram( NumberHistogramBins, true /*reset*/ );
  histogram.SetRange( imagesMinValue, imagesMaxValue );
  
  if ( ! outputData ) 
    {
    outputData = cmtk::TypedArray::SmartPtr( cmtk::TypedArray::Create( cmtk::TYPE_FLOAT, volume->GetNumberOfPixels() ) );
    outputData->SetPaddingValue( CMTK_FLOAT_NAN );
    } 
  
  cmtk::ProgressConsole progressIndicator;
  const int pixelsPerPercent = volume->GetNumberOfPixels() / 100;
  cmtk::Progress::Begin( 0, 100, 1, "Image averaging" );

  std::vector<cmtk::Types::DataItem> pixelData( dataList.size() );
  for ( size_t i = 0; i < volume->GetNumberOfPixels(); ++i ) 
    {
    if ( !(i % pixelsPerPercent) ) 
      cmtk::Progress::SetProgress( i / pixelsPerPercent );

    pixelData.resize( dataList.size() );
    size_t actualSize = 0;
    std::list<cmtk::TypedArray::SmartPtr>::const_iterator dit;
    for ( dit = dataList.begin(); dit != dataList.end(); ++dit )
      {
      cmtk::Types::DataItem v;
      if ( (*dit)->Get( v, i ) ) 
	{
	pixelData[actualSize++] = v;
	}
      
      if ( actualSize )
	{
	pixelData.resize( actualSize );
	cmtk::Types::DataItem avg, var;
	avg = cmtk::MathUtil::Mean<cmtk::Types::DataItem>( pixelData );
	var = cmtk::MathUtil::Variance<cmtk::Types::DataItem>( pixelData, avg );
	switch ( Mode )
	  {
	  case MODE_AVG: 
	  {
	  const cmtk::Types::DataItem avg = cmtk::MathUtil::Mean<cmtk::Types::DataItem>( pixelData );
	  outputData->Set( avg, i );
	  break;
	  }
	  case MODE_VAR:
	  {
	  const cmtk::Types::DataItem avg = cmtk::MathUtil::Mean<cmtk::Types::DataItem>( pixelData );
          const cmtk::Types::DataItem var = cmtk::MathUtil::Variance<cmtk::Types::DataItem>( pixelData, avg );
	  outputData->Set( var, i );
	  break;
	  }
	  case MODE_STDEV:
	  {
          const cmtk::Types::DataItem avg = cmtk::MathUtil::Mean<cmtk::Types::DataItem>( pixelData );
          const cmtk::Types::DataItem var = cmtk::MathUtil::Variance<cmtk::Types::DataItem>( pixelData, avg );
	  outputData->Set( sqrt(var), i );
	  break;
	  }
	  case MODE_ENTROPY:
	  {
	  histogram.Reset();
	  for ( size_t idx = 0; idx < actualSize; ++idx )
	    histogram.IncrementFractional( histogram.ValueToBinFractional( pixelData[idx] ) );
	  
	  outputData->Set( histogram.GetEntropy(), i );
	  break;
	  }
	  case MODE_ZSCORE:
	  {
          const cmtk::Types::DataItem avg = cmtk::MathUtil::Mean<cmtk::Types::DataItem>( pixelData );
          const cmtk::Types::DataItem var = cmtk::MathUtil::Variance<cmtk::Types::DataItem>( pixelData, avg );
	  outputData->Set( avg / sqrt(var), i );
	  break;
	  }
	  }
	} 
      else
	{
	outputData->SetPaddingAt( i );
	}
      }
    }
  
  cmtk::Progress::Done();

  volume->SetData( outputData );
  cmtk::VolumeIO::Write( volume, OutputFileName, Verbose );
  
  return 0;
}
