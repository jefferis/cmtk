/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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
#include <System/cmtkProgressConsole.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkHistogram.h>
#include <Base/cmtkMathUtil.h>
#include <Base/cmtkMathFunctionWrappers.h>

#include <IO/cmtkVolumeIO.h>

#include <math.h>
#include <list>
#include <algorithm>

const char* OutputFileName = "average.nii";

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

cmtk::ScalarDataType DataType = cmtk::TYPE_FLOAT;

std::list<const char*> imagePathList;

void
GetNormalizationCoefficients
( const cmtk::TypedArray* floatingData, const cmtk::Types::DataItem refMean, const cmtk::Types::DataItem refVariance, cmtk::Types::DataItem& scale, cmtk::Types::DataItem& offset )
{
  cmtk::Types::DataItem fltMean, fltVariance;
  floatingData->GetStatistics( fltMean, fltVariance );
  scale = sqrt( refVariance ) / sqrt( fltVariance );
  offset = refMean - scale * fltMean;

  cmtk::DebugOutput( 1 ).GetStream().printf( "Converted grey values: %f +- %f -> %f +- %f\n", fltMean, sqrt( fltVariance ), refMean, sqrt( refVariance ) );
}

int
doMain( const int argc, const char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Average images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] image0 ..." );

    typedef cmtk::CommandLine::Key Key;

    cmtk::CommandLine::EnumGroup<ModeEnum>::SmartPtr modeGroup = cl.AddEnum( "mode", &Mode, "Mode of averaging operation" );
    modeGroup->AddSwitch( Key( "avg" ), MODE_AVG, "Compute average (i.e., mean) image" );
    modeGroup->AddSwitch( Key( "var" ), MODE_VAR, "Compute variance image" );
    modeGroup->AddSwitch( Key( "stdev" ), MODE_STDEV, "Compute standard deviation image" );
    modeGroup->AddSwitch( Key( "zscore" ), MODE_ZSCORE, "Compute z-score image" );
    modeGroup->AddSwitch( Key( "entropy" ), MODE_ENTROPY, "Compute pixel-by-pixel population entropy image" );

    cl.BeginGroup( "Preprocessing", "Data Preprocessing" );
    cl.AddSwitch( Key( 'l', "log" ), &ApplyLog, true, "Apply log to input data" );
    cl.AddSwitch( Key( 'a', "abs" ), &ApplyAbs, true, "Use absolute input values" );
    cl.AddSwitch( Key( 'n', "normalize-mean-stdev" ), &Normalize, true, "Normalize image intensities using means and standard deviations" );
    cl.AddOption( Key( "set-padding-value" ), &PaddingValue, "Define padding value in input images", &Padding );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'o', "outfile-name" ), &OutputFileName, "Output file name" );

    cmtk::CommandLine::EnumGroup<cmtk::ScalarDataType>::SmartPtr typeGroup = cl.AddEnum( "type", &DataType, "Scalar data type of output image." );
    typeGroup->AddSwitch( Key( "float" ), cmtk::TYPE_FLOAT, "Single-precision float." );
    typeGroup->AddSwitch( Key( "double" ), cmtk::TYPE_DOUBLE, "Double-precision float." );
    cl.EndGroup();

    cl.Parse( argc, argv );

    const char* next = cl.GetNext();
    while ( next )
      {
      imagePathList.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException(1);
    }

  cmtk::UniformVolume::SmartPtr volume( NULL );
  cmtk::TypedArray::SmartPtr outputData( NULL );

  bool firstImage = true;
  cmtk::Types::DataItem refMean, refVariance;

  cmtk::Types::DataItemRange imagesValueRange( 0, 0 );

  std::list<cmtk::TypedArray::SmartPtr> dataList;

  std::list<const char*>::const_iterator it;
  for ( it = imagePathList.begin(); it != imagePathList.end(); ++it ) 
    {
    cmtk::UniformVolume::SmartPtr nextVolume( cmtk::VolumeIO::ReadOriented( *it) );
    
    if ( ! nextVolume ) 
      {
      cmtk::StdErr << "ERROR: Could not open image " << *it << "\n";
      throw cmtk::ExitException( 1 );
      }
    
    if ( ! volume )
      {
      volume = nextVolume;
      }
    
    cmtk::TypedArray::SmartPtr data = nextVolume->GetData();
    data->Convert( DataType );

    if ( Padding ) 
      {
      data->SetPaddingValue( PaddingValue );
      }

    if ( ApplyLog )
      {
      data->ApplyFunctionDouble( cmtk::Wrappers::Log );
      }

    if ( ApplyAbs )
      {
      data->ApplyFunctionDouble( cmtk::Wrappers::Abs );
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
    
    const cmtk::Types::DataItemRange dataRange = data->GetRange();
    if ( firstImage )
      {
      imagesValueRange = dataRange;
      }
    else
      {
      imagesValueRange.m_LowerBound = std::min( imagesValueRange.m_LowerBound, dataRange.m_LowerBound );
      imagesValueRange.m_UpperBound = std::min( imagesValueRange.m_UpperBound, dataRange.m_UpperBound );
      }
    
    dataList.push_back( data );
    }

  // this is only used in "Entropy" mode, but we'll instantiate it anyway to save time
  cmtk::Histogram<float> histogram( NumberHistogramBins );
  histogram.SetRange( imagesValueRange );
  
  if ( ! outputData ) 
    {
    outputData = cmtk::TypedArray::SmartPtr( cmtk::TypedArray::Create( DataType, volume->GetNumberOfPixels() ) );
    outputData->SetPaddingValue( cmtk::MathUtil::GetFloatNaN() );
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
      }
      
    if ( actualSize )
      {
      pixelData.resize( actualSize );
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

  cmtk::Progress::Done();

  volume->SetData( outputData );
  cmtk::VolumeIO::Write( *volume, OutputFileName );
  
  return 0;
}

#include "cmtkSafeMain"
