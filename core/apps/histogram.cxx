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
#include <cmtkHistogram.h>
#include <cmtkVolume.h>
#include <cmtkTypedArray.h>

#include <iostream>
#include <fstream>
#include <list>

#include <stdio.h>

bool Verbose = false;

const char *maskFile = NULL;
const char *outFile = NULL;

std::list<const char*> inFileList;

size_t NumberOfBins = 0;

bool NormalizeBins = false;

bool UserMinMaxValue = false;
cmtk::Types::DataItem UserMinValue = 0;
cmtk::Types::DataItem UserMaxValue = 0;
bool Truncate = false;

bool PaddingFlag = false;
bool PaddingValue = 0;

void
ComputeHistogram
( std::list<cmtk::TypedArray::SmartPtr>& dataList, cmtk::Histogram<double>& histogram )
{
  cmtk::Types::DataItem rangeFrom, rangeTo;
  histogram.GetRange( rangeFrom, rangeTo );

  std::list<cmtk::TypedArray::SmartPtr>::const_iterator it = dataList.begin();
  for ( ; it != dataList.end(); ++it )
    {
    const cmtk::TypedArray* data = *it;
    const size_t size = data->GetDataSize();
    
    cmtk::Types::DataItem value;
    for ( size_t offset = 0; offset < size; ++offset )
      if ( data->Get( value, offset ) )
	{
	if ( !Truncate || ((value >= rangeFrom) && (value <= rangeTo)) )
	  histogram.Increment( histogram.ValueToBin( value ) );
	}
    }
  
  if ( NormalizeBins )
    {
    histogram.NormalizeMaximum( 1.0 );
    }
}

void
ComputeHistogram
( std::list<cmtk::TypedArray::SmartPtr>& dataList, const cmtk::TypedArray* mask, const int maskValue, cmtk::Histogram<double>& histogram )
{
  cmtk::Types::DataItem rangeFrom, rangeTo;
  histogram.GetRange( rangeFrom, rangeTo );

  std::list<cmtk::TypedArray::SmartPtr>::const_iterator it = dataList.begin();
  for ( ; it != dataList.end(); ++it )
    {
    const cmtk::TypedArray* data = *it;
    const size_t size = data->GetDataSize();
    
    cmtk::Types::DataItem value, mv;
    for ( size_t offset = 0; offset < size; ++offset )
      if ( data->Get( value, offset ) )
	if ( !Truncate || ((value >= rangeFrom) && (value <= rangeTo)) )
	  if ( mask->Get( mv, offset ) && ( mv == maskValue ) )
	    histogram.Increment( histogram.ValueToBin( value ) );
    }
  
  if ( NormalizeBins )
    {
    histogram.NormalizeMaximum( 1.0 );
    } 
}

void 
WriteHistogram( const cmtk::Histogram<double>& histogram, const char* outfile )
{
  cmtk::Types::DataItem cumulative = 0;
  if ( outfile ) 
    {
    std::ofstream stream( outfile );
    for ( unsigned int bin = 0; bin < histogram.GetNumBins(); ++bin ) 
      {
      cumulative += histogram.GetBin( bin );
      stream << bin << "\t"
	     << histogram.BinToValue( bin ) << "\t"
	     << histogram.GetBin( bin ) << "\t"
	     << cumulative << "\n";
      }
    } 
  else
    {
    for ( unsigned int bin = 0; bin < histogram.GetNumBins(); ++bin ) 
      {
      cumulative += histogram.GetBin( bin );
      std::cout << bin << "\t"
		<< histogram.BinToValue( bin ) << "\t"
		<< histogram.GetBin( bin ) << "\t"
		<< cumulative << "\n";
      }
    }
}

int
main ( const int argc, const char* argv[] ) 
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image histogram" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Create a histogram of image intensities and write as tab-separated text file to standard output" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Statistics and Modeling" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'n', "nbins" ), &NumberOfBins, "Number of histogram bins (default: number of grey levels)" );
    cl.AddSwitch( Key( 'N', "normalize" ), &NormalizeBins, true, "Normalize histogram to maximum 1.0" );

    cl.AddOption( Key( "min" ), &UserMinValue, "User-defined minimum value", &UserMinMaxValue );
    cl.AddOption( Key( "max" ), &UserMaxValue, "User-defined maximum value", &UserMinMaxValue );
    cl.AddSwitch( Key( 't', "truncate" ), &Truncate, true, "Truncate histogram (do not enter values outside range into first/last bin)" );
    cl.AddOption( Key( "pad" ), &PaddingValue, "Image padding value", &PaddingFlag );
    
    cl.AddOption( Key( 'o', "outfile" ), &outFile, "File name pattern for histograms" );
    cl.AddOption( Key( 'm', "mask" ), &maskFile, "File name for multi-valued mask file" );

    cl.Parse();

    const char* next = cl.GetNext();
    while ( next )
      {
      inFileList.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e;
    return 1;
    }
  
  std::list<cmtk::TypedArray::SmartPtr> dataList;

  cmtk::Types::DataItem minData = 0, maxData = 0;
  std::list<const char*>::const_iterator it = inFileList.begin();
  for ( ; it != inFileList.end(); ++it )
    {
    cmtk::Volume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( *it, Verbose ) );
    if ( ! volume ) 
      {
      cmtk::StdErr << "Cannot read image " << *it << "\n";
      exit( 1 );
      }

    cmtk::TypedArray::SmartPtr data = volume->GetData();
    if ( ! data ) 
      {
      cmtk::StdErr << "Cannot read pixel data for image " << *it << "\n";
      exit( 1 );
      }
    
    if ( PaddingFlag )
      {
      data->SetPaddingValue( PaddingValue );
      }
    dataList.push_back( data );
    
    cmtk::Types::DataItem min, max;
    data->GetRange( min, max );
    
    if ( it == inFileList.begin() )
      {
      minData = min;
      maxData = max;
      }
    else
      {
      minData = std::min( minData, min );
      maxData = std::max( maxData, max );
      }
    }
  
  size_t bins = NumberOfBins;
  if ( !bins )
    {
    bins = static_cast<int>( maxData - minData + 1 );
    if ( bins > 256 ) bins = 256;
    }
  
  cmtk::Histogram<double> histogram( bins );
  if ( UserMinMaxValue )
    {
    histogram.SetRange( UserMinValue, UserMaxValue );
    } 
  else
    {
    histogram.SetRange( minData, maxData );
    }
  
  if ( ! maskFile ) 
    {
    histogram.Reset();
    ComputeHistogram( dataList, histogram );
    WriteHistogram( histogram, outFile );
    } 
  else 
    {
    cmtk::Volume::SmartPtr mask( cmtk::VolumeIO::ReadOriented( maskFile, Verbose ) );

    if ( mask ) 
      {
      const cmtk::TypedArray* maskData = mask->GetData();
      
      bool labelFlags[256];
      memset( labelFlags, 0, sizeof( labelFlags ) );
      for ( size_t i = 0; i < maskData->GetDataSize(); ++i )
	{
	cmtk::Types::DataItem l;
	if ( maskData->Get( l, i ) )
	  labelFlags[static_cast<byte>( l )] = true;
	}
      
      char fname[PATH_MAX];
      for ( int mv = 0; mv < 256; ++mv )
	{
	if ( labelFlags[mv] )
	  {
	  histogram.Reset();
	  ComputeHistogram( dataList, mask->GetData(), mv, histogram );
	  if ( snprintf( fname, PATH_MAX, outFile, mv ) > PATH_MAX )
	    {
	    cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	    }
	  else
	    {
	    WriteHistogram( histogram, fname );
	    }
	  }
	}
      }
    }
}
