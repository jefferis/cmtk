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

#include <cmtkTypedStreamStudylist.h>
#include <cmtkVolumeIO.h>
#include <cmtkSplineWarpXform.h>

#include <cmtkHistogram.h>
#include <cmtkValueSequence.h>

#include <cmtkClassStream.h>
#include <cmtkLandmark.h>
#include <cmtkLandmarkList.h>

#include <cmtkMathFunctionWrappers.h>

#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <list>

#include <cmtkReformatVolume.h>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace statistics
{
#endif
bool Verbose = false;

bool Label = false;
bool LogData = false;
bool ExpData = false;

bool WriteAsColumn = false;
bool OutputExpNotation = false;

const char *LandmarksFileName = NULL;
const char *MaskFileName = NULL;
bool MaskIsBinary = false;
std::list<const char*> ImageFileNames;

int NumberOfHistogramBins = 256;

std::list<cmtk::Types::DataItem> percentiles;

void
CallbackAddPercentile( const double arg )
{
  percentiles.push_back( static_cast<cmtk::Types::DataItem>( arg ) );
}

void
AnalyzeLabels( const cmtk::UniformVolume* volume, const cmtk::TypedArray* maskData ) 
{
  const cmtk::TypedArray* data = volume->GetData();
  cmtk::Types::DataItem min, max;
  data->GetRange( min, max );
  
  const unsigned int numberOfLabels = static_cast<unsigned int>( max - min + 1 );

  // Number of label voxels.
  std::vector<unsigned int> count( numberOfLabels );
  std::fill( count.begin(), count.end(), 0 );

  // Number of label surface voxels.
  std::vector<unsigned int> countSurface( numberOfLabels );
  std::fill( countSurface.begin(), countSurface.end(), 0 );

  // Centers-of-mass for each label
  std::vector<cmtk::Vector3D> centerOfMass( numberOfLabels );
  std::fill( centerOfMass.begin(), centerOfMass.end(), cmtk::Vector3D(0,0,0) );

  cmtk::Vector3D v;

  int index = 0;
  cmtk::Types::DataItem value, neighbor, maskValue;
  for ( int z = 0; z < volume->GetDims( cmtk::AXIS_Z ); ++z ) 
    {
    for ( int y = 0; y < volume->GetDims( cmtk::AXIS_Y ); ++y ) 
      {
      for ( int x = 0; x < volume->GetDims( cmtk::AXIS_X ); ++x, ++index ) 
	{
	if ( maskData && !(maskData->Get( maskValue, index ) && (maskValue != 0) ) )
	  continue;
	
	if ( data->Get( value, index ) ) 
	  {
	  const unsigned int labelIdx = static_cast<unsigned int>( value - min );
	  ++count[labelIdx];
	  volume->GetGridLocation( v, x, y, z );
	  centerOfMass[labelIdx] += v;
	  
	  bool isSurface = false;
	  for ( int dz = -1; (dz < 2) && !isSurface; ++dz )
	    for ( int dy = -1; (dy < 2) && !isSurface; ++dy )
	      for ( int dx = -1; (dx < 2) && !isSurface; ++dx )
		if ( dx || dy || dz )
		  if ( (dx+x)>=0 && (dx+x)<volume->GetDims( cmtk::AXIS_X ) && (dy+y)>=0 && (dy+y)<volume->GetDims( cmtk::AXIS_Y ) && (dz+z)>=0 && (dz+z)<volume->GetDims( cmtk::AXIS_Z ) ) 
		    {
		    int offset = (x+dx) + volume->GetDims( cmtk::AXIS_X ) * ( ( y+dy ) + volume->GetDims( cmtk::AXIS_Y ) * (z+dz) );
		    if ( data->Get( neighbor, offset ) && ( neighbor != value ) )
		      isSurface = true;
		    }
	  
	  if ( isSurface )
	    ++countSurface[labelIdx];
	  
	  }
	
	}
      }
    }
  
  cmtk::LandmarkList landmarkList;
  
  if ( Verbose )
    fputs( "idx\t\tcount\t\tsurface\t\tvolume\tCenterOfMass\n", stdout );

  const cmtk::Types::Coordinate voxelVolume = volume->m_Delta[0] * volume->m_Delta[1] * volume->m_Delta[2];

  size_t totalCount = 0;
  for ( unsigned int idx = 0; idx < numberOfLabels; ++idx ) 
    {
    if ( count[idx] ) 
      {
      centerOfMass[idx] *= (1.0 / count[idx]);
      if ( OutputExpNotation )
	fprintf( stdout, "%03d\t%14d\t%14d\t%.4e\t(%e,%e,%e)\n", 
		 (int)(idx + min), count[idx], countSurface[idx], count[idx] * voxelVolume, centerOfMass[idx].XYZ[0], centerOfMass[idx].XYZ[1], centerOfMass[idx].XYZ[2] );
      else
	fprintf( stdout, "%03d\t%14d\t%14d\t%.4f\t(%f,%f,%f)\n", 
		 (int)(idx + min), count[idx], countSurface[idx], count[idx] * voxelVolume, centerOfMass[idx].XYZ[0], centerOfMass[idx].XYZ[1], centerOfMass[idx].XYZ[2] );
      totalCount += count[idx];
      
      cmtk::Landmark* landmark = new cmtk::Landmark();
      char name[4];
      sprintf( name, "%3d", idx );
      landmark->SetName( name );
      landmark->SetLocation( centerOfMass[idx].XYZ );
      landmarkList.push_back( cmtk::Landmark::SmartPtr( landmark ) );
      }
    }
  
  if ( LandmarksFileName ) 
    {
    cmtk::ClassStream stream( LandmarksFileName, cmtk::ClassStream::WRITE );
    stream << landmarkList;
    stream.Close();
    }
  
  cmtk::Types::DataItem entropy = 0;
  if ( totalCount )
    {
    for ( unsigned int idx=0; idx < numberOfLabels; ++idx ) 
      {
      if ( count[idx] ) 
	{
	cmtk::Types::DataItem p = ((cmtk::Types::DataItem) count[idx]) / totalCount;
	entropy += p * log( p );
	}
      }
    }
  
  if ( OutputExpNotation )
    fprintf( stdout, "\nEntropy:\t%.5e\n", entropy );
  else
    fprintf( stdout, "\nEntropy:\t%.5f\n", entropy );
}

void
AnalyzeGrey( const cmtk::UniformVolume* volume, const cmtk::TypedArray* maskData ) 
{
  const cmtk::TypedArray* data = volume->GetData();
  cmtk::Types::DataItem min, max;
  data->GetRange( min, max );  

  cmtk::Histogram<unsigned int> histogram( NumberOfHistogramBins );
  histogram.SetRange( min, max );

  bool maskFlags[256];
  memset( maskFlags, 0, sizeof( maskFlags ) );

  for ( size_t i = 0; i < maskData->GetDataSize(); ++i )
    {
    cmtk::Types::DataItem l;
    if ( maskData->Get( l, i ) )
      maskFlags[static_cast<byte>( l )] = true;
    }

  if ( ! WriteAsColumn )
    fprintf( stdout, "#M\tmin\tmax\tmean\tsdev\tn\tH1\tH2\tsum\n" );
  
  for ( int maskSelect = 0; maskSelect < 256; ++maskSelect )
    {
    histogram.Reset();
    if ( ! maskFlags[maskSelect] ) continue;

    cmtk::Types::DataItem value, maskValue;
    cmtk::ValueSequence<cmtk::Types::DataItem> seq;
    
    size_t index = 0;
    for ( int z = 0; z < volume->GetDims( cmtk::AXIS_Z ); ++z ) 
      {
      for ( int y = 0; y < volume->GetDims( cmtk::AXIS_Y ); ++y ) 
	{
	for ( int x = 0; x < volume->GetDims( cmtk::AXIS_X ); ++x, ++index ) 
	  {
	  if ( maskData && maskData->Get( maskValue, index ) && (maskValue == maskSelect) )
	    {
	    if ( data->Get( value, index ) ) 
	      {
	      seq.Proceed( value );
	      histogram.Increment( histogram.ValueToBin( value ) );
	      }
	    }
	  }
	}
      }

    if ( seq.GetNValues() )
      {
      if ( ! WriteAsColumn )
	{
	if ( OutputExpNotation )
	  fprintf( stdout, "%03d\t%.5f\t%.5e\t%.5e\t%.5e\t%d\t%.5e\t%f\n",
		   maskSelect, seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(), sqrt( seq.GetVariance() ), seq.GetNValues(), histogram.GetEntropy(), seq.GetSum() );
	else
	  fprintf( stdout, "%03d\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%.5f\t%f\n",
		   maskSelect, seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(), sqrt( seq.GetVariance() ), seq.GetNValues(), histogram.GetEntropy(), seq.GetSum() );
	}
      }
    }

  // BUG: need to support percentiles also for masked operation
}

void
AnalyzeGrey( const cmtk::UniformVolume* volume ) 
{
  const cmtk::TypedArray* data = volume->GetData();
  cmtk::Types::DataItem min, max;
  data->GetRange( min, max );  

  if ( ! WriteAsColumn )
    fprintf( stdout, "min\tmax\tmean\tsdev\tn\tEntropy\tsum\n" );
  
  cmtk::Types::DataItem value;
  cmtk::ValueSequence<cmtk::Types::DataItem> seq;
  
  size_t index = 0;
  for ( int z = 0; z < volume->GetDims( cmtk::AXIS_Z ); ++z ) 
    {
    for ( int y = 0; y < volume->GetDims( cmtk::AXIS_Y ); ++y ) 
      {
      for ( int x = 0; x < volume->GetDims( cmtk::AXIS_X ); ++x, ++index ) 
	{
	if ( data->Get( value, index ) ) 
	  {
	  seq.Proceed( value );
	  }
	}
      }
    }

  if ( seq.GetNValues() )
    {
    if ( WriteAsColumn )
      {
      if ( OutputExpNotation )
	fprintf( stdout, "min\t%.5e\nmax\t%.5e\nmean\t%.5e\nsdev\t%.5e\nn\t%d\nEntropy\t%.5e\nsum\t%.5e\n",
		 seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(), 
		 sqrt( seq.GetVariance() ), seq.GetNValues(), data->GetEntropy( true /*fractional*/, NumberOfHistogramBins ), seq.GetSum());
      else
	fprintf( stdout, "min\t%.5f\nmax\t%.5f\nmean\t%.5f\nsdev\t%.5f\nn\t%d\nEntropy\t%.5f\nsum\t%f\n",
		 seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(), 
		 sqrt( seq.GetVariance() ), seq.GetNValues(), data->GetEntropy( true /*fractional*/, NumberOfHistogramBins ), seq.GetSum());
      }
    else
      {
      if ( OutputExpNotation )
	fprintf( stdout, "%.5e\t%.5e\t%.5e\t%.5e\t%d\t%.5e\t%e\n",
		 seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(), 
		 sqrt( seq.GetVariance() ), seq.GetNValues(), data->GetEntropy( true /*fractional*/, NumberOfHistogramBins ), seq.GetSum() );   
      else
	fprintf( stdout, "%.5f\t%.5f\t%.5f\t%.5f\t%d\t%.5f\t%f\n",
		 seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(), 
		 sqrt( seq.GetVariance() ), seq.GetNValues(), data->GetEntropy( true /*fractional*/, NumberOfHistogramBins ), seq.GetSum() );   
      }
    }
    
  if ( percentiles.size() )
    {
    fprintf( stdout, "PERC\tVALUE\n" );
    for ( std::list<cmtk::Types::DataItem>::const_iterator it = percentiles.begin(); it != percentiles.end(); ++it )
      {
      fprintf( stdout, "%.5f\t%.5f\n", (cmtk::Types::DataItem)(*it), (cmtk::Types::DataItem)data->GetPercentile( *it, NumberOfHistogramBins ) );
      }
    }
}

int
main ( int argc, char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image statistics" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Statistical computations on image pixel intensities, i.e., means and standard deviations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] ImageFile0 [ImageFile1 ...]" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Statistics and Modeling" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );
    
    cl.AddSwitch( Key( 'l', "label" ), &Label, true, "Interpret voxel values as labels" );
    cl.AddSwitch( Key( 'g', "gray" ), &Label, false, "Interpret voxel values as gray values (default)" );
    
    cl.AddSwitch( Key( 'L', "log" ), &LogData, true, "Apply log() function to data" );
    cl.AddSwitch( Key( 'e', "exp" ), &ExpData, true, "Apply exp() function to data" );
    
    cl.AddOption( Key( 'n', "num-bins" ), &NumberOfHistogramBins, "Number of histogram bins." );
    
    cl.AddSwitch( Key( 'C', "column" ), &WriteAsColumn, true, "Write statistics in column format rather than line format." );
    cl.AddSwitch( Key( 'E', "e-notation" ), &OutputExpNotation, true, "Write floating point numbers in #e# notation (i.e., 1e-3 instead of 0.001)" ); 
    cl.AddOption( Key( 'c', "com" ), &LandmarksFileName, "Write centers of mass of labeled regions to landmarks file." );
    
    cl.AddOption( Key( 'm', "mask" ), &MaskFileName, "Analyze region based on binary mask file", &MaskIsBinary );
    cl.AddOption( Key( 'M', "Mask" ), &MaskFileName, "Analyze regions separately based on mask file (label field)." );
    
    cl.AddCallback( Key( 'p', "percentile" ), CallbackAddPercentile, "Add value to list of percentile to compute." );
    
    cl.Parse();
        
    const char* next = cl.GetNext();
    while ( next )
      {
      ImageFileNames.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e;
    exit( 1 );
    }
  
  cmtk::TypedArray::SmartPtr maskData( NULL );
  if ( MaskFileName ) 
    {
    cmtk::UniformVolume::SmartPtr mask( cmtk::VolumeIO::ReadOriented( MaskFileName, Verbose ) );
    if ( ! mask ) 
      {
      cmtk::StdErr << "ERROR: could not read mask file " << MaskFileName << "\n";
      exit( 1 );
      }
    maskData = mask->GetData();
    if ( ! maskData ) 
      {
      cmtk::StdErr << "ERROR: could not read data from mask file " << MaskFileName << "\n";
      exit( 1 );
      }
    
    if ( MaskIsBinary )
      {
      maskData->Binarize();
      }
    }
  
  std::list<const char*>::const_iterator it = ImageFileNames.begin();
  for ( ; it != ImageFileNames.end(); ++it )
    {
    const char* imageFileName = *it;
    cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( imageFileName, Verbose ) );
    if ( ! volume ) 
      {
      cmtk::StdErr << "ERROR: could not read image file " << imageFileName << "\n";
      continue;
      }
    
    cmtk::TypedArray::SmartPtr volumeData = volume->GetData();
    if ( ! volumeData ) 
      {
      cmtk::StdErr << "ERROR: could not read data from image file " << imageFileName << "\n";
      continue;
      }
    
    fprintf( stdout, "Statistics for image %s\n", imageFileName );
    
    if ( Label ) 
      {
      AnalyzeLabels( volume, maskData );
      } 
    else 
      {
      if ( LogData )
	volumeData->ApplyFunctionDouble( cmtk::Wrappers::Log );
      if ( ExpData )
	volumeData->ApplyFunctionDouble( cmtk::Wrappers::Exp );
      if ( maskData )
	AnalyzeGrey( volume, maskData );
      else
	AnalyzeGrey( volume );
      }
    }
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace statistics
} // namespace apps
} // namespace cmtk
#endif

