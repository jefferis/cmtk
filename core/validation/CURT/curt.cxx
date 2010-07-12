/*
//
//  Copyright 2010 SRI International
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

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>

#include <cmtkTypedArrayFunctionHistogramMatching.h>

#include <vector>
#include <algorithm>

class
IndexValue
{
public:
  /// Constructor.
  IndexValue( const size_t index = 0, const double value = 0 ) : m_Index( index ), m_Value( value ) {};

  /// Pixel index.
  size_t m_Index;

  /// Pixel value.
  double m_Value;
};

bool
operator<( const IndexValue& x, const IndexValue& y )
{
  return x.m_Value < y.m_Value;
}

int
main( const int argc, const char* argv[] )
{
  bool verbose = false;

  const char* pathFix = NULL;
  const char* pathMov = NULL;

  std::vector<std::string> pathsLbls;

  bool padZero = false;
  bool collapse = false;

  try
    {
    cmtk::CommandLine cl( argc, argv, cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Completely Useless Registration Tool" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "For experimental demonstration ONLY! This program coregisters two intensity images by simple permutation and reformats one or more matching label maps accordingly." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Validation" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Be verbose" )->SetProperties( cmtk::CommandLine::PROPS_NOXML );

    cl.BeginGroup( "Preprocessing", "Input Image Preprocessing" );
    cl.AddSwitch( Key( "pad" ), &padZero, true, "Pad (ignore) zero-filled areas in both images" );
    cl.AddSwitch( Key( "collapse" ), &collapse, true, "Collapse transformation for equal reference image values to the same floating image pixel." );
    cl.EndGroup();

    cl.AddParameter( &pathFix, "FixedImage", "Fixed image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &pathMov, "MovingImage", "Moving image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameterVector( &pathsLbls, "LabelImages", "Label image paths" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    
    cl.Parse();
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  cmtk::UniformVolume::SmartPtr refImage( cmtk::VolumeIO::ReadOriented( pathFix, verbose ) );
  cmtk::UniformVolume::SmartPtr fltImage( cmtk::VolumeIO::ReadOriented( pathMov, verbose ) );

  if ( padZero )
    {
    refImage->GetData()->SetPaddingValue( 0 );
    fltImage->GetData()->SetPaddingValue( 0 );
    }

  const size_t nPixelsRef = refImage->GetNumberOfPixels();
  std::vector<IndexValue> refIndexValue( nPixelsRef );
  for ( size_t i = 0; i < nPixelsRef; ++i )
    refIndexValue[i] = IndexValue( i, refImage->GetDataAt( i ) );
  std::sort( refIndexValue.begin(), refIndexValue.end() );

  const size_t nPixelsFlt = fltImage->GetNumberOfPixels();
  std::vector<IndexValue> fltIndexValue( nPixelsFlt );
  for ( size_t i = 0; i < nPixelsFlt; ++i )
    fltIndexValue[i] = IndexValue( i, fltImage->GetDataAt( i ) );
  std::sort( fltIndexValue.begin(), fltIndexValue.end() );

  std::vector<size_t> refIndexLookup( nPixelsRef );
  for ( size_t i = 0; i < nPixelsRef; ++i )
    {
    refIndexLookup[refIndexValue[i].m_Index] = i;
    }

  const double factor = double( nPixelsFlt ) / double( nPixelsRef );
  std::vector<size_t> lookup( nPixelsRef );

  if ( collapse )
    {
    for ( size_t iRef = 0; iRef < nPixelsRef; ++iRef )
      {
      const size_t rangeStart = iRef;
      for ( const double rangeValue = refIndexValue[iRef].m_Value; (iRef < nPixelsRef) && (refIndexValue[iRef].m_Value == rangeValue); ++iRef );
      
      const size_t translation = fltIndexValue[static_cast<size_t>( 0.5 * factor * (rangeStart + iRef))].m_Index;
      for ( size_t ii = rangeStart; ii < iRef; ++ii )
	{
	lookup[refIndexValue[ii].m_Index] = translation;
	}
      }
    }
  else
    {
    for ( size_t iRef = 0; iRef < nPixelsRef; ++iRef )
      {
      lookup[iRef] = fltIndexValue[static_cast<size_t>(refIndexLookup[iRef]*factor)].m_Index;
      }
    }
  
  for ( size_t iRef = 0; iRef < nPixelsRef; ++iRef )
    {
    refImage->SetDataAt( fltImage->GetDataAt( lookup[iRef] ), iRef );
    }
  
  cmtk::VolumeIO::Write( *refImage, "reformat.nii", verbose );

  for ( size_t l = 0; l < pathsLbls.size(); ++l )
    {
    cmtk::UniformVolume::SmartPtr lblImage( cmtk::VolumeIO::ReadOriented( pathsLbls[l].c_str(), verbose ) );
    
    for ( size_t iRef = 0; iRef < nPixelsRef; ++iRef )
      {
      refImage->SetDataAt( lblImage->GetDataAt( lookup[iRef] ), iRef );
      }
    
    char output[PATH_MAX];
    snprintf( output, PATH_MAX, "labels%d.nii", static_cast<int>( l ) );
    cmtk::VolumeIO::Write( *refImage, output, verbose );    
    }
}
