/*
//
//  Copyright 2010-2011 SRI International
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
#include <System/cmtkConsole.h>
#include <System/cmtkFileUtils.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkValueSequence.h>

#include <IO/cmtkVolumeIO.h>

#include <vector>
#include <algorithm>

#include <QImage>
#include <QImageWriter>
#include <QColor>
#include <QRgb>

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
  const char* pathFix = NULL;
  const char* pathMov = NULL;
  const char* pathOut = NULL;
  const char* pathDif = NULL;

  try
    {
    cmtk::CommandLine cl( cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Completely Useless Registration Tool" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "For experimental demonstration ONLY! This program coregisters two intensity images by simple permutation and reformats one or more matching label maps accordingly." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Validation" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddParameter( &pathFix, "FixedImage", "Fixed image path" );
    cl.AddParameter( &pathMov, "MovingImage", "Moving image path" );
    cl.AddParameter( &pathOut, "OutputImage", "Output image path" );
    cl.AddParameter( &pathDif, "DiffImage", "Difference image path" );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  // Read input images
  QImage fixImage( pathFix );
  QImage movImage( pathMov );
  QImage outImage( movImage );
  QImage difImage( movImage );

  // Get sorted list of all fixed image pixels
  const size_t nPixelsFix = fixImage.width() * fixImage.height();
  std::vector<IndexValue> fixIndexValue( nPixelsFix );

  for ( size_t i = 0; i < nPixelsFix; ++i )
    {
    const float value = QColor( fixImage.pixel( i / fixImage.width(), i % fixImage.width() ) ).red();
    fixIndexValue[i] = IndexValue( i, value );
    }
  std::sort( fixIndexValue.begin(), fixIndexValue.end() );
  
  // Get sorted list of all moving pixels
  const size_t nPixelsMov = movImage.width() * movImage.height();
  std::vector<IndexValue> movIndexValue( nPixelsMov );

  for ( size_t i = 0; i < nPixelsMov; ++i )
    {
    const float value = QColor( movImage.pixel( i / movImage.width(), i % movImage.width() ) ).red();
    movIndexValue[i] = IndexValue( i, value );
    }
  std::sort( movIndexValue.begin(), movIndexValue.end() );
  
  const double factor = double( nPixelsMov ) / double( nPixelsFix );
  
  for ( size_t i = 0; i < nPixelsFix; ++i )
    {
    const size_t movIndex = movIndexValue[static_cast<size_t>(0.5+i*factor)].m_Index;
    const size_t fixIndex = fixIndexValue[i].m_Index;
    
    const int px = QColor( movImage.pixel( movIndex / movImage.width(), movIndex % movImage.width() ) ).red();
    outImage.setPixel( fixIndex / fixImage.width(), fixIndex % fixImage.width(), px );

    const int px2 = QColor( fixImage.pixel( fixIndex / fixImage.width(), fixIndex % fixImage.width() ) ).red();
    difImage.setPixel( fixIndex / fixImage.width(), fixIndex % fixImage.width(), (px-px2)/2+128 );
    }
  
  outImage.save( pathOut );
  difImage.save( pathDif );
}
