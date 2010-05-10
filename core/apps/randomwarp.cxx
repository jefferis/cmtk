/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include <cmtkVolumeIO.h>
#include <cmtkSplineWarpXform.h>
#include <cmtkClassStream.h>
#include <cmtkMathUtil.h>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_TIME_H
#  include <time.h>
#endif

int
main( const int argc, const char *argv[] )
{
  if ( argc < 4 )
    {
    cmtk::StdErr << "randomwarp imageFile gridSpacing sigma outlist0 [outlist1...]\n";
    exit( 2 );
    }

  const char *filenameIn = argv[1];
  const cmtk::Types::Coordinate gridSpacing = atof( argv[2] );
  const cmtk::Types::Coordinate sigma = atof( argv[3] );

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( filenameIn, true ) );
  if ( ! volume ) return 1;

  cmtk::SplineWarpXform::SmartPtr warp( new cmtk::SplineWarpXform( volume, gridSpacing ) );

  // seed RNG with (supposedly) random system time
  cmtk::MathUtil::NormalRandom( sigma, static_cast<unsigned int>( time( NULL ) ) );

  for ( int outFileIdx = 4; outFileIdx < argc; ++outFileIdx ) 
    {
    const char *filenameOut = argv[outFileIdx];
    
    warp->InitControlPoints();
    cmtk::Types::Coordinate* coeff = warp->m_Parameters;
    for ( unsigned int idx=0; idx < warp->m_NumberOfParameters; ++idx, ++coeff ) 
      {
      *coeff += cmtk::MathUtil::NormalRandom( sigma );
      }
    
    cmtk::ClassStream outStream( filenameOut, "studylist", 
			      cmtk::ClassStream::WRITE );
    if ( outStream.IsValid() )
      {
      outStream.Begin( "studylist" );
      outStream.WriteInt( "num_sources", 2 );
      outStream.End();
      
      outStream.Begin( "source" );
      outStream.WriteString( "studyname", filenameIn );
      outStream.End();
      
      outStream.Begin( "source" );
      outStream.WriteString( "studyname", filenameIn );
      outStream.End();
      
      outStream.Close();
      }
    
    outStream.Open( filenameOut, "registration", cmtk::ClassStream::WRITE );
    if ( outStream.IsValid() ) 
      {
      outStream.Begin( "registration" );
      outStream.WriteString( "reference_study", filenameIn );
      outStream.WriteString( "model_study", filenameIn );
      outStream << warp;
      outStream.Close();
      }
    }
  
  return 0;
}

