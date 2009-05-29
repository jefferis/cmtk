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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkRectangle_h_included_
#define __cmtkRectangle_h_included_

#include <cmtkconfig.h>

#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// N-dimensional rectangles.
template<class C,int N>
class Rectangle 
{
public:
  /// Beginning of rectangle.
  C From[N];

  /// End of rectangle.
  C To[N];

  /// Default constructor: do nothing.
  Rectangle() {}

  /// Explicit init constructor.
  Rectangle( const C* from, const C* to ) 
  {
    this->Set( from, to );
  }

  /// Copy constructor.
  Rectangle( const Rectangle<C,N>& src ) 
  {
    this->Set( src.From, src.To );
  }

  /// Set beginning and end of rectangle.
  Rectangle<C,N>& Set( const C* from, const C* to ) 
  {
    memcpy( From, from, sizeof( From ) );
    memcpy( To, to, sizeof( To ) );
    return *this;
  }

  /// Set beginning of rectangle.
  Rectangle<C,N>& SetFrom( const C* from ) 
  {
    memcpy( From, from, sizeof( From ) );
    return *this;
  }

  /// Set end of rectangle.
  Rectangle<C,N>& SetTo( const C* to ) 
  {
    memcpy( To, to, sizeof( To ) );
    return *this;
  }

  /// Validate rectangle and switch upper and lower boundaries if necessary.
  Rectangle<C,N>& Validate() 
  {
    for ( int n = 0; n < N; ++n )
      if ( From[n] > To[n] ) 
	MathUtil::Swap( From[n], To[n] );
    return *this;
  }

  /// Crop rectangle.
  Rectangle<C,N>& Crop( const C* from, const C* to ) 
  {
    for ( int n = 0; n < N; ++n ) 
      {
      From[n] = ( From[n] < from[n] ) ? from[n] : From[n];
      To[n] = ( To[n] > to[n] ) ? to[n] : To[n];
      }
    return *this;
  }

  /// Expand rectangle by factor about its center.
  Rectangle<C,N>& Expand( const C factor ) 
  {
    for ( int n = 0; n < N; ++n ) 
      {
      const C length = (To[n] - From[n]) * factor;
      const C center = (To[n] + From[n]) / 2;
      
      From[n] = center - length / 2;
      To[n] = center + length / 2;
      }
    return *this;
  }

  /// Expand rectangle by margin on all sides.
  Rectangle<C,N>& ExpandMargin( const int margin ) 
  {
    for ( int n = 0; n < N; ++n ) 
      {
      From[n] -= margin;
      To[n] += margin;
      }
    return *this;
  }
  
  /// Get center of region.
  template<class T> void GetCenter( T center[N] ) const 
  {
    for ( int n = 0; n < N; ++n )
      center[n] = static_cast<T>( .5 * static_cast<T>(From[n] + To[n]) );
  }

  /// Get center of integer region.
  template<class T> void GetCenterInt( T center[N] ) const 
  {
    for ( int n = 0; n < N; ++n )
      center[n] = static_cast<T>( .5 * static_cast<T>(From[n] + To[n] - 1) );
  }

  /// Recturn volume of rectangle.
  C GetVolume() const 
  {
    C volume = 1;
    for ( int n = 0; n < N; ++n ) volume *= (To[n]-From[n]);
    return volume;
  }

  /// Assignment operator.
  Rectangle<C,N>& operator=( const Rectangle<C,N>& src ) 
  {
    this->Set( src.From, src.To );
    return *this;
  }
};

typedef Rectangle<int,2> IntROI2D;
typedef Rectangle<float,2> FloatROI2D;
typedef Rectangle<Types::Coordinate,2> CoordinateROI2D;

typedef Rectangle<int,3> IntROI3D;
typedef Rectangle<Types::Coordinate,3> CoordinateROI3D;

//@}

} // namespace cmtk

#endif
