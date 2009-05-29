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

#ifndef __cmtkScaleHistogramValueTrait_h_included_
#define __cmtkScaleHistogramValueTrait_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Histogram kernel scaling traits.
template<class T> 
class ScaleHistogramValueTrait
{
public:
  /// In general return original value as scaled value.
  static T Scale( const float value ) 
  {
    return value;
  }
};

template<> 
class ScaleHistogramValueTrait<int>
{
public:
  /// For integer bins, scale with factor 256 for added resolution.
  static int Scale( const float value )
  {
    return static_cast<int>( 256 * value );
  }
};

template<> 
class ScaleHistogramValueTrait<unsigned int>
{
public:
  /// For integer bins, scale with factor 256 for added resolution.
  static unsigned int Scale( const float value )
  {
    return static_cast<unsigned int>( 256 * value );
  }
};

//@}

} // namespace cmtk

#endif
