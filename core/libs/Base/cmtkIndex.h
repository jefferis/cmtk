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

#ifndef __cmtkIndex_h_included_
#define __cmtkIndex_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{
/// Class for n-dimensional image index.
template<size_t NDIM>
class Index
{
public:
  /// This class.
  typedef Index<NDIM> Self;

  /// Default constructor.
  Index() {}

  /// Constructor from const int array.
  Index( const int (&indexArray)[NDIM] )
  {
    for ( size_t i = 0; i < NDIM; ++i )
      this->m_Index[i] = indexArray[i];
  }

  /// Access operator.
  int& operator[]( const size_t idx )
  {
    return this->m_Index[idx];
  }

  /// Constant access operator.
  const int& operator[]( const size_t idx ) const
  {
    return this->m_Index[idx];
  }
  
private:
  /// The actual index array.
  int m_Index[NDIM];
};

} // namespace cmtk

#endif // #ifndef __cmtkIndex_h_included_
