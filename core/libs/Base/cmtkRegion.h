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

#ifndef __cmtkRegion_h_included_
#define __cmtkRegion_h_included_

#include <cmtkconfig.h>

#include <cmtkIndex.h>

namespace
cmtk
{
/// Class for n-dimensional image index.
template<size_t NDIM>
class Region
{
public:
  /// This class.
  typedef Region<NDIM> Self;

  /// Index type.
  typedef Index<NDIM> IndexType;

  /// Default constructor.
  Region() {}

  /// Constructor from two index, from and to.
  Region( const IndexType& fromIndex, const IndexType& toIndex )
  {
    this->m_RegionFrom = fromIndex;
    this->m_RegionTo = toIndex;
  }
  
  /// Copy constructor.
  Region( const Self& other )
  {
    *this = other;
  }

  /// Assignment operator.
  Self& operator=( const Self& other )
  {
    this->m_RegionFrom = other.m_RegionFrom;
    this->m_RegionTo = other.m_RegionTo;
    return *this;
  }

private:
  /// Beginning index.
  IndexType m_RegionFrom;

  /// End index.
  IndexType m_RegionTo;
};

} // namespace cmtk

#endif // #ifndef __cmtkRegion_h_included_
