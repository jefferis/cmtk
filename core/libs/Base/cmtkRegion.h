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
#include <fstream>

namespace
cmtk
{
/// Class for n-dimensional image index.
template<size_t NDIM,typename T=int>
class Region
{
public:
  /// This class.
  typedef Region<NDIM,T> Self;

  /// Index type.
  typedef Index<NDIM,T> IndexType;

  /// Default constructor.
  Region() {}

  /// Constructor from two index, from and to.
  Region( const IndexType& fromIndex, const IndexType& toIndex )
  {
    this->m_RegionFrom = fromIndex;
    this->m_RegionTo = toIndex;
  }
  
  /// Get "from".
  IndexType& From()
  {
    return this->m_RegionFrom;
  }

  /// Get const "from".
  const IndexType& From() const
  {
    return this->m_RegionFrom;
  }

  /// Get "from".
  IndexType& To()
  {
    return this->m_RegionTo;
  }

  /// Get const "from".
  const IndexType& To() const
  {
    return this->m_RegionTo;
  }

  /// Compute region size (e.g., number of pixels for grid regions).
  T Size() const
  {
    T size = 0;
    for ( size_t i = 0; i < NDIM; ++i )
      size *= (this->m_RegionTo[i]-this->m_RegionFrom[i]);
    return size;
  }
  
private:
  /// Beginning index.
  IndexType m_RegionFrom;

  /// End index.
  IndexType m_RegionTo;
};

/// Stream input operator.
template<size_t NDIM,typename T>
std::ofstream& operator<<( std::ofstream& stream, const Region<NDIM,T>& region )
{
  return stream << region.From() << region.To();
}

/// Stream output operator.
template<size_t NDIM,typename T>
std::ofstream& operator>>( std::ofstream& stream, Region<NDIM,T>& region )
{
  return stream >> region.From() >> region.To();
}

} // namespace cmtk

#endif // #ifndef __cmtkRegion_h_included_
