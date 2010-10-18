/*
//
//  Copyright 2010 SRI International
//
//  Copyright 2010 Torsten Rohlfing
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

#ifndef __cmtkRegionIndexIterator_h_included_
#define __cmtkRegionIndexIterator_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkRegion.h>
#include <Base/cmtkFixedVector.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

namespace
cmtk
{
/// Class for n-dimensional image index.
template<size_t NDIM,typename T=int>
class RegionIndexIterator
{
public:
  /// This class.
  typedef RegionIndexIterator<NDIM,T> Self;

  /// Region type.
  typedef Region<NDIM,T> RegionType;

  /// Index type.
  typedef FixedVector<NDIM,T> IndexType;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Constructor from two index, from and to.
  RegionIndexIterator( const typename Self::RegionType& region )
    : m_Region( region )
  {}

  /// Increment operator.
  Self& operator++()
  {
    int idx = NDIM-1;
    while ( idx >= 0 )
      {
      if ( (++this->m_Index[idx]) >= this->m_Region.To()[idx] )
	{
	this->m_Index[idx] = this->m_Region.From()[idx];
	--idx;
	}
      }
  }

  /// Get index.
  const typename Self::IndexType& Index() const
  {
    return this->m_Index;
  }

private:
  /// Beginning index.
  typename Self::RegionType m_Region;

  /// End index.
  typename Self::IndexType m_Index;
};

} // namespace cmtk

#endif // #ifndef __cmtkRegionIndexIterator_h_included_
