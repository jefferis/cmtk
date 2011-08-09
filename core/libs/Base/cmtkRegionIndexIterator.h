/*
//
//  Copyright 2010-2011 SRI International
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
template<typename TRegion>
class RegionIndexIterator
{
public:
  /// This class.
  typedef RegionIndexIterator<TRegion> Self;

  /// Region type.
  typedef TRegion RegionType;

  /// Region dimension.
  static const size_t Dimension = RegionType::Dimension;

  /// Region scalar type.
  typedef typename RegionType::ScalarType ScalarType;

  /// Index type.
  typedef FixedVector<Self::Dimension,Self::ScalarType> IndexType;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Constructor from two index, from and to.
  RegionIndexIterator( const typename Self::RegionType& region )
    : m_Region( region ),
      m_Index( region.From() )
  {
    // "End" index is one after last valid element.
    this->m_End = this->m_Region.From();
    this->m_End[Self::Dimension-1] = this->m_Region.To()[Self::Dimension-1];    
  }

  /// Increment operator.
  Self& operator++()
  {
    for ( size_t idx = 0; idx < Self::Dimension; ++idx)
      {
      if ( (++this->m_Index[idx]) >= this->m_Region.To()[idx] )
	{
	if ( idx+1 < Self::Dimension )
	  this->m_Index[idx] = this->m_Region.From()[idx];
	}
      else
	break;
      }
    return *this;
  }
  
  /// Get index.
  const typename Self::IndexType& Index() const
  {
    return this->m_Index;
  }

  /// Region "begin" index.
  const typename Self::IndexType begin() const
  {
    return this->m_Region.From();
  }
  
  /// Region "end" index.
  const typename Self::IndexType& end() const
  {
    return this->m_End;
  }

  /// Assign index.
  Self& operator=( const typename Self::IndexType& index )
  {
    this->m_Index = index;
    return *this;
  }

  /// Index equality.
  bool operator==( const typename Self::IndexType& index )
  {
    return (this->m_Index == index);
  }

  /// Index inequality.
  bool operator!=( const typename Self::IndexType& index )
  {
    return !(this->m_Index == index);
  }

private:
  /// Iterated region.
  typename Self::RegionType m_Region;

  /// End index (i.e., first non-valid index).
  typename Self::IndexType m_End;

  /// Current index.
  typename Self::IndexType m_Index;
};

} // namespace cmtk

#endif // #ifndef __cmtkRegionIndexIterator_h_included_
