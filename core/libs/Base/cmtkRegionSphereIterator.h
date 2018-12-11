/*
//
//  Copyright 2010-2012 SRI International
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

#ifndef __cmtkRegionSphereIterator_h_included_
#define __cmtkRegionSphereIterator_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkRegion.h>
#include <Base/cmtkFixedVector.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

#include <list>

namespace
cmtk
{
/// Iterator for spherical region of n-dimensional image.
template<typename TRegion>
class RegionSphereIterator
{
public:
  /// This class.
  typedef RegionSphereIterator<TRegion> Self;

  /// Region type.
  typedef TRegion RegionType;

  /// Region dimension.
  static const size_t Dimension = RegionType::Dimension;

  /// Region scalar type.
  typedef typename RegionType::ScalarType ScalarType;

  /// Index type.
  typedef FixedVector<Self::Dimension,typename Self::ScalarType> IndexType;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Internal index list type.
  typedef std::list<typename Self::IndexType> IndexListType;

protected:
  /// Protected default constructor - only to be used by derived class.
  RegionSphereIterator() {}
  
public:
  /// Constructor with radius only (center is zero-index).
  explicit RegionSphereIterator( const typename Self::IndexType radius /*!< Radius of the sphere in each dimension.*/ )
  {
    typename Self::IndexType center( typename Self::IndexType::Init( 0 ) );
    this->Populate( radius, center, 0, 1.0 );
  }

  /// Constructor with radius and (not necessarily) non-zero center
  explicit RegionSphereIterator( const typename Self::IndexType radius /*!< Radius of the sphere in each dimension.*/, const typename Self::IndexType center /*!< Center index of the sphere. */ )
  {
    this->Populate( radius, center, 0, 1.0 );
  }

  /// Increment operator.
  typename Self::IndexListType::const_iterator& operator++()
  {
    return ++this->m_IndexListIterator;
  }
  
  /// Get index.
  const typename Self::IndexType& Index() const
  {
    return *(this->m_IndexListIterator);
  }

  /// Iterator assignment.
  const Self& operator=( const typename Self::IndexListType::const_iterator& it )
  {
    this->m_IndexListIterator = it;
    return *this;
  }
  
  /// Region "begin" index.
  const typename Self::IndexListType::const_iterator begin() const
  {
    return this->m_IndexList.begin();
  }
  
  /// Region "end" index.
  const typename Self::IndexListType::const_iterator end() const
  {
    return this->m_IndexList.end();
  }

  /// Equality operator.
  bool operator==( const typename Self::IndexListType::const_iterator& it )
  {
    return it == this->m_IndexListIterator;
  }

  /// Inequality operator.
  bool operator!=( const typename Self::IndexListType::const_iterator& it )
  {
    return it != this->m_IndexListIterator;
  }
  
protected:
  /// Pre-computed list of grid indexes on the sphere.
  typename Self::IndexListType m_IndexList;

  /// Current position in index list.
  typename Self::IndexListType::const_iterator m_IndexListIterator;

  /// Recursively populate the list of indexes.
  virtual void Populate( const typename Self::IndexType& radius /*!< Sphere radius in index steps by dimension.*/, const typename Self::IndexType& center /*!< Sphere center. */, const size_t dim /*!< Next dimension. */, 
			 const double remainSquare /*!< Remaining proportion of total squared sphere radius. */)
  {
    if ( remainSquare >= 0 )
      {
      typename Self::IndexType index = center;
      const int radiusThisDimension = static_cast<int>( sqrt( remainSquare ) * radius[dim] );
      if ( dim < Self::Dimension )
	{
	this->Populate( radius, index, dim+1, remainSquare );

	for ( int r = 1; r <= radiusThisDimension; ++r )
	  {
	  const double newRemainSquare = remainSquare - MathUtil::Square(1.0 * r / radius[dim] );

	  index[dim] = center[dim]+r;
	  this->Populate( radius, index, dim+1, newRemainSquare );
	  
	  index[dim] = center[dim]-r;
	  this->Populate( radius, index, dim+1, newRemainSquare );
	  }
	}
      else
	{
	this->m_IndexList.push_back( center );
	}
      }
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkRegionSphereIterator_h_included_
