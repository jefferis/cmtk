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

#ifndef __cmtkRegionSphereSurfaceIterator_h_included_
#define __cmtkRegionSphereSurfaceIterator_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkRegionSphereIterator.h>

namespace
cmtk
{
/// Iterator for surface of a spherical region of n-dimensional image.
template<typename TRegion>
class RegionSphereSurfaceIterator : public RegionSphereIterator<TRegion>
{
public:
  /// This class.
  typedef RegionSphereSurfaceIterator<TRegion> Self;

  /// Parent class.
  typedef RegionSphereIterator<TRegion> Superclass;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Re-declare index list type.
  typedef typename Superclass::IndexListType IndexListType;

  /// Constructor with radius only (center is zero-index).
  explicit RegionSphereSurfaceIterator( const typename Self::IndexType radius /*!< Radius of the sphere in each dimension.*/ )
  {
    typename Self::IndexType center( typename Self::IndexType::Init( 0 ) );
    this->Populate( radius, center, 0, 1.0 );
  }

  /// Constructor with radius and (not necessarily) non-zero center
  explicit RegionSphereSurfaceIterator( const typename Self::IndexType radius /*!< Radius of the sphere in each dimension.*/, const typename Self::IndexType center /*!< Center index of the sphere. */ )
  {
    this->Populate( radius, center, 0, 1.0 );
  }

  /// Iterator assignment.
  const Self& operator=( const typename Self::IndexListType::const_iterator& it )
  {
    this->m_IndexListIterator = it;
    return *this;
  }
  
protected:
  /// Recursively populate the list of indexes.
  virtual void Populate( const typename Self::IndexType& radius /*!< Sphere radius in index steps by dimension.*/, const typename Self::IndexType& center /*!< Sphere center. */, const size_t dim /*!< Next dimension. */, 
			 const double remainSquare /*!< Remaining proportion of total squared sphere radius. */)
  {
    if ( remainSquare > 0 )
      {
      typename Self::IndexType index = center;
      const int radiusThisDimension = static_cast<int>( sqrt( remainSquare ) * radius[dim] );
      if ( dim < Self::Dimension-1 )
	{
	this->Populate( radius, center, dim+1, remainSquare );
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
	for ( int r = -radiusThisDimension; r <= radiusThisDimension; r += 2*radiusThisDimension )
	  {
	  index[dim] = center[dim]+r;
	  this->m_IndexList.push_back( index );
	  }
	}
      }
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkRegionSphereSurfaceIterator_h_included_
