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

namespace
cmtk
{

/// Stream input operator.
template<size_t NDIM,typename T>
std::ofstream& operator<<( std::ofstream& stream, const Region<NDIM,T>& region )
{
  return stream << region.From() << region.To();
}

/// Stream output operator.
template<size_t NDIM,typename T>
std::ifstream& operator>>( std::ifstream& stream, Region<NDIM,T>& region )
{
  return stream >> region.From() >> region.To();
}

template<size_t NDIM,typename T>
const Region<NDIM,T> 
Region<NDIM,T>::GetFaceRegion( const int dim, const bool upper ) const
{
  Self::IndexType from = this->m_RegionFrom;
  Self::IndexType to = this->m_RegionTo;

  if ( upper )
    {
    from[dim] = to[dim]-1;
    }
  else
    {
    to[dim] = from[dim]+1;
    }

  return Self( from, to );
}

} // namespace cmtk
