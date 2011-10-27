/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

/** \addtogroup Base */
//@{

/** Vector addition operator.
 * Two vectors are added elementwise. A newly created vector object is 
 * returned. This operator is rather inefficient in terms of allocation and
 * destruction of objects; use -= instead if possible.
 */
template<class T>
Vector<T> operator+ ( const Vector<T>& p, const Vector<T>& delta )
{
  assert ( p.Dim == delta.Dim );

  T* Result = Memory::ArrayC::Allocate<T>( p.Dim );
#pragma omp parallel for if (p.Dim>1e4)
  for ( size_t i=0; i<p.Dim; ++i )
    Result[i] = p.Elements[i] + delta.Elements[i];

  return Vector<T>( p.Dim, Result );
}

/** Vector subtraction operator.
 * Two vectors are subtracted elementwise. A newly created vector object is 
 * returned. This operator is rather inefficient in terms of allocation and
 * destruction of objects; use += instead if possible.
 */
template<class T>
inline Vector<T> operator- 
( const Vector<T>& p, const Vector<T>& delta )
{
  assert ( p.Dim == delta.Dim );

  T* Result = Memory::ArrayC::Allocate<T>( p.Dim );
#pragma omp parallel for if (p.Dim>1e4)
  for ( size_t i=0; i<p.Dim; ++i )
    Result[i] = p.Elements[i] - delta.Elements[i];

  return Vector<T>( p.Dim, Result );
}

/** Scalar-to-vector multiplication operator.
 * Every element of a vector is multiplies by the same scalar factor. The 
 * result is returned as an automatically created object.  This operator is 
 * rather inefficient in terms of allocation and destruction of objects; use
 * *= instead if possible.
 */
template<class T>
Vector<T> operator* ( const T c, const Vector<T>& p ) 
{
  T* Result = Memory::ArrayC::Allocate<T>( p.Dim );
#pragma omp parallel for if (p.Dim>1e4)
  for ( size_t i=0; i<p.Dim; ++i )
    Result[i] = c * p.Elements[i];

  return Vector<T>( p.Dim, Result );
}

/** Coordinatewise multiplication operator.
 * Two vectors are multiplied element by element. The result is returned as an
 * automatic variable.
 */
template<class T>
Vector<T> CoordMult ( const Vector<T>& p, const Vector<T>& q ) 
{
  assert ( p.Dim == q.Dim );

  T* Result = Memory::ArrayC::Allocate<T>( p.Dim );
#pragma omp parallel for if (p.Dim>1e4)
  for ( size_t i=0; i<p.Dim; ++i )
    Result[i] = p.Elements[i] * q.Elements[i];

  return Vector<T>( p.Dim, Result );
}

/** Scalar product.
 * This operator computes the standard scalar product of two vectors over the
 * same primitive type. As only a primitive object is returned as the result of
 * this operator, it is time- and memory-efficient.
 */
template<class T> 
inline T operator* ( const Vector<T>& p, const Vector<T>& q ) 
{
  assert ( p.Dim == q.Dim );

  T Result = 0;
#pragma omp parallel for if (p.Dim>1e4)
  for ( int i=0; i<static_cast<int>( p.Dim ); ++i )
    Result += p.Elements[i] * q.Elements[i];

  return Result;
}

} // namespace cmtk
