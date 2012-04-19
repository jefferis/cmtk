/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include <Base/cmtkDataTypeTraits.h>

namespace
cmtk
{

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>::FixedSquareMatrix()
{
  memset( this->m_Matrix, 0, sizeof( this->m_Matrix ) );
  for ( size_t i = 0; i < Self::Dimension; ++i )
    this->m_Matrix[i][i] = 1.0;
}

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>::FixedSquareMatrix( const Self& other )
{
  memcpy( this->m_Matrix, other.m_Matrix, sizeof( this->m_Matrix ) );
}

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>&
FixedSquareMatrix<NDIM,TSCALAR>::Set( const typename Self::ScalarType *const values )
{
  memcpy( this->m_Matrix, values, sizeof( this->m_Matrix ) );
  return *this;
}

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>&
FixedSquareMatrix<NDIM,TSCALAR>::Fill( const typename Self::ScalarType& value )
{
  for ( size_t i = 0; i < Self::Dimension; ++i )
    for ( size_t j = 0; j < Self::Dimension; ++j )
      this->m_Matrix[i][j] = value;
  return *this;
}

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>::GetTranspose() const
{
  Self transpose;
  for ( size_t i = 0; i < Self::Dimension; ++i ) 
    {
    for ( size_t j = 0; j < Self::Dimension; ++j )
      transpose[i][j] = this->m_Matrix[j][i];
}
  return transpose;
}
  
template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>::GetInverse() const
{
  Self inverse;

  typename Self::ScalarType matrix[NDIM][NDIM];
  memcpy( matrix, this->m_Matrix, sizeof( matrix ) );
  
  typename Self::ScalarType rowBuff[Self::Dimension];
  for ( size_t col = 0; col< Self::Dimension; ++col ) 
    {
    size_t pivIdx = col;
    typename Self::ScalarType pivAbs = fabs( matrix[col][col] );
    
    for ( size_t row = col+1; row< Self::Dimension; ++row ) 
      {
      const typename Self::ScalarType nextAbs = fabs( matrix[row][col] );
      if (nextAbs > pivAbs ) 
	{
	pivIdx = row;
	pivAbs = nextAbs;
	}
      }
    
    if ( col != pivIdx )
      {
      memcpy( rowBuff, matrix[col], sizeof(rowBuff) );
      memcpy( matrix[col], matrix[pivIdx], sizeof(rowBuff) );
      memcpy( matrix[pivIdx], rowBuff, sizeof(rowBuff) );

      memcpy( rowBuff, inverse[col], sizeof(rowBuff) );
      memcpy( inverse[col], inverse[pivIdx], sizeof(rowBuff) );
      memcpy( inverse[pivIdx], rowBuff, sizeof(rowBuff) );
      }
    
    for ( size_t c=0; c< Self::Dimension; ++c ) 
      {
      if (c>col )
	matrix[col][c] /= matrix[col][col];
      inverse[col][c] /= matrix[col][col];
      }
    matrix[col][col] = DataTypeTraits<TSCALAR>::One();
    
    for ( size_t row = 0; row< Self::Dimension; ++row ) 
      {
      if (row != col ) 
	{
	for ( size_t c=0; c< Self::Dimension; ++c ) 
	  {
	  if ( c>col ) 
	    matrix[row][c] -= matrix[row][col] * matrix[col][c];
	  inverse[row][c] -= matrix[row][col] * inverse[col][c];
	  }
	matrix[row][col] = 0;
	}
      }
    }

  return inverse;
}

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>& 
FixedSquareMatrix<NDIM,TSCALAR>::operator*=( const Self& other )
{
  return (*this = ((*this) * other));
}

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>& 
FixedSquareMatrix<NDIM,TSCALAR>::operator*=( const typename Self::ScalarType& scalar )
{
  for ( size_t j=0; j < Self::Dimension; ++j ) 
    {
    for ( size_t i=0; i < Self::Dimension; ++i ) 
      {
      this->m_Matrix[i][j] *= scalar;
      }
    }
  return *this;
}

template<size_t NDIM,class TSCALAR>
const FixedSquareMatrix<NDIM,TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>::operator*
( const Self& other ) const
{
  Self result( NULL );

  for ( size_t j=0; j < Self::Dimension; ++j ) 
    {
    for ( size_t i=0; i < Self::Dimension; ++i ) 
      {
      result[i][j] = this->m_Matrix[i][0] * other.m_Matrix[0][j];
      for ( size_t k=1; k < Self::Dimension; ++k )
	result[i][j] += this->m_Matrix[i][k] * other.m_Matrix[k][j];
      }
    }
  
  return result;
}

template<size_t NDIM,class TSCALAR>
FixedSquareMatrix<NDIM,TSCALAR>& 
FixedSquareMatrix<NDIM,TSCALAR>::operator=( const Self& other )
{
  memcpy( this->m_Matrix, other.m_Matrix, sizeof( this->m_Matrix ) );
  return *this;
}

template<size_t NDIM,class TSCALAR>
TSCALAR
FixedSquareMatrix<NDIM,TSCALAR>::FrobeniusNorm() const
{
  typename Self::ScalarType norm = DataTypeTraits<TSCALAR>::Zero();
  for ( size_t i = 0; i < Self::Dimension; ++i ) 
    {
    for ( size_t j = 0; j < Self::Dimension; ++j )
      norm += MathUtil::Square( this->m_Matrix[i][j] );
    }
  return sqrt( norm );
}

} // namespace cmtk
