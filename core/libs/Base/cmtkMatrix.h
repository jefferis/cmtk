/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkMatrix_h_included_
#define __cmtkMatrix_h_included_

#include <cmtkconfig.h>

#include <string.h>
#include <vector>
#include <iostream>

#include <cmtkSmartPtr.h>
#include <cmtkMemory.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/// Rekursive matrix template.
template<class TElement,size_t NDim>
class Matrix
{
public:
  /// This class.
  typedef Matrix<TElement,NDim> Self;

  /// Superclass.
  typedef Matrix<TElement,NDim-1> Superclass;

  /// Public constructor.
  Matrix( const size_t (&dims)[NDim] )
    : m_SubMatrixArray( dims[0] )
  {
    size_t nItems = dims[0];
    for ( size_t i = 1; i < NDim; ++i )
      nItems *= dims[i];

    TElement* data = Memory::AllocateArray<TElement>( nItems );
  }

  /// Destructor.
  ~Matrix() {};

  /// Element pointer type.
  typedef typename Superclass::ElementPointerType* ElementPointerType;

  typename Self::ElementPointerType operator[]( const size_t idx )
  {
    return this->m_SubMatrixArray[idx];
  }
  
  const typename Self::ElementPointerType operator[]( const size_t idx ) const
  {
    return this->m_SubMatrixArray[idx];
  }

protected:
  /// Recursive constructor.
  Matrix() {};
  
private:
  /// Vector of pointers to lower-dimensional sub-matrices.
  std::vector<typename Self::ElementPointerType> m_SubMatrixArray;
}; // class Matrix

template<class TElement>
class Matrix<TElement,1>
{
};

/// Two-dimensional matrix template.
template<class T>
class Matrix2D :
  /// For access, make this a vector of pointers.
  public std::vector<T*>
{
public:
  /// Superclass.
  typedef std::vector<T*> Superclass;

  /// This class.
  typedef Matrix2D<T> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Row vector type.
  typedef std::vector<T*> RowVectorType;  
  
  /// Default constructor.
  Matrix2D()
    : Superclass( 1 )
  { 
    this->m_NumberOfColumns = 0;
    this->m_NumberOfRows = 0;
    this->m_NumberOfElements = 0;
    
    (*this)[0] = NULL;
  }
  
  /// Constructor: allocate and create cross-references.
  Matrix2D( const size_t dims1, const size_t dims0, const T* data = NULL )
    : Superclass( dims1 )
  { 
    this->m_NumberOfColumns = dims0;
    this->m_NumberOfRows = dims1;
    this->m_NumberOfElements = dims0 * dims1;

    (*this)[0] = Memory::AllocateArray<T>(  this->m_NumberOfElements  );
    for ( size_t i = 1; i < this->m_NumberOfRows; ++i )
      (*this)[i] = (*this)[i-1] + this->m_NumberOfColumns;
    
    if ( data )
      memcpy( (*this)[0], data, this->m_NumberOfElements * sizeof( T ) );
  }
  
  /// Copy constructor.
  Matrix2D( const Matrix2D<T>& other ) :
    Superclass( other.size() )
  { 
    this->m_NumberOfColumns = other.m_NumberOfColumns;
    this->m_NumberOfRows = other.m_NumberOfRows;
    this->m_NumberOfElements = other.m_NumberOfElements;

    (*this)[0] = Memory::AllocateArray<T>(  this->m_NumberOfElements  );
    for ( size_t i = 1; i < this->m_NumberOfRows; ++i )
      (*this)[i] = (*this)[i-1] + this->m_NumberOfColumns;
    
    memcpy( (*this)[0], other[0], this->m_NumberOfElements * sizeof( T ) );
  }

  /// Destructor: free allocated array.
  ~Matrix2D() 
  {
    if ( (*this)[0] )
      {
      delete[] (*this)[0];
      (*this)[0] = NULL;
      }
  }

  /// Get number of rows.
  size_t GetNumberOfRows() const 
  { 
    return this->m_NumberOfRows;
  }

  /** Get number of columns.
   * Get this from underlying Array.
   */
  size_t GetNumberOfColumns() const 
  { 
    return this->m_NumberOfColumns; 
  }

  /// Resize the matrix.
  void Resize( const size_t numberOfRows, const size_t numberOfColumns )
  {
    if ( (numberOfColumns != this->m_NumberOfColumns) ||
	 (numberOfRows != this->m_NumberOfRows) )
      {
      if ( (*this)[0] )
	{
	delete[] (*this)[0];
	(*this)[0] = NULL;
	}
      
      this->m_NumberOfColumns = numberOfColumns;
      this->m_NumberOfRows = numberOfRows;
      this->m_NumberOfElements = numberOfColumns * numberOfRows;
      
      this->Superclass::resize( numberOfRows );
      (*this)[0] = Memory::AllocateArray<T>(  this->m_NumberOfElements  );
      for ( size_t i = 1; i < numberOfRows; ++i )
	(*this)[i] = (*this)[i-1] + numberOfColumns;
      }
  }
  
  /// Reset all values to zero.
  void SetAllToZero() 
  {
    memset( (*this)[0], 0, this->m_NumberOfElements * sizeof( T ) );
  }
  
  /// Set all values.
  void SetAll( const T value) 
  {
    for ( size_t i = 0; i < this->m_NumberOfElements; ++i ) 
      {
      (*this)[0][i] = value;
      }
  }
  
  /// Copy another matrix.
  Matrix2D<T>& operator= ( const Matrix2D<T>& other ) 
  {
    this->Resize( other.GetNumberOfColumns(), other.GetNumberOfRows() );
    memcpy( (*this)[0], other[0], this->m_NumberOfElements * sizeof( T ) );
    return *this;
  }
  
private:
  /// Size of the allocated array.
  size_t m_NumberOfElements;

  /// Number of rows.
  size_t m_NumberOfColumns;

  /// Number of rows.
  size_t m_NumberOfRows;
};

/// Three-dimensional matrix template.
template<class T>
class Matrix3D :
  /// For access, make this a 2-D matrix of pointers.
  public Matrix2D<T*>
{
public:
  /// This class.
  typedef Matrix3D<T> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef Matrix2D<T*> Superclass;

  /// Constructor: allocate and create cross-references.
  Matrix3D<T>
  ( const size_t dims2, const size_t dims1, const size_t dims0 )
    : Matrix2D<T*>( dims2, dims1 )
  {
    this->m_NumberOfPlanes = dims0;
    this->m_NumberOfElements = dims0 * dims1 * dims2;
    (*this)[0][0] = Memory::AllocateArray<T>( this->m_NumberOfElements );

    for ( size_t j = 0; j < this->GetNumberOfRows(); ++j )
      for ( size_t i = 0; i < this->GetNumberOfColumns(); ++i )
	if ( i && j )
	  {
	  (*this)[i][j] = (*this)[0][0] + this->GetNumberOfRows() * ( i + this->GetNumberOfColumns() * j );
	  }
  }

  /// Return number of planes
  size_t GetNumberOfPlanes() const
  {
    return this->m_NumberOfPlanes;
  }
  
  /// Resize the matrix.
  void Resize( const size_t numberOfRows, const size_t numberOfColumns, const size_t numberOfPlanes )
  {
    if ( ( numberOfColumns != this->GetNumberOfColumns() ) ||
	 ( numberOfRows != this->GetNumberOfRows() ) ||
	 ( numberOfPlanes != this->GetNumberOfPlanes() ) )
      {
      if ( (*this)[0][0] )
	{
	delete[] (*this)[0][0];
	(*this)[0][0] = NULL;
	}
      
      this->m_NumberOfPlanes = numberOfPlanes;
      this->m_NumberOfElements = numberOfPlanes * numberOfRows * numberOfColumns;

      this->Superclass::Resize( numberOfRows, numberOfColumns );
      (*this)[0][0] = Memory::AllocateArray<T>( this->m_NumberOfElements );
      
      for ( size_t j = 0; j < this->GetNumberOfRows(); ++j )
	for ( size_t i = 0; i < this->GetNumberOfColumns(); ++i )
	  if ( i && j )
	    {
	    (*this)[i][j] = (*this)[0][0] + this->GetNumberOfPlanes() * ( i + this->GetNumberOfColumns() * j );
	    }
      }
  }
  
  /// Reset all values to zero.
  void SetAllToZero() 
  {
    memset( (*this)[0][0], 0, this->m_NumberOfElements * sizeof( T ) );
  }
  
  /// Set all values.
  void SetAll( const T value) 
  {
    for ( size_t i = 0; i < this->m_NumberOfElements; ++i ) 
      {
      (*this)[0][0][i] = value;
      }
  }
  
  /// Copy another matrix.
  Matrix2D<T>& operator= ( const Matrix2D<T>& other ) 
  {
    this->Resize( other.GetNumberOfColumns(), other.GetNumberOfRows(), other.GetNumberOfPlanes() );
    memcpy( (*this)[0], other[0], this->m_NumberOfElements * sizeof( T ) );
    return *this;
  }
  
private:
  /// Planes in the 3D matrix.
  size_t m_NumberOfPlanes;

  /// Number of matrix elements.
  size_t m_NumberOfElements;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMatrix_h_included_
