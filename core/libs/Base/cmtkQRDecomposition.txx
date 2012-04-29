/*
//
//  Copyright 1997-2012 Torsten Rohlfing
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

#include "Numerics/sevd.h"
#include "Numerics/qr.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class TFloat>
QRDecomposition<TFloat>
::QRDecomposition( const typename Self::MatrixType& matrix )
{
  this->m_Rows = matrix.NumberOfRows();
  this->m_Cols = matrix.NumberOfColumns();

  /* Copy matrix into this->m_CompactQR */
  this->m_CompactQR.setbounds(0, static_cast<int>( this->m_Rows-1 ), 0, static_cast<int>( this->m_Cols-1 ) );
  for ( size_t j = 0; j < this->m_Rows; j++ )
    for ( size_t i = 0; i < this->m_Cols; i++ )
      this->m_CompactQR(i,j) = static_cast<double>( matrix[i][j] );
  
  /* Run AlgLib QR decomposition */
  rmatrixqr( this->m_CompactQR, this->m_Rows, this->m_Cols, this->m_Tau );
}

template<class TFloat> template<size_t NDIM> 
QRDecomposition<TFloat>
::QRDecomposition( const FixedSquareMatrix<NDIM,TFloat>& matrix )
{
  this->m_Rows = this->m_Cols = NDIM;

  /* Copy matrix into this->m_CompactQR */
  this->m_CompactQR.setbounds(0, static_cast<int>( this->m_Rows-1 ), 0, static_cast<int>( this->m_Cols-1 ) );
  for ( size_t j = 0; j < this->m_Rows; j++ )
    for ( size_t i = 0; i < this->m_Cols; i++ )
      this->m_CompactQR(i,j) = static_cast<double>( matrix[i][j] );
  
  /* Run AlgLib QR decomposition */
  rmatrixqr( this->m_CompactQR, this->m_Rows, this->m_Cols, this->m_Tau );
}

/// Get the Q factor 
template<class TFloat>
Matrix2D<TFloat>&
QRDecomposition<TFloat>
::GetQ() 
{
  if ( ! this->m_Q ) 
    {
    this->m_Q = typename Self::MatrixType::SmartPtr ( new typename Self::MatrixType( this->m_Rows, this->m_Cols ) );

    /* Extract Q from this->m_CompactQR */
    ap::real_2d_array tmp_ap_matrix;
    rmatrixqrunpackq( this->m_CompactQR, this->m_Rows, this->m_Cols, this->m_Tau, this->m_Cols, tmp_ap_matrix );

    for ( int j = 0; j < this->m_Rows; j++ )
      for ( int i = 0; i < this->m_Cols; i++ )
        (*this->m_Q)[i][j] = tmp_ap_matrix(i,j);
    }
  return *(this->m_Q);
}

/// Get the R factor 
template<class TFloat>
Matrix2D<TFloat>&
QRDecomposition<TFloat>
::GetR()
{
  if ( ! this->m_R ) 
    {
    this->m_R = typename Self::MatrixType::SmartPtr ( new typename Self::MatrixType( this->m_Rows, this->m_Cols ) );
    
    /* Extract R from this->m_CompactQR */ 
    ap::real_2d_array tmp_ap_matrix;
    rmatrixqrunpackr( this->m_CompactQR, this->m_Rows, this->m_Cols, tmp_ap_matrix );
    for ( size_t j = 0; j < this->m_Rows; j++ )
      for ( size_t i = 0; i < this->m_Cols; i++ )
        (*this->m_R)[i][j] = tmp_ap_matrix(i,j);
    }

  return *(this->m_R);
}

} // namespace cmtk
