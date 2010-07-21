/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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
::QRDecomposition( const Matrix2D<TFloat>& matrix )
{
  m = matrix.GetNumberOfRows();
  n = matrix.GetNumberOfColumns();

  /* Copy matrix into compactQR
   */

  compactQR.setbounds(0, (int)matrix.GetNumberOfRows(), 0, (int)matrix.GetNumberOfColumns());
  for ( int j = 0; j < m; j++ )
    for ( int i = 0; i < n; i++ )
      compactQR(i,j) = (double)(1.0 * matrix[i][j]);
  
  /* Run AlgLib QR decomposition
   */
  rmatrixqr( compactQR, m, n, tau );

  /* Initialized Q and R objects
   */
  Q = matrixPtr ( new matrix2D( m, n ) );
  R = matrixPtr ( new matrix2D( m, n ) );

  /* Set these to true once the Q or R 
   * matrix has been extracted from the compact
   * alglib representation
   */ 
  extractedQ = extractedR = false;

}

/// Get the Q factor 
template<class TFloat>
SmartPointer< Matrix2D<TFloat> >
QRDecomposition<TFloat>
::GetQ() 
{
  if ( ! extractedQ ) 
    {
    /* Extract Q from compactQR
      */
    ap::real_2d_array tmp_ap_matrix;
    rmatrixqrunpackq( compactQR, m, n, tau, n, tmp_ap_matrix );

    for ( int j = 0; j < m; j++ )
      for ( int i = 0; i < n; i++ )
        (*Q)[i][j] = tmp_ap_matrix(i,j);

    extractedQ = true;
    }
  return this->Q;
}

/// Get the R factor 
template<class TFloat>
SmartPointer< Matrix2D<TFloat> >
QRDecomposition<TFloat>
::GetR()
{
  if ( ! extractedR ) 
    {
    /* Extract R from compactQR
      */ 
    ap::real_2d_array tmp_ap_matrix;
    rmatrixqrunpackr( compactQR, m, n, tmp_ap_matrix );
    for ( int j = 0; j < m; j++ )
      for ( int i = 0; i < n; i++ )
        (*R)[i][j] = tmp_ap_matrix(i,j);

    extractedR = true;
    }

  return this->R;
}

} // namespace cmtk
