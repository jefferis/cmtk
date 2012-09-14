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

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkMatrix.h>

#include <System/cmtkConsole.h>

#include "Numerics/sevd.h"
#include "Numerics/spddet.h"
#include "Numerics/svd.h"

#include <math.h>
#include <algorithm>
#include <vector>

namespace
cmtk
{

void 
MathUtil::SVD( Matrix2D<double>& U, std::vector<double>& W, Matrix2D<double>& V )
{
  const size_t m = U.NumberOfRows();
  const size_t n = U.NumberOfColumns();

  W.resize( n );
  V.Resize( n, n );

  ap::real_2d_array apA;
  apA.setbounds(0, m-1, 0, n-1);
  for (size_t j = 0; j < n; j++)
    for (size_t i = 0; i < m; i++)
      apA(i,j) = U[i][j];

  ap::real_1d_array w;
  ap::real_2d_array u;
  ap::real_2d_array vt;

  rmatrixsvd( apA, m, n, 
              true /* U needed */, 
              true /* V needed */, 
              2 /*max-level memory usage */, 
              w, u, vt);

  /* Put u in U */
  for (size_t j = 0; j < n; j++)
    for (size_t i = 0; i < m; i++)
      U[i][j] = u(i,j);
  
  /* Put w in W */
  for (size_t i = 0; i < n; i++)
    W[i] = w(i);
  
  /* Un-transpose vt and put it in V */
  for (size_t j = 0; j < n; j++)
    for (size_t i = 0; i < n; i++)
      V[i][j] = vt(j,i);
}

/** TODO: move this someplace more logical than the linear-algebra module
 */
void
MathUtil::SVDLinearRegression( const Matrix2D<double>& U, const std::vector<double>& W, const Matrix2D<double>& V, const std::vector<double>& b, std::vector<double>& lm_params )
{
  const size_t m = U.NumberOfRows();
  const size_t n = U.NumberOfColumns();
  
  lm_params.resize( n );

  // From alglib linear regression:
  // Take the inverses of the singular values, setting the inverse 
  // to 0 if the sv is close to 0 (tolerance controlled by epstol)
  double epstol = 1000;
  ap::real_1d_array svi;
  svi.setbounds( 0, n-1 );
  for( size_t i = 0; i < n; i++ )
    if( W[i] > epstol*ap::machineepsilon * W[0] )
      svi(i) = 1 / W[i];
    else
      svi(i) = 0;
  
  // Calculate linear model parameters following Heath, Ch. 3.6
  // (Scientific Computing: An Introductory Survey, 2nd Ed., 2002)
  for ( size_t i = 0; i < n; i++ )
    lm_params[i] = 0.0;
  double ut_times_b; 
  
  for ( size_t i = 0; i < n; i++ )
    {
    ut_times_b = 0.0;
    for ( size_t j = 0; j < m; j++ )
      ut_times_b += U[j][i] * b[j];
    ut_times_b *= svi(i);
    for ( size_t j = 0; j < n; j++ )
      lm_params[j] += ut_times_b * V[j][i]; 
    }
}

/////////////////////////////////////////////////////////////////////
// HELPERS
/////////////////////////////////////////////////////////////////////

template<class T> 
T
MathUtil::CholeskyDeterminant
(const Matrix2D<T>& matrix, int n)
{    
  ap::real_2d_array apMatrix;
  apMatrix.setbounds(0, n-1, 0, n-1);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      apMatrix(i,j) = (double)(1.0 * matrix[i][j]);
  spdmatrixcholesky( apMatrix, n, false );
  T determinant = (T) spdmatrixcholeskydet( apMatrix, n );
  return determinant;
}

template double MathUtil::CholeskyDeterminant<double>(const Matrix2D<double>& matrix, int n);
template float MathUtil::CholeskyDeterminant<float>(const Matrix2D<float>& matrix, int n);

} // namespace cmtk

