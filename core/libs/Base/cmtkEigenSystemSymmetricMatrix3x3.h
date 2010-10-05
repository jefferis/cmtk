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

#ifndef __cmtkEigenSystemSymmetricMatrix3x3_h_included_
#define __cmtkEigenSystemSymmetricMatrix3x3_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMatrix3x3.h>
#include <Base/cmtkFixedVector.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/** Compute the eigenvectors and eigenvalues of a symmetric 3x3 matrix.
   *  (Eigen decomposition code for symmetric 3x3 matrices, copied from the public
   *   domain Java Matrix library JAMA by Connelly Barnes. ) Eigenvectors and eigenvalues
   *  are returned in sorted order, by ascending absolute values of the eigenvalues.
   */  
template<class TFloat>
class EigenSystemSymmetricMatrix3x3
{
public:
  /// Constructor: compute eigensystem of given matrix.
  EigenSystemSymmetricMatrix3x3( const Matrix3x3<TFloat>& matrix, /*!< Symmetric 3x3 matrix for which we are computing the eigenvalues and eigenvectors.*/ 
				 const bool sortAbsolute = true /*!< Flag for sorting by absolute eigenvalues (default) vs. sorting by actual eigenvalues.*/ );
  
  /// Get n-th eigenvector.
  const FixedVector<3,TFloat> GetNthEigenvector( const size_t n ) const
  {
    return FixedVector<3,TFloat>( this->m_Eigenvectors[n] );
  }
  
  /// Get n-th eigenvalue.
  TFloat GetNthEigenvalue( const size_t n ) const
  {
    return this->m_Eigenvalues[n];
  }

protected:
  /// Eigenvector matrix.
  TFloat m_Eigenvectors[3][3];

  /// Eigenvalues vector.
  TFloat m_Eigenvalues[3];

private:
  /// Helper function that computes the Euclidean length of (x,y).
  static TFloat hypot2( const TFloat& x, const TFloat& y);

  /** Symmetric Householder reduction to tridiagonal form.
   */
  static void tred2(TFloat V[3][3], TFloat d[3], TFloat e[3]);

  /* Symmetric tridiagonal QL algorithm.
   */
  static void tql2(TFloat V[3][3], TFloat d[3], TFloat e[3], const bool sortAbsolute = true );
};

//@}

} // namespace cmtk

#include "cmtkEigenSystemSymmetricMatrix3x3.txx"

#endif // #ifndef __cmtkEigenSystemSymmetricMatrix3x3_h_included_
