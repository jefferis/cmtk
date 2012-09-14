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
//  $Revision: 4427 $
//
//  $LastChangedDate: 2012-06-12 11:23:22 -0700 (Tue, 12 Jun 2012) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkEigenSystemSymmetricMatrix_h_included_
#define __cmtkEigenSystemSymmetricMatrix_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkSymmetricMatrix.h>
#include <Base/cmtkVector.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/** Compute the eigenvectors and eigenvalues of a symmetric, square matrix of arbitrary size.
  */  
template<class TFloat>
class EigenSystemSymmetricMatrix
{
public:
  /// This class.
  typedef EigenSystemSymmetricMatrix<TFloat> Self;

  /// Constructor: compute eigensystem of given matrix.
  EigenSystemSymmetricMatrix( const SymmetricMatrix<TFloat>& matrix /*!< Symmetric  matrix for which we are computing the eigenvalues and eigenvectors.*/ );
  
  /// Get n-th eigenvector.
  const Vector<TFloat> GetNthEigenvector( const size_t n ) const
  {
    return Vector<TFloat>( this->m_Eigenvectors[n] );
  }
  
  /// Get n-th eigenvector.
  const TFloat EigenvectorElement( const size_t n, const size_t i ) const
  {
    return this->m_Eigenvectors[n][i];
  }
  
  /// Get n-th eigenvalue.
  TFloat GetNthEigenvalue( const size_t n ) const
  {
    return this->m_Eigenvalues[n];
  }

  /// Get vector of eigenvalues.
  std::vector<TFloat> GetEigenvalues() const
  {
    return this->m_Eigenvalues;
  }

private:
  /// Eigenvector matrix.
  std::vector< Vector<TFloat> > m_Eigenvectors;

  /// Eigenvalues vector.
  std::vector<TFloat> m_Eigenvalues;
};

//@}

} // namespace cmtk

#include "cmtkEigenSystemSymmetricMatrix.txx"

#endif // #ifndef __cmtkEigenSystemSymmetricMatrix_h_included_
