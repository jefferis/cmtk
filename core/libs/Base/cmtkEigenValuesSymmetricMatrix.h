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

#ifndef __cmtkEigenValuesSymmetricMatrix_h_included_
#define __cmtkEigenValuesSymmetricMatrix_h_included_

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
class EigenValuesSymmetricMatrix
{
public:
  /// This class.
  typedef EigenValuesSymmetricMatrix<TFloat> Self;

  /// Constructor: compute eigensystem of given matrix.
  EigenValuesSymmetricMatrix( const SymmetricMatrix<TFloat>& matrix /*!< Symmetric  matrix for which we are computing the eigenvalues and eigenvectors.*/ );
  
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

protected:
  /// Eigenvalues vector.
  std::vector<TFloat> m_Eigenvalues;
};

//@}

} // namespace cmtk

#include "cmtkEigenValuesSymmetricMatrix.txx"

#endif // #ifndef __cmtkEigenValuesSymmetricMatrix_h_included_
