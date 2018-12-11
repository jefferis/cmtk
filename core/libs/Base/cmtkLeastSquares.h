/*
//
//  Copyright 2012 SRI International
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

#ifndef __cmtkLeastSquares_h_included_
#define __cmtkLeastSquares_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMatrix.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Solve linear least-squares problem.
 * This class solves a least-squares problem of the form Ax=y for x, where A is a n-by-m design matrix, y is an n-dimensional measurement vector, 
 * and x is the m-dimensional parameter vector.
 *
 * Solution of the least-squares problem is implemented via Singular Value Decomposition of the design matrix. Multiple problems using the same
 * design matrix but different measurement vectors can be solved without repeating the SVD.
 */
template<class TScalar>
class LeastSquares
{
public:
  /// This class.
  typedef LeastSquares<TScalar> Self;

  /// The scalar data type.
  typedef TScalar ScalarType;

  /// Constructor.
  LeastSquares( const Matrix2D<typename Self::ScalarType>& designMatrix )
    : m_MatrixU( designMatrix ),
      m_MatrixV( designMatrix.NumberOfColumns(), designMatrix.NumberOfColumns() ),
      m_VectorW( designMatrix.NumberOfColumns() )
  {
    MathUtil::SVD( this->m_MatrixU, this->m_VectorW, this->m_MatrixV );
  }

  /** Compute parameter vector for a given measurement data vector.
   *\return Parameter vector that minimizes the least squares fitting error.
   */
  std::vector<typename Self::ScalarType> Solve( const std::vector<typename Self::ScalarType>& measurements /*!< The n-dimensional measurement data vector.*/ ) const
  {
    std::vector<typename Self::ScalarType> parameters( this->m_MatrixU.NumberOfRows() );
    MathUtil::SVDLinearRegression( this->m_MatrixU, this->m_VectorW, this->m_MatrixV, measurements, parameters );
    return parameters;
  }
  
private:
  /// Matrix U returned from singular value decomposition.
  Matrix2D<typename Self::ScalarType> m_MatrixU;

  /// Matrix V returned from singular value decomposition.
  Matrix2D<typename Self::ScalarType> m_MatrixV;

  /// Vector W returned from singular value decomposition.
  std::vector<typename Self::ScalarType> m_VectorW;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLeastSquares_h_included_
