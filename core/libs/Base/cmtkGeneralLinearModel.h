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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkGeneralLinearModel_h_included_
#define __cmtkGeneralLinearModel_h_included_

#include <cmtkconfig.h>

#include <cmtkTypedArray.h>
#include <cmtkMatrix.h>
#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Pixelwise linear modeling and t statistics of data.
 * \note This class formerly contained a method for 
 * getting the covariance matrix of an SVD.  It has 
 * been removed due to obsolete implementation 
 * (Numerical Recipes) and un-use
 */
class GeneralLinearModel
{
public:
  /// This class.
  typedef GeneralLinearModel Self;

  /// Smart pointer type.
  typedef SmartPointer<Self> SmartPtr;

  /** Constructor.
   * This will take care of SVD of the design matrix and perform all necessary
   * pre-computations for the actual modeling.
   */
  GeneralLinearModel( const size_t nParameters, const size_t nData, const double* designMatrix );
  
  /// Destructor.
  ~GeneralLinearModel();

  /** Get singular value of the SVD decomposition of the design matrix.
   *\param n Index of the model parameter [0..NParameters-1].
   *\return The singular value for parameter n.
   */
  double GetSingularValue( const size_t n ) const 
  {
    return (*(this->W))[n];
  }

  /** Get the parameter correlation matrix from design matrix.
   */
  Matrix2D<double>* GetCorrelationMatrix() const;

  /** Model y[] distribution and return model parameters a[].
   *\param y A vector of TypedArray smart pointers. Each object in this
   * vector points to a pixel array from a different subject in a population.
   *\param normalizeParameters If this flag is set (default), then the
   * linear model parameters are normalized w.r.t. the maghnitudes of the
   * respective measurements.
   */
  void FitModel
  ( std::vector<TypedArray::SmartPtr> y, 
    const bool normalizeParameters = true );

  /// Get pointer to n-th model parameter array.
  TypedArray::SmartPtr& GetModel( const size_t n )
  {
    return this->Model[n];
  }

  /// Get normalization factor for parameter number 'p'.
  double GetNormFactor( const size_t p )
  {
    // do not normalize constant part
    if ( this->VariableSD[p] > 0 ) 
      return this->VariableSD[p];
    else
      return 1.0;
  }

  /// Get pointer to n-th parameter t statistics array.
  TypedArray::SmartPtr& GetTStat( const size_t n )
  {
    return this->TStat[n];
  }

  /// Get pointer to n-th parameter t statistics array.
  TypedArray::SmartPtr& GetFStat()
  {
    return this->FStat;
  }

private:
  /// Initialize results arrays with the correct number of pixels.
  void InitResults( const size_t nPixels );

  /// Number of model parameters.
  size_t NParameters;

  /// Number of data items.
  size_t NData;

  /// Design matrix.
  Matrix2D<double> DesignMatrix;

  /// Matrix U of the design matrix SVD.
  Matrix2D<double>* U;

  /// Array of partial design matrices.
  std::vector< Matrix2D<double>* > Up;

  /// Matrix V the design matrix SVD.
  Matrix2D<double>* V;

  /// SVD of partial design matrices.
  std::vector< Matrix2D<double>* > Vp;

  /// Vector W (workspace).
  std::vector<double>* W;

  /// Workspace vectors for partial regressions.
  std::vector< std::vector<double>* > Wp;

  /// Means of variables.
  std::vector<double> VariableMean;

  /// Standard deviations of variables.
  std::vector<double> VariableSD;

  /// Computed model coefficients.
  std::vector<TypedArray::SmartPtr> Model;

  /// Computed model t statistics coefficients.
  std::vector<TypedArray::SmartPtr> TStat;

  /// Computed model F statistics.
  TypedArray::SmartPtr FStat;

  /// Model fitting thread function.
  static CMTK_THREAD_RETURN_TYPE FitModelThreadFunc( void* args );

  /// Thread parameters for model fitting function.
  class FitModelThreadArgs : public ThreadParameters<Self>
  {
  public:
    /// Vector of pointer to input images.
    std::vector<TypedArray::SmartPtr> m_ImageVector;
    
    /// Flag for model parameter normalization.
    bool m_NormalizeParameters;
  };
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkGeneralLinearModel_h_included_
