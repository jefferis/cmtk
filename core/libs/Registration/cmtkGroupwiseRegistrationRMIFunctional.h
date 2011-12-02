/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkGroupwiseRegistrationRMIFunctional_h_included_
#define __cmtkGroupwiseRegistrationRMIFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkGroupwiseRegistrationFunctionalXformTemplate.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkThreads.h>
#include <System/cmtkMutexLock.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkXform.h>
#include <Base/cmtkMathUtil.h>
#include <Base/cmtkMatrix.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for groupwise registration.
 */
template<class TXform>
class GroupwiseRegistrationRMIFunctional : 
  public GroupwiseRegistrationFunctionalXformTemplate<TXform>
{
public:
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalXformTemplate<TXform> Superclass;

  /// Type of this class.
  typedef GroupwiseRegistrationRMIFunctional<TXform> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Transformation type.
  typedef TXform XformType;

  /// Smart pointer to transformation type.
  typedef typename XformType::SmartPtr XformPointer;

  /// Constructor.
  GroupwiseRegistrationRMIFunctional();

  /// Destructor.
  virtual ~GroupwiseRegistrationRMIFunctional();

  /** Set template grid.
   */
  virtual void SetTemplateGrid( UniformVolume::SmartPtr& templateGrid /*!< The template grid that defines size and resolution for the implicit registration template.*/, 
				const int downsample = 1 /*!< Downsampling factor */, const bool useTemplateData = false /*!< Flag to use template pixel data, not just grid, in registration */ );

  /** Compute functional value and gradient.
   *\param v Parameter vector.
   *\param g The extimated gradient of the functional is stored in this vector.
   *\param step Step size for finite difference gradient approximation. Default
   *  is 1 mm.
   *\return Const function value for given parameters.
   */
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 );

  /// Evaluate functional with currently set parameters.
  virtual typename Self::ReturnType Evaluate();
  
protected:
  /// Covariance matrix type.
  typedef Matrix2D<typename Self::ReturnType> CovarianceMatrixType;

  /// Actual covariance matrix.
  typename Self::CovarianceMatrixType m_CovarianceMatrix;

  /// Type for vectors of sums and products.
  typedef std::vector<long int> SumsAndProductsVectorType;

  /// Sum of products matrix.
  SumsAndProductsVectorType m_SumOfProductsMatrix;

  /// Sums vector.
  SumsAndProductsVectorType m_SumsVector;

  /// Compute metric from partial matrices using temporary matrix storage.
  typename Self::ReturnType GetMetric
  ( const SumsAndProductsVectorType& sumOfProductsMatrix, 
    const SumsAndProductsVectorType& sumsVector,
    const unsigned int totalNumberOfSamples,
    typename Self::CovarianceMatrixType& covarianceMatrix ) const;

  /// Sum of products matrix.
  std::vector<SumsAndProductsVectorType> m_ThreadSumOfProductsMatrix;

  /// Sums vector.
  std::vector<SumsAndProductsVectorType> m_ThreadSumsVector;

  /** Total number of samples that went into CC computation.
   */
  unsigned int m_TotalNumberOfSamples;

  /// Mutex for schared data structures.
  MutexLock m_MutexLock;

  /// Update probabilistic sample table..
  virtual bool Wiggle();

private:
  /// Thread parameters with no further data.
  typedef ThreadParameters<Self> ThreadParametersType;

  /// Thread parameter for entropy evaluation.
  class EvaluateThreadParameters : 
    /// Inherit from generic thread parameter class.
    public ThreadParametersType
  {
  };
  
  /// Evaluate functional with currently set parameters.
  static void EvaluateThread( void *const threadParameters, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Evaluate functional with currently set parameters with probabilistic sampling.
  static void EvaluateProbabilisticThread( void *const threadParameters, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkGroupwiseRegistrationRMIFunctional_h_included_
