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

#ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctional_h_included_
#define __cmtkEntropyMinimizationIntensityCorrectionFunctional_h_included_

#include <cmtkconfig.h>

#include <Segmentation/cmtkEntropyMinimizationIntensityCorrectionFunctionalBase.h>

#include <Base/cmtkPolynomial.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkDataGrid.h>

#include <System/cmtkThreadPool.h>
#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{
/** Functional to correct MR intensity bias by miniming image entropy.
 */
template<unsigned int NOrderAdd,unsigned int NOrderMul>
class EntropyMinimizationIntensityCorrectionFunctional :
  /// Inherit from base class.
  public EntropyMinimizationIntensityCorrectionFunctionalBase
{
public:
  /// This class type.
  typedef EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul> Self;
  
  /// Pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass type.
  typedef EntropyMinimizationIntensityCorrectionFunctionalBase Superclass;

  /// Additive correction polynomial.
  typedef Polynomial<NOrderAdd,Types::Coordinate> PolynomialTypeAdd;

  /// Multiplicative correction polynomial.
  typedef Polynomial<NOrderMul,Types::Coordinate> PolynomialTypeMul;

  /// Constructor.
  EntropyMinimizationIntensityCorrectionFunctional() 
  {
    this->m_ParameterVector.SetDim( this->ParamVectorDim() );
    this->m_ParameterVector.SetAll( 0.0 );

    ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
    this->m_NumberOfThreads = threadPool.GetNumberOfThreads();
    this->m_MonomialsPerThread = std::max( (int)PolynomialTypeAdd::NumberOfMonomials, (int)PolynomialTypeMul::NumberOfMonomials );
    this->m_MonomialsVec = Memory::AllocateArray<Types::Coordinate>( this->m_NumberOfThreads * this->m_MonomialsPerThread  );
  }

  /// Virtual destructor.
  virtual ~EntropyMinimizationIntensityCorrectionFunctional() 
  {
    Memory::DeleteArray( this->m_MonomialsVec );
  }

  /// Total number of parameters is number of additive plus number of multiplicative parameters.
  static const size_t m_NumberOfParameters = PolynomialTypeAdd::NumberOfMonomials + PolynomialTypeMul::NumberOfMonomials;

  /// Return parameter vector length.
  virtual size_t ParamVectorDim() const
  {
    return this->m_NumberOfParameters;
  }

  /// Return parameter stepping.
#pragma GCC diagnostic ignored "-Wtype-limits"
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    if ( idx < PolynomialTypeAdd::NumberOfMonomials )
      return (this->m_InputImageRange / 256) * this->m_StepSizeAdd[idx] * mmStep;
    else
      return (this->m_InputImageRange / 256) * this->m_StepSizeMul[idx-PolynomialTypeAdd::NumberOfMonomials] * mmStep;
  }

  /// Get number of additive monomials.
  virtual size_t GetNumberOfMonomialsAdd() const 
  {
    return PolynomialTypeAdd::NumberOfMonomials;
  }

  /// Get number of multiplicative monomials.
  virtual size_t GetNumberOfMonomialsMul() const 
  {
    return PolynomialTypeMul::NumberOfMonomials;
  }

  /// Copy parameters to the two correction polynomials.
  virtual void SetParamVector( CoordinateVector& v )
  {
    this->m_ParameterVector = v;
    for ( int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
      this->m_CoefficientsAdd[i] = v[i] * this->m_MulCorrectionAdd[i];
    
    size_t ofs = PolynomialTypeAdd::NumberOfMonomials;
    for ( int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i, ++ofs )
      this->m_CoefficientsMul[i] = v[ofs] * this->m_MulCorrectionMul[i];
  }

  /// Extract parameter vector from the two correction polynomials.
  virtual void GetParamVector( CoordinateVector& v )
  {
    v = this->m_ParameterVector;
  }

  /** Fast implementation of gradient evaluation.
   * This function only updates either the additive of the multiplicative bias 
   * field, depending on which gradient component is being determined.
   */
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step );

protected:
  /// Keep a copy of the current parameter vector.
  CoordinateVector m_ParameterVector;

  /// Additive correction term coefficients.
  Types::Coordinate m_StepSizeAdd[1+PolynomialTypeAdd::NumberOfMonomials];
  
  /// Multiplicative correction term coefficients.
  Types::Coordinate m_StepSizeMul[1+PolynomialTypeMul::NumberOfMonomials];

  /// Additive correction term coefficients.
  Types::Coordinate m_CoefficientsAdd[1+PolynomialTypeAdd::NumberOfMonomials];
  
  /// Multiplicative correction term coefficients.
  Types::Coordinate m_CoefficientsMul[1+PolynomialTypeMul::NumberOfMonomials];

  /// Additive correction constants for additive polynomials.
  Types::Coordinate m_AddCorrectionAdd[1+PolynomialTypeAdd::NumberOfMonomials];

  /// Multiplicative correction constants for additive polynomials.
  Types::Coordinate m_MulCorrectionAdd[1+PolynomialTypeAdd::NumberOfMonomials];

  /// Additive correction constants for multiplicative polynomials.
  Types::Coordinate m_AddCorrectionMul[1+PolynomialTypeMul::NumberOfMonomials];

  /// Multiplicative correction constants for additive polynomials.
  Types::Coordinate m_MulCorrectionMul[1+PolynomialTypeMul::NumberOfMonomials];

private:
  /** Number of parallel threads.
   * Because the constructor allocated temporary storage for monomial
   * computation, we need to know in advance how many threads there will
   * be, and make sure we never exceed this number later when calling the
   * thread functions.
   */
  size_t m_NumberOfThreads;

  /// Maximum number of monomials per thread (maximum of additive and multiplicative polynomial number of monomials).
  size_t m_MonomialsPerThread;

  /// Temporary storage for evaluating monomials.
  Types::Coordinate* m_MonomialsVec;

  /// Update polynomial correctionfactors from input image.
  virtual void UpdateCorrectionFactors();

  /// Jointly update both bias images.
  virtual void UpdateBiasFields( const bool foregroundOnly = true );

  /// Thread function: jointly update both bias images.
  static void UpdateBiasFieldsThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Thread function: jointly update both bias images.
  static void UpdateBiasFieldsAllThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Update additive bias image.
  virtual void UpdateBiasFieldAdd( const bool foregroundOnly = true );

  /// Thread function: update foreground additive bias images.
  static void UpdateBiasFieldAddThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Thread function: update complete additive bias images.
  static void UpdateBiasFieldAddAllThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Update multiplicative bias image.
  virtual void UpdateBiasFieldMul( const bool foregroundOnly = true );

  /// Thread function: update foreground multiplicative bias images.
  static void UpdateBiasFieldMulThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Thread function: update complete multiplicative bias images.
  static void UpdateBiasFieldMulAllThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );
};

/// Create functional templated over polynomial degrees.
template<unsigned int NDegreeMul>
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd );

/// Create functional templated over polynomial degrees.
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul );

/** Create functional templated over polynomial degrees with initialization from old functional.
 * This function creates a new functional and copies the polynomial coefficients from an existing
 * functional of equal or lower polynomial degrees into the correct locations of the new functional's
 * parameter vector. This is for incremental computation.
 */
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul,
  EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr oldFunctional );

//@}

} // namespace cmtk

#include "cmtkEntropyMinimizationIntensityCorrectionFunctional.txx"

#endif // #ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctional_h_included_

