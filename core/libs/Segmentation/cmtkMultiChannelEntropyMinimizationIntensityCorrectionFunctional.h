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
//  $Revision: 3517 $
//
//  $LastChangedDate: 2011-10-27 12:42:15 -0700 (Thu, 27 Oct 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkMultiChannelEntropyMinimizationIntensityCorrectionFunctional_h_included_
#define __cmtkMultiChannelEntropyMinimizationIntensityCorrectionFunctional_h_included_

#include <cmtkconfig.h>

#include <Segmentation/cmtkEntropyMinimizationIntensityCorrectionFunctional.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkDataGrid.h>

#include <System/cmtkSmartPtr.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{
/** Functional to correct MR intensity bias by miniming image entropy.
 */
template<unsigned int NOrderAdd,unsigned int NOrderMul>
class MultiChannelEntropyMinimizationIntensityCorrectionFunctional :
  /// Inherit from base class.
  public Functional
{
public:
  /// This class type.
  typedef MultiChannelEntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul> Self;
  
  /// Class for single-channel functional type.
  typedef EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul> SingleChannelFunctionalType;

  /// Pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass type.
  typedef Functional Superclass;

  /// Constructor.
  MultiChannelEntropyMinimizationIntensityCorrectionFunctional( std::vector<UniformVolume::SmartConstPointer>& inputImages ) 
    : m_InputImages( inputImages )
  {
    this->m_SingleChannelFunctionals.resize( this->m_InputImages.size() );
    for ( size_t idx = 0; idx < this->m_SingleChannelFunctionals.size(); ++idx )
      {
      this->m_SingleChannelFunctionals[idx] = SingleChannelFunctionalType::SmartPtr( new SingleChannelFunctionalType );
      this->m_SingleChannelFunctionals[idx]->SetInputImage( this->m_InputImages[idx] );
      }
  }

  /// Set foreground mask for all input images.
  virtual void SetForegroundMask( const UniformVolume& foregroundMask )
  {
    for ( size_t idx = 0; idx < this->m_SingleChannelFunctionals.size(); ++idx )
      {
      this->m_SingleChannelFunctionals[idx]->SetForegroundMask( foregroundMask );
      }
  }

  /// Virtual destructor.
  virtual ~MultiChannelEntropyMinimizationIntensityCorrectionFunctional() {}

  /// Return parameter vector length.
  virtual size_t ParamVectorDim() const
  {
    return SingleChannelFunctionalType::NumberOfParameters * this->m_InputImages.size();
  }

  /// Return parameter stepping.
#pragma GCC diagnostic ignored "-Wtype-limits"
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    // delegate to corresponding single-channel functional, although they are really all the same.
    return this->m_SingleChannelFunctionals[idx / SingleChannelFunctionalType::NumberOfParameters]->GetParamStep( idx % SingleChannelFunctionalType::NumberOfParameters );
  }

  /// Copy parameters to the two correction polynomials.
  virtual void SetParamVector( CoordinateVector& v )
  {
    this->m_ParameterVector = v;
    for ( size_t idx = 0; idx < this->m_SingleChannelFunctionals.size(); ++idx )
      {
      }
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

private:
  /// Keep a copy of the current parameter vector.
  CoordinateVector m_ParameterVector;

  /// Vector of input images.
  std::vector<UniformVolume::SmartConstPtr> m_InputImages;

  /// Vector of single-channel functionals.
  std::vector<typename SingleChannelFunctionalType::SmartPtr> m_SingleChannelFunctionals;
};

/// Create functional templated over polynomial degrees.
template<unsigned int NDegreeMul>
MultiChannelEntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateMultiChannelEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd );

/// Create functional templated over polynomial degrees.
MultiChannelEntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateMultiChannelEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul );

/** Create functional templated over polynomial degrees with initialization from old functional.
 * This function creates a new functional and copies the polynomial coefficients from an existing
 * functional of equal or lower polynomial degrees into the correct locations of the new functional's
 * parameter vector. This is for incremental computation.
 */
MultiChannelEntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateMultiChannelEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul,
  MultiChannelEntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr oldFunctional );

//@}

} // namespace cmtk

#include "cmtkMultiChannelEntropyMinimizationIntensityCorrectionFunctional.txx"

#endif // #ifndef __cmtkMultiChannelEntropyMinimizationIntensityCorrectionFunctional_h_included_

