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

#ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDA_h_included_
#define __cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDA_h_included_

#include <cmtkconfig.h>

#include <cmtkEntropyMinimizationIntensityCorrectionFunctionalBase.h>

#include "cmtkDeviceMemoryCUDA.h"
#include "cmtkUniformVolumeCUDA.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{
/// Base class for GPU implementation entropy-minimzation MR bias correction functional using CUDA.
class EntropyMinimizationIntensityCorrectionFunctionalCUDA
    /// Inherit non-GPU base class.
  : public EntropyMinimizationIntensityCorrectionFunctionalBase
{
public:
  /// This class type.
  typedef EntropyMinimizationIntensityCorrectionFunctionalCUDA Self;

  /// Pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass type.
  typedef EntropyMinimizationIntensityCorrectionFunctionalBase Superclass;

  /// Constructor.
  EntropyMinimizationIntensityCorrectionFunctionalCUDA( const size_t degreeAdd = 0, const size_t degreeMul = 0 ) : m_PolyDegreeAdd( degreeAdd ), m_PolyDegreeMul( degreeMul ) {}

  /// Virtual destructor.
  virtual ~EntropyMinimizationIntensityCorrectionFunctionalCUDA() {}

  /// Set input image.
  virtual void SetInputImage( UniformVolume::SmartConstPtr& inputImage );

  /// Set additive bias field.
  virtual void SetBiasFieldAdd( const UniformVolume& biasFieldAdd );

  /// Set multiplicative bias field.
  virtual void SetBiasFieldMul( const UniformVolume& biasFieldMul );

  /// Set multiplicative bias field polynomial degree.
  virtual void SetDegreeMul( const size_t degree )
  {
    this->m_PolyDegreeMul = degree;
  }

  /// Set multiplicative bias field polynomial degree.
  virtual void SetDegreeAdd( const size_t degree )
  {
    this->m_PolyDegreeAdd = degree;
  }

  /// Get number of additive monomials.
  virtual size_t GetNumberOfMonomialsAdd() const;

  /// Get number of multiplicative monomials.
  virtual size_t GetNumberOfMonomialsMul() const;

  /// Return parameter vector length.
  virtual size_t ParamVectorDim() const
  {
    return GetNumberOfMonomialsAdd() + GetNumberOfMonomialsMul();
  }
  
protected:
  /// Additive bias field polynomial degree.
  size_t m_PolyDegreeAdd;

  /// Multiplicative bias field polynomial degree.
  size_t m_PolyDegreeMul;

  /// Input image in CUDA memory.
  UniformVolumeCUDA::SmartPtr m_InputImageCUDA;

  /// Number of image pixels.
  size_t m_NumberOfPixels;

  /// Additive bias field in device memory.
  DeviceMemoryCUDA<float>::SmartPtr m_BiasFieldAddCUDA;

  /// Multiplicative bias field.
  DeviceMemoryCUDA<float>::SmartPtr m_BiasFieldMulCUDA;

  /// Output image data.
  DeviceMemoryCUDA<float>::SmartPtr m_OutputDataCUDA;

  /// Update polynomial correctionfactors from input image.
  virtual void UpdateCorrectionFactors();

  /// Jointly update both bias images.
  virtual void UpdateBiasFields( const bool foregroundOnly = true );

  /// Update additive bias image.
  virtual void UpdateBiasFieldAdd( const bool foregroundOnly = true );

  /// Update additive bias image.
  virtual void UpdateBiasFieldMul( const bool foregroundOnly = true );

  /// Update output image.
  virtual void UpdateOutputImage( const bool foregroundOnly = true );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDA_h_included_

