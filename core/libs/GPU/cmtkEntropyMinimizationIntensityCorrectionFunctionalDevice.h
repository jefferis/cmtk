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

#ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_h_included_
#define __cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_h_included_

#include <cmtkconfig.h>

#include <cmtkEntropyMinimizationIntensityCorrectionFunctional.h>

#include "cmtkDeviceMemory.h"
#include "cmtkDeviceUniformVolume.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{
/// Base class for GPU implementation entropy-minimzation MR bias correction functional using Device.
template<unsigned int NOrderAdd,unsigned int NOrderMul>
class EntropyMinimizationIntensityCorrectionFunctionalDevice
    /// Inherit non-GPU base class.
  : public EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
{
public:
  /// This class type.
  typedef EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul> Self;

  /// Pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass type.
  typedef EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul> Superclass;

  /// Virtual destructor.
  virtual ~EntropyMinimizationIntensityCorrectionFunctionalDevice() {}

  /// Set input image.
  virtual void SetInputImage( UniformVolume::SmartConstPtr& inputImage );

protected:
  /// Number of image pixels.
  size_t m_NumberOfPixels;

  /// Input image in device memory.
  DeviceUniformVolume::SmartPtr m_InputImageDevice;

  /// Output image data.
  DeviceMemory<float>::SmartPtr m_OutputDataDevice;

  /// Update output image.
  virtual void UpdateOutputImage( const bool foregroundOnly = true );
};

/// Create functional templated over polynomial degrees.
template<unsigned int NDegreeMul>
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctionalDevice
( const unsigned int polynomialDegreeAdd );

/// Create functional templated over polynomial degrees.
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctionalDevice
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul );

/** Create functional templated over polynomial degrees with initialization from old functional.
 * This function creates a new functional and copies the polynomial coefficients from an existing
 * functional of equal or lower polynomial degrees into the correct locations of the new functional's
 * parameter vector. This is for incremental computation.
 */
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctionalDevice
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul,
  EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr oldFunctional );

//@}

} // namespace cmtk

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice.txx"

#endif // #ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_h_included_

