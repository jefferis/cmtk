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

#ifndef __cmtkLeastSquaresPolynomialIntensityBiasField_h_included_
#define __cmtkLeastSquaresPolynomialIntensityBiasField_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkDataGrid.h>

#include <System/cmtkException.h>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Least-squares fit of a polynomial intensity bias field.
 */
class LeastSquaresPolynomialIntensityBiasField
{
public:
  /// This class.
  typedef LeastSquaresPolynomialIntensityBiasField Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart const pointer to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Exception thrown if there are no non-zero mask pixels.
  class EmptyMaskException : public Exception {};

  /// Constructor.
  LeastSquaresPolynomialIntensityBiasField( const UniformVolume& image /*!< Image for which bias field is estimated.*/, 
					    const std::vector<bool>& mask /*!< Mask vector - one bool per image pixel. Only pixels with "true" mask entry will be considered for bias estimation.*/, 
					    const int degree /*!< Polynomial degree of the estimated bias field.*/ );

  /// Get estimated bias field data.
  TypedArray::SmartPtr GetBiasData()
  {
    return this->m_BiasData;
  }

  /// Get estimated bias field data.
  TypedArray::SmartConstPtr GetBiasData() const
  {
    return this->m_BiasData;
  }

  /// Get estimated bias field-corrected data.
  TypedArray::SmartPtr GetCorrectedData()
  {
    return this->m_CorrectedData;
  }

  /// Get estimated bias field-corrected data.
  TypedArray::SmartConstPtr GetCorrectedData() const
  {
    return this->m_CorrectedData;
  }

private:
  /// The estimated bias field data.
  TypedArray::SmartPtr m_BiasData;

  /// The corrected image data.
  TypedArray::SmartPtr m_CorrectedData;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLeastSquaresPolynomialIntensityBiasField_h_included_

