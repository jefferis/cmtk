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

#ifndef __cmtkLabelCombinationMultiClassSTAPLE_h_included_
#define __cmtkLabelCombinationMultiClassSTAPLE_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypedArray.h>
#include <System/cmtkSmartPtr.h>
#include <Base/cmtkMatrix.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Multi-class STAPLE label combination.
  * This class implements combination of multiple-label images using the
  * multi-class STAPLE algorithm.
  */
class
LabelCombinationMultiClassSTAPLE
{
public:
  /// This class.
  typedef LabelCombinationMultiClassSTAPLE Self;

  /// Real value type for internal computations.
  typedef double RealValueType;

  /// Confusion matrix type.
  typedef Matrix2D<RealValueType> ConfusionMatrixType;
  
  /// Constructor: compute label combination.
  LabelCombinationMultiClassSTAPLE( const std::vector<TypedArray::SmartPtr>& data /*!< Array of typed arrays with input data.*/,
				    const int maxIterations /*!< Maximum number of STAPLE iterations.*/ );

  /// Get result.
  TypedArray::SmartPtr& GetResult()
  {
    return this->m_Result;
  }

private:
  /// Resulting data array.
  TypedArray::SmartPtr m_Result;

  /// Array of prior probabilities per class.
  std::vector<RealValueType> m_Priors;

  /// Array of confusion matrices.
  std::vector<Self::ConfusionMatrixType> m_Confusion;

  /// Array of updated confusion matrices.
  std::vector<Self::ConfusionMatrixType> m_ConfusionNew;
};

}; // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationMultiClassSTAPLE_h_included_
