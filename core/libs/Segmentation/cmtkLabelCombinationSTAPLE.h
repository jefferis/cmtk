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

#ifndef __cmtkLabelCombinationSTAPLE_h_included_
#define __cmtkLabelCombinationSTAPLE_h_included_

#include <cmtkconfig.h>

#include <cmtkTypedArray.h>
#include <cmtkSmartPtr.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Binary STAPLE label combination.
  * This class implements combination of binary images using the
  * binary STAPLE algorithm. If multi-class images are provided as inputs,
  * all values other than zero are interpreted as binary "1" for the
  * purpose of the algorithm.
  */
class
LabelCombinationSTAPLE
{
public:
  /// Constructor: compute label combination.
  LabelCombinationSTAPLE( const std::vector<TypedArray::SmartPtr>& data, //!< Array of typed arrays with input data.
			  const int maxIterations, //!< Maximum number of STAPLE iterations. 
			  const ScalarDataType resultType = TYPE_FLOAT //!< Primitive data type for results.
    );

  /// Get result.
  TypedArray::SmartPtr& GetResult()
  {
    return this->m_Result;
  }

  /// Get one p value.
  float GetPValue( const size_t i ) const
  {
    return this->m_VecP[i];
  }

  /// Get one q value.
  float GetQValue( const size_t i ) const
  {
    return this->m_VecQ[i];
  }

private:
  /// Resulting data array.
  TypedArray::SmartPtr m_Result;

  /// p-Values.
  std::vector<float> m_VecP;

  /// q-Values.
  std::vector<float> m_VecQ;
};

}; // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationSTAPLE_h_included_
