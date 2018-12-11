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

#ifndef __cmtkLabelCombinationVoting_h_included_
#define __cmtkLabelCombinationVoting_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypedArray.h>
#include <System/cmtkSmartPtr.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Label voting image combination.
 * This class implements combination of multiple multi-class or binary label images
 * using label voting. Each pixel in the output image is assigned the label that the
 * majority of input images assign to that pixel. Pixels with tied voting are assigned
 * a value of 256. The output image is allocated as 16bit short data to accommodate
 * this overflow value.
 *\attention All labels must be between 0 and 255.
 */
class
LabelCombinationVoting
{
public:
  /// Constructor: compute label combination.
  LabelCombinationVoting( const std::vector<TypedArray::SmartPtr>& data );

  /// Get result.
  TypedArray::SmartPtr& GetResult()
  {
    return this->m_Result;
  }

private:
  /// Resulting data array.
  TypedArray::SmartPtr m_Result;
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationVoting_h_included_
