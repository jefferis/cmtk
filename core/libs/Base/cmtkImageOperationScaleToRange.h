/*
//
//  Copyright 2010-2011 SRI International
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

#ifndef __cmtkImageOperationScaleToRange_h_included_
#define __cmtkImageOperationScaleToRange_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

namespace cmtk {

/** \addtogroup Base */
//@{

/// Image operation: scale image values to given range.
class ImageOperationScaleToRange
    /// Inherit from image operation base class.
    : public ImageOperation {
 public:
  /// Constructor.
  ImageOperationScaleToRange(const Types::DataItemRange &toRange)
      : m_ToRange(toRange) {}

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr Apply(
      cmtk::UniformVolume::SmartPtr &volume);

  /// Create a new lower thresholding operation.
  static void New(const char *range);

 private:
  /// Start of range we're scaling to.
  Types::DataItemRange m_ToRange;
};

//@}

}  // namespace cmtk

#endif  // #ifndef __cmtkImageOperationScaleToRange_h_included_
