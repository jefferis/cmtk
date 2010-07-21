/*
//
//  Copyright 2009-2010 SRI International
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

#ifndef __cmtkImageOperationCropThreshold_h_included_
#define __cmtkImageOperationCropThreshold_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkImageOperation.h"

namespace
cmtk
{

/// Image operation: crop by threshold
class ImageOperationCropThreshold
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationCropThreshold( const double threshold, const bool writeRegion = false, const bool writeXform = false ) : m_Threshold( threshold ), m_WriteRegion( writeRegion ), m_WriteXform( writeXform ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a new crop operation.
  static void New( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationCropThreshold( threshold ) ) );
  }
  
  /// Create a new crop operation with region output.
  static void NewWriteRegion( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationCropThreshold( threshold, true, false ) ) );
  }
  
  /// Create a new crop operation with transformation output.
  static void NewWriteXform( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationCropThreshold( threshold, false, true ) ) );
  }
  
private:
  /// Cropping threshold.
  double m_Threshold;

  /// Flag for writing region to standard output.
  bool m_WriteRegion;

  /// Flag for writing transformation to standard output.
  bool m_WriteXform;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationCropThreshold_h_included_
