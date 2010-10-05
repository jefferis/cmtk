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

#ifndef __cmtkImageOperationThreshold_h_included_
#define __cmtkImageOperationThreshold_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

namespace
cmtk
{

/// Image operation: thresholding
class ImageOperationThreshold
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationThreshold( const double threshold, const bool above = false, const bool toPadding = false, const bool binarize = false ) 
    : m_Threshold( threshold ), m_Above( above ), m_ToPadding( toPadding ), m_Binarize( binarize ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a new lower thresholding operation.
  static void NewBelow( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationThreshold( threshold, false, false, false ) ) );
  }
  
  /// Create a new upper thresholding operation.
  static void NewAbove( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationThreshold( threshold, true, false, false ) ) );
  }
  
  /// Create a new lower thresholding operation to padding.
  static void NewBelowToPadding( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationThreshold( threshold, false, true, false ) ) );
  }
  
  /// Create a new upper thresholding operation to padding.
  static void NewAboveToPadding( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationThreshold( threshold, true, true, false ) ) );
  }
  
  /// Create a new binarization thresholding operation.
  static void NewBinarize( const double threshold )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationThreshold( threshold, false, false, true ) ) );
  }
  
private:
  /// Threshold.
  double m_Threshold;

  /// Flag for using the threshold as an upper, rather than lower, threshold.
  bool m_Above;

  /// Flag for setting values beyond threshold to padding, rather than to threshold value.
  bool m_ToPadding;

  /// Flag for binarization.
  bool m_Binarize;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationThreshold_h_included_
