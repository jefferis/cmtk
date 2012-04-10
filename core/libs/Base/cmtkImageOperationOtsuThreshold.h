/*
//
//  Copyright 2009-2012 SRI International
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

#ifndef __cmtkImageOperationOtsuThreshold_h_included_
#define __cmtkImageOperationOtsuThreshold_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: Otsu threshold binarization
class ImageOperationOtsuThreshold
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationOtsuThreshold( const size_t nBins = 1024 /*!< Number of histogram bins for threshold computation.*/ ) : m_NumberOfBins( nBins ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a new thresholding operation.
  static void New()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationOtsuThreshold() ) );
  }
  
  /// Create a new thresholding operation with explicit number of histogram bins.
  static void NewBins( const long int nBins = 1024 /*!< Number of histogram bins for threshold computation.*/ )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationOtsuThreshold( nBins ) ) );
  }
  
private:
  /// Number of histogram bins for threshold computation.
  size_t m_NumberOfBins;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationOtsuThreshold_h_included_
