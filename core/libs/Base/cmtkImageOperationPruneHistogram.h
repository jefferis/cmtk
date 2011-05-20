/*
//
//  Copyright 2009-2011 SRI International
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
//  $Revision: 2902 $
//
//  $LastChangedDate: 2011-02-24 12:12:46 -0800 (Thu, 24 Feb 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkImageOperationPruneHistogram_h_included_
#define __cmtkImageOperationPruneHistogram_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: thresholding
class ImageOperationPruneHistogram
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationPruneHistogram( const long int numberOfBins, const bool high = false, const bool low = false ) 
    : m_NumberOfBins( numberOfBins ), m_High( high ), m_Low( low ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create an upper plus lower thresholding operation.
  static void New( const long int nBins )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationPruneHistogram( nBins, true, true ) ) );
  }
  
  /// Create an upper thresholding operation.
  static void NewHigh( const long int nBins )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationPruneHistogram( nBins, true, false ) ) );
  }
  
  /// Create a lower thresholding operation.
  static void NewLow( const long int nBins )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationPruneHistogram( nBins, false, true ) ) );
  }
  
private:
  /// PruneHistogram.
  long int m_NumberOfBins;

  /// Flag for using the threshold as upper threshold.
  bool m_High;

  /// Flag for using the threshold as lower threshold.
  bool m_Low;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationPruneHistogram_h_included_
