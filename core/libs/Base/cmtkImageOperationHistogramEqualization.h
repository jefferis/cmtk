/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkImageOperationHistogramEqualization_h_included_
#define __cmtkImageOperationHistogramEqualization_h_included_

#include <cmtkconfig.h>

#include <cmtkImageOperation.h>

namespace
cmtk
{

/// Image operation: histogram equalization with optional number of bins.
class ImageOperationHistogramEqualization
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Default number of bins for histogram equalization.
  static const size_t DefaultNumberOfBins = 1024;

  /// Constructor.
  ImageOperationHistogramEqualization( const size_t nBins ) : m_NumberOfBins( nBins ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a histogram equalization operation with default number of bins.
  static void New();
  
  /// Create a histogram equalization operation with user-supplied number of bins.
  static void NewBins( const long int nBins);
  
private:
  /// Number of histogram bins.
  size_t m_NumberOfBins;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationHistogramEqualization_h_included_
