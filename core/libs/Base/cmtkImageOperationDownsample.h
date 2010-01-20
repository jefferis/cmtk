/*
//
//  Copyright 2009 SRI International
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

#ifndef __cmtkImageOperationDownsample_h_included_
#define __cmtkImageOperationDownsample_h_included_

#include <cmtkconfig.h>

#include <cmtkImageOperation.h>

namespace
cmtk
{

/// Image operation: grid downsampling.
class ImageOperationDownsample
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationDownsample( const int factorX, const int factorY, const int factorZ ) : m_FactorX( factorX ), m_FactorY( factorY ), m_FactorZ( factorZ ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a new downsampler.
  static void New( const char* arg );
  
private:
  /// Downsampling factor in X direction.
  int m_FactorX;

  /// Downsampling factor in Y direction.
  int m_FactorY;

  /// Downsampling factor in Z direction.
  int m_FactorZ;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationDownsample_h_included_
