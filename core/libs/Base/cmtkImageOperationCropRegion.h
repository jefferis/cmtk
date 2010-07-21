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

#ifndef __cmtkImageOperationCropRegion_h_included_
#define __cmtkImageOperationCropRegion_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkImageOperation.h"

namespace
cmtk
{

/// Image operation: crop to region.
class ImageOperationCropRegion
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationCropRegion( const DataGrid::RegionType& region ) 
  {
    this->m_Region = region; 
  }
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a new downsampler.
  static void New( const char* arg );
  
private:
  /// Cropping region: x0,y0,z0,x1,y1,z1
  DataGrid::RegionType m_Region;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationCropRegion_h_included_
