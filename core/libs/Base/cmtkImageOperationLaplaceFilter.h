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

#ifndef __cmtkImageOperationLaplaceFilter_h_included_
#define __cmtkImageOperationLaplaceFilter_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkUniformVolumeLaplaceFilter.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: Laplacian edge enhancement filter.
class ImageOperationLaplaceFilter
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationLaplaceFilter() {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->SetData( UniformVolumeLaplaceFilter( volume ).Get() );
    return volume;
  }
  
  /// Create a new filter.
  static void New()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationLaplaceFilter ) );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationLaplaceFilter_h_included_
