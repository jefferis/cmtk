/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkRegistration2d2d_h_included_
#define __cmtkRegistration2d2d_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkMatrix3x3.h"
#include "Base/cmtkScalarImage.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Class for registration of two 2D images.
class Registration2d2d
{
public:
  /// Register two 2D images.
  static void Register( CoordinateMatrix3x3& matrix, ScalarImage::SmartPtr& refImage, ScalarImage::SmartPtr& fltImage );

  /// Register two 2D images, where the second is cropped to a given ROI.
  static void Register( CoordinateMatrix3x3& matrix, ScalarImage::SmartPtr& refImage, ScalarImage::SmartPtr& fltImage, const ScalarImage::RegionType* fltROI );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRegistration2d2d_h_included_
