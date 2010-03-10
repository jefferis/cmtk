/*
//
//  Copyright 2009,2010 SRI International
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

#ifndef __cmtkImageOperationMedialSkeleton_h_included_
#define __cmtkImageOperationMedialSkeleton_h_included_

#include <cmtkconfig.h>

#include <cmtkImageOperation.h>
#include <cmtkUniformDistanceMap.h>

namespace
cmtk
{

/// Compute medical skeleton of a (binary) mask image.
class ImageOperationMedialSkeleton
/// Inherit generic image operation.
  : public ImageOperation
{
public:
  /// Distance map type.
  typedef UniformDistanceMap<double> DistanceMapType;

  /// Constructor.
  ImageOperationMedialSkeleton( const bool signedDistance ) {}

  /// Apply this operation to an image in place.
  virtual UniformVolume::SmartPtr Apply( UniformVolume::SmartPtr& volume );

  /// Create new medial skeleton operation.
  static void New()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationMedialSkeleton() ) );
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationMedialSkeleton_h_included_
