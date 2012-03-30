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

#ifndef __cmtkImageOperationHoughTransformSphere_h_included_
#define __cmtkImageOperationHoughTransformSphere_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkUniformVolumeHoughTransformSphere.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: grid downsampling.
class ImageOperationHoughTransformSphere
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationHoughTransformSphere( const Types::Coordinate& radius ) : m_Radius( radius ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->SetData( UniformVolumeHoughTransformSphere( volume ).Get( this->m_Radius ) );
    return volume;
  }
  
  /// Create a new filter based on sigma parameter.
  static void New( const Types::Coordinate radius )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationHoughTransformSphere( radius ) ) );
  }
  
private:
  /// Radius of the spheres detected.
  Types::Coordinate m_Radius;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationHoughTransformSphere_h_included_
