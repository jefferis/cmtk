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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkImageOperationFlip_h_included_
#define __cmtkImageOperationFlip_h_included_

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: flip.
class ImageOperationFlip
  /// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationFlip( const int normalAxis ) : m_NormalAxis( normalAxis ) {}

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->ApplyMirrorPlane( this->m_NormalAxis );
    return volume;
  }

  /// Create x flip object.
  static void NewX()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationFlip( cmtk::AXIS_X ) ) );
  }

  /// Create y flip object.
  static void NewY()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationFlip( cmtk::AXIS_Y ) ) );
  }

  /// Create y flip object.
  static void NewZ()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationFlip( cmtk::AXIS_Z ) ) );
  }

private:
  /// The normal axis of the flip.
  int m_NormalAxis;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationFlip_h_included_
