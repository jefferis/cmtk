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

#ifndef __cmtkUniformVolumePainter_h_included_
#define __cmtkUniformVolumePainter_h_included_

#include <cmtkconfig.h>

#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for painting in uniform volume objects.
 * This class provides operations to draw simple geometric objects into UniformVolume objects.
 * This is useful, for example, to create electronic phantom images.
 */
class UniformVolumePainter
{
public:
  /// This class.
  typedef UniformVolumePainter Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Enum for coordinate mode.
  typedef enum
  {
    /// Absolute coordinates: all coordinates are absoluate values within the volume field of view.
    COORDINATES_ABSOLUTE = 0,
    /// Relative coordinates: all coordinates are in range [0,1] which is mapped to volume field of view for each dimension.
    COORDINATES_RELATIVE = 1,
    /// Indexed grid coordinates: all coordinates are grid indexes.
    COORDINATES_INDEXED = 2
  } CoordinateModeEnum;
  

  /// Constructor: link to target volume.
  UniformVolumePainter( UniformVolume::SmartPtr& volume, const CoordinateModeEnum coordinateMode = Self::COORDINATES_INDEXED ) : m_Volume( volume ), m_CoordinateMode( coordinateMode ) {}

  /// Draw a sphere.
  void DrawSphere( const UniformVolume::CoordinateVectorType& center, const Types::Coordinate radius, const Types::DataItem value );

  /// Draw a box.
  void DrawBox( const UniformVolume::CoordinateVectorType& boxFrom, const UniformVolume::CoordinateVectorType& boxTo, const Types::DataItem value );

private:
  /// Pointer to target volume.
  UniformVolume::SmartPtr m_Volume;

  /// Coordinate mode for all operations.
  Self::CoordinateModeEnum m_CoordinateMode;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkUniformVolumePainter_h_included_
