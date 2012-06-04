/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkTransformedVolumeAxes_h_included_
#define __cmtkTransformedVolumeAxes_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkVector3D.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkParametricPlane.h>
#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class that represents pre-transformed axes samples for 3D volumes.
 * This class generates an array of 3D coordinate vectors for each of the
 * three coordinate axes X, Y, and Z. These arrays contain direction vectors
 * that directly point to the coordinate of every sample on the respective
 * axis. These vectors are computed with respect to the object's coordinate
 * transformation as well as a second Volume's transformation and coordinate
 * offset. The Z-axis array contains this Volume's coordinate offset as 
 * defined by the transformation translation component as well.
 * The vectors from the arrays can therefore directly be used for probing
 * the other Volume using the ProbeNoXform member function.
 */
class TransformedVolumeAxes
{
public:
  /// This class.
  typedef TransformedVolumeAxes Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /** Constructor using affine transformation.
   *\param volume The volume for which we are constructing an axes hash.
   *\param xform Coordinate transformation. If this pointer is NULL, identity is assumed.
   *\param deltas If this parameter is given, it is used as a pointer to a
   * three-element coordinate array defining the voxel size in the other volume
   * to obtain fractional voxel indices rather than actual coordinates.
   *
   * Alternatively, this parameter can also be used to provide the volume size (Volume::Size), which creates
   * normalized coordinates (0..1) for each volume axis.
   *
   *\param otherOrigin If this parameter is given, it is used as a pointer to a
   * three-element coordinate array defining an m_Origin vector for the transformed
   * coordinates.
   *\return Pointer to an array of three pointers. Each of these points to an
   * array of Vector3D objects that contain the vector hashes for the X-, Y-, 
   * and Z-axis.
   */
  TransformedVolumeAxes( const UniformVolume& volume, const AffineXform* xform = NULL, const Types::Coordinate* deltas = NULL, const Types::Coordinate* otherOrigin = NULL );
  
  /** Constructor using mirror plane.
   *\param volume The volume whose axes we are transforming.
   *\param mirrorPlane Mirror plane with respect to which the coordinates of
   * this volume and thus all hash values are mirrored.
   *\param deltas If this parameter is given, it is used as a pointer to a
   * 3 element coordinate array defining the voxel size in the other volume
   * to obtain fractional voxel indices rather than actual coordinates. 
   *
   * Alternatively, this parameter can also be used to provide the volume size (Volume::Size), which creates
   * normalized coordinates (0..1) for each volume axis.
   */
  TransformedVolumeAxes( const UniformVolume& volume, const ParametricPlane& mirrorPlane, const Types::Coordinate* deltas = NULL );

  /// Free all storage.
  ~TransformedVolumeAxes();

  /// Access operator.
  const Vector3D* operator[]( const size_t index ) const
  {
    return this->m_Hash[index];
  }

  /// Get dimensions.
  const FixedVector<3,int>& Dims() const
  {
    return this->m_Dims;
  }

private:
  /// Array of pointers to transformed axes points.
  FixedArray<3,UniformVolume::SpaceVectorType*> m_Hash;

  /// Dimensions of the transformed grid: numbers of samples per axis.
  FixedVector<3,int> m_Dims;

  /// Create the actual hash: allocate and fill according to given offset and delta vectors.
  void MakeHash( const UniformVolume& volume, const UniformVolume::SpaceVectorType& offset, const UniformVolume::SpaceVectorType& dX, const UniformVolume::SpaceVectorType& dY, const UniformVolume::SpaceVectorType& dZ );

};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTransformedVolumeAxes_h_included_
