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

#ifndef __cmtkSplineWarpXformUniformVolume_h_included_
#define __cmtkSplineWarpXformUniformVolume_h_included_

#include <cmtkconfig.h>

#include <cmtkXformUniformVolume.h>

#include <cmtkUniformVolume.h>
#include <cmtkSplineWarpXform.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Pre-compute transformation for grid locations in a uniform volume.
 */
class SplineWarpXformUniformVolume :
  /// Inherit from class to prevent copying.
  public XformUniformVolume
{
public:
  /// This class.
  typedef SplineWarpXformUniformVolume Self;

  /// Parent class.
  typedef XformUniformVolume Superclass;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  SplineWarpXformUniformVolume( const UniformVolume& volume, const SplineWarpXform::SmartConstPtr& xform );
  
  /// Virtual destructor.
  virtual ~SplineWarpXformUniformVolume() {}
  
  /** Get transformed location of linked grid pixel.
   */
  virtual void GetTransformedGrid( Vector3D& v, const int idxX, const int idxY, const int idxZ ) const;
  
  /** Get transformed locations of a series (scanline) of linked grid pixels.
   */
  virtual void GetTransformedGridSequence( Vector3D *const v, const size_t numPoints, const int idxX, const int idxY, const int idxZ ) const;
  
private:
  /// The linked transformation.
  const SplineWarpXform::SmartConstPtr m_Xform;

  /// Register axes points of the volume to be deformed.
  void RegisterVolume( const UniformVolume& volume );

  /// Register a single axis of the uniform volume to be deformed.
  void RegisterVolumeAxis ( const int, const Types::Coordinate delta, const Types::Coordinate origin, const int, const Types::Coordinate, std::vector<int>& g, 
			    std::vector<Types::Coordinate>& spline, std::vector<Types::Coordinate>& dspline );
  
  /**@name Precomputed grid indices.
   * These arrays hold the precomputed grid indices of the deformed grid's
   * voxels with respect to the control point grid of this deformation.
   */
  //@{
  /// x-axis.
  std::vector<int> gX;
  /// y-axis.
  std::vector<int> gY;
  /// z-axis.
  std::vector<int> gZ;
  //@}

  /**@name Precomputed spline coefficients.
   * These arrays hold the precomputed spline coefficients for deforming the
   * voxel locations in the associated deformed grid.
   */
  //@{
  /// x-axis.
  std::vector<Types::Coordinate> splineX;
  /// y-axis.
  std::vector<Types::Coordinate> splineY;
  /// z-axis.
  std::vector<Types::Coordinate> splineZ;
  //@}

  /**@name Precomputed derivative spline coefficients.
   * These arrays hold the precomputed derivatives of the spline coefficients.
   * This allows for rapid evaluation of the Jacobian determinant.
   */
  //@{
  /// x-axis.
  std::vector<Types::Coordinate> dsplineX;
  /// y-axis.
  std::vector<Types::Coordinate> dsplineY;
  /// z-axis.
  std::vector<Types::Coordinate> dsplineZ;
  //@}

  /// Relative offsets of all control points in a 4 x 4 x 4 neighborhood.
  int GridPointOffset[48];
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkSplineWarpXformUniformVolume_h_included_
