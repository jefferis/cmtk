/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkSplineWarpXform.h>
#include <cmtkUniformVolume.h>

#include <vector>

namespace
cmtk
{

/// Class for largely pre-computed application of spline warp transformation to a uniform volume.
class SplineWarpXformUniformVolume
{
public:
  /// Constructor.
  SplineWarpXformUniformVolume( const SplineWarpXform::SmartPtr& warp, const UniformVolume::SmartPtr& volume );

  /// Get a grid point from the deformed grid.
  void GetTransformedGrid( Vector3D& v, const int idxX, const int idxY, const int idxZ ) const;
  
  /// Get a sequence of grid points from the deformed grid. 
  void GetTransformedGridSequence( Vector3D *const v, const int numPoints, const int idxX, const int idxY, const int idxZ ) const;
  
private:
  /// The spline warp transformation object.
  SplineWarpXform::SmartPtr m_Warp;

  /// The uniform volume object.
  UniformVolume::SmartPtr m_Volume;

  /// Register one axis of the volume grid with this object.
  void RegisterVolumeAxis( const int dim, const Types::Coordinate delta, const Types::Coordinate origin, const int cpgDim, const Types::Coordinate invCpgSpacing,
			   std::vector<int>& g, std::vector<Types::Coordinate>& spline, std::vector<Types::Coordinate>& dspline );
  
  /// Dimensions of the volume image linked to this transformation.
  int m_VolumeDims[3];

  /**@name Precomputed grid indices.
   * These arrays hold the precomputed grid indices of the deformed grid's
   * voxels with respect to the control point grid of this deformation.
   */
  //@{
  /// x-axis.
  std::vector<int> m_GridX;
  /// y-axis.
  std::vector<int> m_GridY;
  /// z-axis.
  std::vector<int> m_GridZ;
  //@}

  /**@name Precomputed spline coefficients.
   * These arrays hold the precomputed spline coefficients for deforming the
   * voxel locations in the associated deformed grid.
   */
  //@{
  /// x-axis.
  std::vector<Types::Coordinate> m_SplineX;
  /// y-axis.
  std::vector<Types::Coordinate> m_SplineY;
  /// z-axis.
  std::vector<Types::Coordinate> m_SplineZ;
  //@}

  /**@name Precomputed derivative spline coefficients.
   * These arrays hold the precomputed derivatives of the spline coefficients.
   * This allows for rapid evaluation of the Jacobian determinant.
   */
  //@{
  /// x-axis.
  std::vector<Types::Coordinate> m_DerivSplineX;
  /// y-axis.
  std::vector<Types::Coordinate> m_DerivSplineY;
  /// z-axis.
  std::vector<Types::Coordinate> m_DerivSplineZ;
  //@}
};

} // namespace cmtk

#endif  // #ifndef __cmtkSplineWarpXformUniformVolume_h_included_
