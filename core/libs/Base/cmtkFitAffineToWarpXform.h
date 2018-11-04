/*
//
//  Copyright 2012 SRI International
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

#ifndef __cmtkFitAffineToWarpXform_h_included_
#define __cmtkFitAffineToWarpXform_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkMatrix3x3.h>
#include <Base/cmtkWarpXform.h>

namespace cmtk {

/** \addtogroup Base */
//@{

/** Fit affine transformation to nonrigid (B-spline or deformation field)
 * transformation.
 */
class FitAffineToWarpXform {
 public:
  /// This class.
  typedef FitAffineToWarpXform Self;

  /// Constructor.
  FitAffineToWarpXform(WarpXform::SmartConstPtr warp);

  /// Fit affine transformation.
  AffineXform::SmartPtr Fit();

 private:
  /// Input nonrigid warp transformation.
  WarpXform::SmartConstPtr m_WarpXform;

  /** Compute rotation, scale, and shear matrix using previously computed
   * translation. We are using simple pseudoinverse rather than procrustes
   * because we do not care whether the result is rigid (det = 1). In fact, if
   * the underlying transformation is not rigid but full affine, then that is
   * exactly what we want the output to be.
   */
  static Matrix3x3<Types::Coordinate> GetMatrix(const WarpXform &
                                                    warpXform /*!< Reference to
                                                                 current warp
                                                                 transformation.*/
                                                ,
                                                const cmtk::FixedVector<
                                                    3, cmtk::Types::Coordinate>
                                                    &cFrom /*!< Centroid in
                                                              "from" space
                                                              previously
                                                              computed by
                                                              GetCentroids
                                                              member function.*/
                                                ,
                                                const cmtk::
                                                    FixedVector<
                                                        3,
                                                        cmtk::Types::
                                                            Coordinate> &cTo /*!< Centroid in "to" space previously computed by GetCentroids member function.*/);
};

}  // namespace cmtk

#endif  // #ifndef __cmtkFitAffineToWarpXform_h_included_
