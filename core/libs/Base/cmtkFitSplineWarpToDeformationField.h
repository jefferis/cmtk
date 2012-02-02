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

#ifndef __cmtkFitSplineWarpToDeformationField_h_included_
#define __cmtkFitSplineWarpToDeformationField_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkDeformationField.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkCubicSpline.h>
#include <Base/cmtkRegion.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Fit B-spline-based free-form deformation to pixel-wise deformation field.
 * This class implements the algorithm from: N. J. Tustison, B. B. Avants, and J. C. Gee, "Directly manipulated free-form deformation image registration," IEEE Transactions on Image Processing, vol. 18, no. 3, pp. 624-635, 2009;
 * http://dx.doi.org/10.1109/TIP.2008.2010072
 */
class FitSplineWarpToDeformationField
{
public:
  /// This class.
  typedef FitSplineWarpToDeformationField Self;

  /// Constructor.
  FitSplineWarpToDeformationField( DeformationField::SmartConstPtr dfield, const bool absolute /*!< Flag for absolute transformation vs. relative deformation field */ );

  /// Fit spline warp based on final grid dimensions.
  SplineWarpXform::SmartPtr Fit( const SplineWarpXform::ControlPointIndexType& finalDims, const int nLevels );

  /// Fit spline warp based on final grid spacing.
  SplineWarpXform::SmartPtr Fit( const Types::Coordinate finalSpacing /*!< Final control point spacing of the fitted B-spline free-form deformation*/, 
				 const int nLevels = 1 /*!< Number of levels for optional multi-resolution fit (default: single-resolution fit)*/  );

private:
  /** Flag for absolute vs. relative deformation fields.
   * If this is true, the deformation field is an absolute transformation. 
   * If this is false, the deformation field is a relative offset field.
   */
  bool m_DFieldIsAbsolute;

  /// Input deformation field.
  DeformationField::SmartConstPtr m_DeformationField;

  /// Deformation field residuals, i.e., pixel-wise difference between B-spline transformation and deformation field.
  std::vector< FixedVector<3,Types::Coordinate> > m_Residuals;

  /// Deformation field coverage, i.e., field of fiew.
  Region<3,Types::Coordinate> m_DeformationFieldFOV;

  /// Get deformation grid region corresponding to given coordinate region.
  DataGrid::RegionType GetDeformationGridRange( const UniformVolume::CoordinateRegionType& region ) const;

  /// Compute residuals, i.e., pixel-wise difference between B-spline transformation and deformation field.
  void ComputeResiduals( const SplineWarpXform& splineWarp );

  /// Fit spline warp based on initial warp object.
  void FitSpline( SplineWarpXform& splineWarp, const int nLevels );
};

} // namespace

#endif // #ifndef __cmtkFitSplineWarpToDeformationField_h_included_
