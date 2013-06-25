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

#ifndef __cmtkFitSplineWarpToLandmarks_h_included_
#define __cmtkFitSplineWarpToLandmarks_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkLandmarkPairList.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkCubicSpline.h>
#include <Base/cmtkRegion.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Fit B-spline-based free-form deformation to a list of matched landmarks.
 *\see This class implements the algorithm from: N. J. Tustison, B. B. Avants, and J. C. Gee, "Directly manipulated free-form deformation image registration," IEEE Transactions on Image Processing, vol. 18, no. 3, pp. 624-635, 2009;
 * http://dx.doi.org/10.1109/TIP.2008.2010072
 *\see The implementation itself is more closely following S. Lee, G. Wolberg, and S. Y. Shin, “Scattered data interpolation with multilevel B-splines,” IEEE Transactions on Visualization and Computer Graphics, 
 * vol. 3, no. 3, pp. 228-244, 1997. http://dx.doi.org/10.1109/2945.620490
 *
 *\attention This class does not compile with Solaris Studio, unless STLport support is activated via the "-library=stlport4" compiler command line option.
 */
class FitSplineWarpToLandmarks
{
public:
  /// This class.
  typedef FitSplineWarpToLandmarks Self;

  /// Class for parameters to the fitting algorithm.
  class Parameters
  {
  public:
    /// Default constructor. 
    Parameters() : m_Levels( 1 ), m_IterationsPerLevel( 100 ), m_ResidualThreshold( 0 ) {}

    /// Number of levels in the multi-resolution fitting.
    int m_Levels;

    /// Number of update iterations per level in the multi-resolution fitting.
    int m_IterationsPerLevel;

    /// Threshold for relative RMS residual improvement. Iteration terminates if (rmsAfterUpdate-rmsBeforeUpdate)/rmsBeforeUpdate < threshold.
    Types::Coordinate m_ResidualThreshold;
  };

  /// Constructor.
  FitSplineWarpToLandmarks( const LandmarkPairList& landmarkList );

  /// Fit spline warp based on final grid dimensions.
  SplineWarpXform::SmartPtr Fit( const SplineWarpXform::SpaceVectorType& domain /*!< Domain of the deformation field. This should be the size of the fixed image grid to be used with the resulting deformation */,
  				 const SplineWarpXform::ControlPointIndexType& finalDims /*!< Final spline control point grid dimensions.*/, 
				 const AffineXform* initialAffine = NULL /*!< Optional affine transformation to initialize the spline control points.*/,
				 const Self::Parameters& parameters = Self::DefaultParameters /*!< Fitting parameters.*/ );

  /// Fit spline warp based on final grid spacing.
  SplineWarpXform::SmartPtr Fit( const SplineWarpXform::SpaceVectorType& domain /*!< Domain of the deformation field. This should be the size of the fixed image grid to be used with the resulting deformation */,
				 const Types::Coordinate finalSpacing /*!< Final control point spacing of the fitted B-spline free-form deformation*/, 
				 const AffineXform* initialAffine = NULL /*!< Optional affine transformation to initialize the spline control points.*/,
				 const Self::Parameters& parameters = Self::DefaultParameters /*!< Fitting parameters.*/ );
  
private:
  /// Default parameters.
  static Self::Parameters DefaultParameters;

  /// Input landmarks.
  std::vector<LandmarkPair> m_Landmarks;

  /// Spline grid index per landmark.
  std::vector< FixedVector<3,int> > m_LandmarksGrid;

  /// Spline coeffiecints per landmark.
  std::vector< FixedArray<3, FixedVector<4,Types::Coordinate> > > m_LandmarksSpline;

  /// Deformation field residuals, i.e., pixel-wise difference between B-spline transformation and deformation field.
  std::vector< SplineWarpXform::SpaceVectorType > m_Residuals;

  /** Compute residuals, i.e., pixel-wise difference between B-spline transformation and deformation field.
   *\return Root-of-mean-squared residual over all landmarks.
   */
  Types::Coordinate ComputeResiduals( const SplineWarpXform& splineWarp );

  /// Fit spline warp based on initial warp object.
  void FitSpline( SplineWarpXform& splineWarp /*!< Fitted spline warp is returned here.*/, const Self::Parameters& parameters = Self::DefaultParameters /*!< Fitting parameters.*/ );
};

} // namespace

#endif // #ifndef __cmtkFitSplineWarpToLandmarks_h_included_
