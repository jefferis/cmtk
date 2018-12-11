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

#ifndef __cmtkFitAffineToXformList_h_included_
#define __cmtkFitAffineToXformList_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkFitToXformListBase.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkMatrix3x3.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Fit affine transformation to series of concatenated, possibly numerically inverted, transformations.
 */
class FitAffineToXformList
  : protected FitToXformListBase
{
public:
  /// This class.
  typedef FitAffineToXformList Self;

  /// Base class.
  typedef FitToXformListBase Superclass;

  /// Constructor.
  FitAffineToXformList( const UniformVolume& sampleGrid /*!< Discrete pixel grid where the spline transformation is sampled and residuals are minimized.*/,
			const XformList& xformList /*!< List of concatenated transformation that the spline transformation is fitted to.*/,
			const bool absolute = true /*!< Flag fitting absolute transformation vs. relative deformation field */ ) : Superclass( sampleGrid, xformList, absolute ) {}

  /// Fit affine transformation.
  AffineXform::SmartPtr Fit( const bool fitRigid = false /*!< If this flag is set, a rigid transformation is fitted, otherwise a full affine transformation */ );
  
private:
  /** Compute rotation, scale, and shear matrix by pseudinverse using previously computed centroid translation.
   * We are using simple pseudoinverse rather than procrustes because we do not care whether
   * the result is rigid (det = 1). In fact, if the underlying transformation is not
   * rigid but full affine, then that is exactly what we want the output to be.
   */
  Matrix3x3<Types::Coordinate> GetMatrixAffinePseudoinverse( const cmtk::FixedVector<3,cmtk::Types::Coordinate>& cFrom /*!< Centroid in "from" space previously computed by GetCentroids member function.*/,
							     const cmtk::FixedVector<3,cmtk::Types::Coordinate>& cTo /*!< Centroid in "to" space previously computed by GetCentroids member function.*/ );
  
  /** Compute rotation matrix by SVD using previously computed centroid translation.
   */
  Matrix3x3<Types::Coordinate> GetMatrixRigidSVD( const cmtk::FixedVector<3,cmtk::Types::Coordinate>& cFrom /*!< Centroid in "from" space previously computed by GetCentroids member function.*/,
						  const cmtk::FixedVector<3,cmtk::Types::Coordinate>& cTo /*!< Centroid in "to" space previously computed by GetCentroids member function.*/ );
};

} // namespace

#endif // #ifndef __cmtkFitAffineToXformList_h_included_
