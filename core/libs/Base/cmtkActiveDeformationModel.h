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

#ifndef __cmtkActiveDeformationModel_h_included_
#define __cmtkActiveDeformationModel_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkActiveShapeModel.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>

#include <list>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/// Active deformation model.
template<class W>
class ActiveDeformationModel :
  /// This is an active shape model of sorts.
  public ActiveShapeModel,
  /// This class is also a deformation (as defined by template parameter).
  public W
{
public:
  /// Smart pointer to spline ADM.
  typedef SmartPointer< ActiveDeformationModel<W> > SmartPtr;

  /** Build deformation model from list of deformations.
   * All deformations in the list must be of the same type, have the same
   * arrangement of control points, and must be defined on the same
   * domain (reference image).
   */
  ActiveDeformationModel( const std::list< SmartPointer<W> > deformationList, const unsigned int numberOfModes, const bool includeScaleInModel = true,
			  const bool includeReferenceInModel = true );

  /// Compose deformation from mean deformation and modes of variation.
  W* Compose( const Types::Coordinate* weights = NULL );

  /** Decompose a deformation into mean and modes of this model.
   *\param input Input deformation.
   *\param weights Weights of the modes that make up the given input vector.
   * This parameter is optional. If not given, no weights will be returned.
   *\return The value of the multivariate Gaussian PDF represented by this 
   * model atr the location of the input vector.
   */
  float Decompose( const W* input, Types::Coordinate *const weights = NULL ) const;


private:
  /** Make ADM sample points from the undeformed c.p.g. of a given deformation.
   * This is for inclusion of the reference individual into the model by
   * putting the identity transformation into the sample set.
   */
  Types::Coordinate* MakeSamplePointsReference( const W* deformation );

  /** Make ADM sample points from a given deformation.
   * This function applies the inverse affine transformation of the given
   * deformation to all of its control points. This brings all deformations
   * into the same reference system, which is a necessary condition for model
   * building. Note that otherwise all control point positions are in the
   * coordinate system of the respective floating image, which varies from
   * deformation to deformation.
   */
  Types::Coordinate* MakeSamplePoints( const W* deformation );

  /// Flag whether to include scale factors in model.
  bool IncludeScaleInModel;

  /// Flag whether to include reference image in model.
  bool IncludeReferenceInModel;
};

/// Spline warp active deformation model.
typedef ActiveDeformationModel<SplineWarpXform> SplineActiveDeformationModel;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkActiveDeformationModel_h_included_
