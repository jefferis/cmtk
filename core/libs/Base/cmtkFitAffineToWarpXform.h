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
#include <Base/cmtkWarpXform.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Fit affine transformation to nonrigid (B-spline or deformation field) transformation.
 */
class FitAffineToWarpXform
{
public:
  /// This class.
  typedef FitAffineToWarpXform Self;

  /// Constructor.
  FitAffineToWarpXform( WarpXform::SmartConstPtr warp );

  /// Fit affine transformation.
  AffineXform::SmartPtr Fit();
  
private:
  /// Input nonrigid warp transformation.
  WarpXform::SmartConstPtr m_WarpXform;

  /// Compute mean translation vector.
  FixedVector<3,Types::Coordinate> GetMeanTranslation( const WarpXform& warpXform /*!< Reference to current warp transformation.*/ ) const;
};

} // namespace

#endif // #ifndef __cmtkFitAffineToWarpXform_h_included_
