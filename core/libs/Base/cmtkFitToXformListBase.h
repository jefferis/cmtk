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

#ifndef __cmtkFitToXformListBase_h_included_
#define __cmtkFitToXformListBase_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkMatrix3x3.h>
#include <Base/cmtkXformList.h>
#include <Base/cmtkImageTemplate.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Fit affine transformation to nonrigid (B-spline or deformation field) transformation.
 */
class FitToXformListBase
{
public:
  /// This class.
  typedef FitToXformListBase Self;

  /// Constructor.
  FitToXformListBase( const UniformVolume& sampleGrid /*!< Discrete pixel grid where the fitted transformation is sampled and residuals are minimized.*/,
		      const XformList& xformList /*!< List of concatenated transformation that the affine transformation is fitted to.*/,
		      const bool absolute = true /*!< Flag fitting absolute transformation vs. relative deformation field */ );
  
protected:
  /// Sampled transformation field.
  ImageTemplate<Xform::SpaceVectorType> m_XformField;
  
  /// Bit flags to mark pixels where the transformation is valid or not (e.g., due to failed numerical inversion).
  std::vector<bool> m_XformValidAt;
};

} // namespace

#endif // #ifndef __cmtkFitToXformListBase_h_included_
