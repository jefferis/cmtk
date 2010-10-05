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

#ifndef __cmtkAffineCongealingFunctional_h_included_
#define __cmtkAffineCongealingFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkCongealingFunctional.h>

#include <Base/cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for affine congealing.
 * This functional evaluates Lilla Zollei's entropy criterion for massively
 * groupwise image registration.
 *
 *\References
 *
 * [1] L . Zöllei, E. Learned-Miller, E. Grimson, W.M. Wells III: "Efficient 
 *     Population Registration of 3D Data", ICCV 2005, Computer Vision for
 *     Biomedical Image Applications; Beijing, China
 */
typedef CongealingFunctional<AffineXform> AffineCongealingFunctional;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkAffineCongealingFunctional_h_included_
