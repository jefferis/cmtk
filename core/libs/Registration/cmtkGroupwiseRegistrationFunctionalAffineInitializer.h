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

#ifndef __cmtkGroupwiseRegistrationFunctionalAffineInitializer_h_included_
#define __cmtkGroupwiseRegistrationFunctionalAffineInitializer_h_included_

#include <cmtkconfig.h>

#include <cmtkGroupwiseRegistrationFunctionalBase.h>

#include <cmtkSmartPtr.h>

#include <cmtkUniformVolume.h>
#include <cmtkAffineXform.h>
#include <cmtkHistogram.h>

#include <vector>

#include <cmtkClassStream.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Affine initialization of groupwise registration functionals.
 */
class GroupwiseRegistrationFunctionalAffineInitializer
{
public:
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalAffineInitializer Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Initialize affine transformations.
   */
  static void InitializeXforms
  ( GroupwiseRegistrationFunctionalBase& functional, //!<The functional to initialize.
    const bool alignCenters = true, //!< If set, the centers of all target images will be aligned with the center of the template grid via translations.
    const bool alignCenterOfMass = false, ///!< If set, target images will be aligned via translations according to their centers of mass.
    const bool initScales = false //!< If set, approximate initial scaling factors will be computed based on image centers of mass and moments.
    );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkGroupwiseRegistrationFunctionalAffineInitializer_h_included_
