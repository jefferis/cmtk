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

#ifndef __cmtkImagePairSymmetricNonrigidRegistrationFunctional_h_included_
#define __cmtkImagePairSymmetricNonrigidRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkFunctional.h"
#include "Base/cmtkSplineWarpXform.h"

#include "Base/cmtkMacros.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Symmtric-consistent elastic registration functional.
class ImagePairSymmetricNonrigidRegistrationFunctional :
  /** Inherit from generic functional. */
  public Functional
{
public:
  /// This class.
  typedef ImagePairSymmetricNonrigidRegistrationFunctional Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef Functional Superclass;

  /// Set inverse consistency weight.
  virtual void SetInverseConsistencyWeight( const Self::ReturnType ) = 0;
  
  /// Set adaptive parameter fixing flag.
  virtual void SetAdaptiveFixParameters( const bool ) = 0;

  /// Set adaptive parameter fixing flag.
  virtual void SetAdaptiveFixThreshFactor( const Self::ReturnType ) = 0;

  /// Set Jacobian constraint weight.
  virtual void SetJacobianConstraintWeight( const Self::ReturnType ) = 0;
  
  /// Set smoothness constraint weight.
  virtual void SetGridEnergyWeight( const Self::ReturnType ) = 0;

  /// Set warp for forward and backward functional.
  virtual void SetWarpXform( SplineWarpXform::SmartPtr& warpFwd, SplineWarpXform::SmartPtr& warpBwd ) = 0;

  /// Constructor function.
  static ImagePairSymmetricNonrigidRegistrationFunctional* Create
  ( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume, const Interpolators::InterpolationEnum interpolation );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSymmetricNonrigidRegistrationFunctional_h_included_
