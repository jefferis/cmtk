/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkImagePairSymmetricAffineRegistrationFunctional_h_included_
#define __cmtkImagePairSymmetricAffineRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkFunctional.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkUniformVolume.h>

#include <Base/cmtkMacros.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Symmtric affine registration functional for simultaneous forward/inverse transformation estimation.
class ImagePairSymmetricAffineRegistrationFunctional :
  /** Inherit from generic functional. */
  public Functional
{
public:
  /// This class.
  typedef ImagePairSymmetricAffineRegistrationFunctional Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef Functional Superclass;

  /// Constructor.
  ImagePairSymmetricAffineRegistrationFunctional( AffineXform::SmartPtr& affineXform ) : m_FwdXform( affineXform ) {};

  /// Constructor function.
  static ImagePairSymmetricAffineRegistrationFunctional* Create
  ( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume, const Interpolators::InterpolationEnum interpolation, AffineXform::SmartPtr& affineXform );

protected:
  /// Forward transformation.
  AffineXform::SmartPtr m_FwdXform;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSymmetricAffineRegistrationFunctional_h_included_
