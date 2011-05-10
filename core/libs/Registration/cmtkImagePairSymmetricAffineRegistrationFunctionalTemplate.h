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

#ifndef __cmtkImagePairSymmetricAffineRegistrationFunctionalTemplate_h_included_
#define __cmtkImagePairSymmetricAffineRegistrationFunctionalTemplate_h_included_

#include <cmtkconfig.h>

#include "cmtkImagePairSymmetricAffineRegistrationFunctionalTemplate.h"

#include <Base/cmtkFunctional.h>
#include <Base/cmtkMacros.h>

#include <Registration/cmtkImagePairAffineRegistrationFunctionalTemplate.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Template for symmtric affine registration functional for simultaneous forward/inverse transformation estimation.
template<class VM>
class ImagePairSymmetricAffineRegistrationFunctionalTemplate :
  /** Inherit from non-template base functional class. */
  public ImagePairSymmetricAffineRegistrationFunctional
{
public:
  /// This class.
  typedef ImagePairSymmetricAffineRegistrationFunctionalTemplate<VM> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef ImagePairSymmetricAffineRegistrationFunctional Superclass;

  /// The forward functional.
  ImagePairAffineRegistrationFunctionalTemplate<VM> FwdFunctional;

  /// The backward functional.
  ImagePairAffineRegistrationFunctionalTemplate<VM> BwdFunctional;

  /// Constructor.
  ImagePairSymmetricAffineRegistrationFunctionalTemplate( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating, const Interpolators::InterpolationEnum interpolation, AffineXform::SmartPtr& affineXform )
    : Superclass( affineXform ),
      FwdFunctional( reference, floating, interpolation, affineXform ),
      BwdFunctional( floating, reference, interpolation, affineXform->GetInverse() )
  {}

  /// Set flag and value for forcing values outside the floating image.
  virtual void SetForceOutside
  ( const bool flag = true, const Types::DataItem value = 0 )
  {
    this->FwdFunctional.SetForceOutside( flag, value );
    this->BwdFunctional.SetForceOutside( flag, value );
  }

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v )
  {
    this->FwdFunctional.GetParamVector( v );
  }

  /// Evaluate functional value and gradient.
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 );
    
  /// Evaluate functional value.
  virtual typename Self::ReturnType EvaluateAt ( CoordinateVector& v ) 
  {
    this->m_FwdXform->SetParamVector( v );

    CoordinateVector vInv;
    this->m_FwdXform->GetInverse()->GetParamVector( vInv );
    
    return this->FwdFunctional.EvaluateAt( v ) + this->BwdFunctional.EvaluateAt( vInv );
  }

  /// Evaluate functional with current parameter vector.
  virtual typename Self::ReturnType Evaluate () 
  {
    return this->FwdFunctional.Evaluate() + this->BwdFunctional.Evaluate();
  }
  
  /// Get parameter stepping in milimeters.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    return 0.5 * (this->FwdFunctional.GetParamStep( idx, mmStep ) + this->FwdFunctional.GetParamStep( idx, mmStep ) );
  }
  
  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const
  {
    return this->FwdFunctional.ParamVectorDim();
  }
  
  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const 
  {
    return this->FwdFunctional.VariableParamVectorDim();
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSymmetricAffineRegistrationFunctionalTemplate_h_included_
