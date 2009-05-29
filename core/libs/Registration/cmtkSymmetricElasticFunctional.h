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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkSymmetricElasticFunctional_h_included_
#define __cmtkSymmetricElasticFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkFunctional.h>

#ifdef CMTK_BUILD_SMP
#  include <cmtkParallelElasticFunctional.h>
#else
#  include <cmtkVoxelMatchingElasticFunctional.h>
#endif

#include <cmtkMacros.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Symmtric-consistent elastic registration functional.
class SymmetricElasticFunctional :
  /** Inherit from generic functional. */
  public Functional
{
public:
  /// This class.
  typedef SymmetricElasticFunctional Self;

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
  
  /// Set Jacobian constraint weight.
  virtual void SetRigidityConstraintWeight( const Self::ReturnType ) = 0;
  
  /// Set smoothness constraint weight.
  virtual void SetGridEnergyWeight( const Self::ReturnType ) = 0;

  /// Set warp for forward and backward functional.
  virtual void SetWarpXform( WarpXform::SmartPtr& warpFwd, WarpXform::SmartPtr& warpBwd ) = 0;
};

/// Template for symmtric-consistent elastic registration functional.
template<class VM, class W>
class SymmetricElasticFunctional_Template :
  /** Inherit from non-template base functional class. */
  public SymmetricElasticFunctional
{
public:
  /// This class.
  typedef SymmetricElasticFunctional_Template<VM,W> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef SymmetricElasticFunctional Superclass;

  /// The forward functional.
#ifdef CMTK_BUILD_SMP
  ParallelElasticFunctional<VM,W> FwdFunctional;
#else
  VoxelMatchingElasticFunctional_Template<VM,W> FwdFunctional;
#endif

  /// The backward functional.
#ifdef CMTK_BUILD_SMP
  ParallelElasticFunctional<VM,W> BwdFunctional;
#else
  VoxelMatchingElasticFunctional_Template<VM,W> BwdFunctional;
#endif

  /// Constructor.
  SymmetricElasticFunctional_Template( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
    : FwdFunctional( reference, floating ),
      BwdFunctional( floating, reference )
  {}

  /// Set inverse consistency weight.
  virtual void SetInverseConsistencyWeight( const typename Self::ReturnType inverseConsistencyWeight ) 
  {
    this->FwdFunctional.SetInverseConsistencyWeight( inverseConsistencyWeight );
    this->BwdFunctional.SetInverseConsistencyWeight( inverseConsistencyWeight );
  }
  
  /// Set adaptive parameter fixing flag.
  virtual void SetAdaptiveFixParameters( const bool adaptiveFixParameters ) 
  {
    this->FwdFunctional.SetAdaptiveFixParameters( adaptiveFixParameters );
    this->BwdFunctional.SetAdaptiveFixParameters( adaptiveFixParameters );
  }

  /// Set adaptive parameter fixing threshold.
  virtual void SetAdaptiveFixThreshFactor( const typename Self::ReturnType threshFactor ) 
  {
    this->FwdFunctional.SetAdaptiveFixThreshFactor( threshFactor );
    this->BwdFunctional.SetAdaptiveFixThreshFactor( threshFactor );
  }

  /// Set Jacobian constraint weight.
  virtual void SetJacobianConstraintWeight( const typename Self::ReturnType jacobianConstraintWeight ) 
  {
    this->FwdFunctional.SetJacobianConstraintWeight( jacobianConstraintWeight );
    this->BwdFunctional.SetJacobianConstraintWeight( jacobianConstraintWeight );
  }
  
  /// Set rigidity constraint weight.
  virtual void SetRigidityConstraintWeight( const typename Self::ReturnType rigidityConstraintWeight ) 
  {
    this->FwdFunctional.SetRigidityConstraintWeight( rigidityConstraintWeight );
    this->BwdFunctional.SetRigidityConstraintWeight( rigidityConstraintWeight );
  }
  
  /// Set smoothness constraint weight.
  virtual void SetGridEnergyWeight( const typename Self::ReturnType gridEnergyWeight ) 
  {
    this->FwdFunctional.SetGridEnergyWeight( gridEnergyWeight );
    this->BwdFunctional.SetGridEnergyWeight( gridEnergyWeight );
  }
  
  /// Set warp for forward and backward functional.
  virtual void SetWarpXform( WarpXform::SmartPtr& warpFwd, WarpXform::SmartPtr& warpBwd );

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v )  
  {
    CoordinateVector vFwd, vBwd;
    this->FwdFunctional.GetParamVector( vFwd );
    this->BwdFunctional.GetParamVector( vBwd );

    v.SetDim( vFwd.Dim + vBwd.Dim );
    v.CopyToOffset( vFwd );
    v.CopyToOffset( vBwd, vFwd.Dim );
  }

  /// Evaluate functional value and gradient.
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 );
    
  /// Evaluate functional value.
  virtual typename Self::ReturnType EvaluateAt ( CoordinateVector& v ) 
  {
    CoordinateVector vFwd( this->FwdFunctional.ParamVectorDim(), v.Elements, false /*freeElements*/ );
    CoordinateVector vBwd( this->BwdFunctional.ParamVectorDim(), v.Elements+this->FwdFunctional.ParamVectorDim(), false /*freeElements*/ );
    return this->FwdFunctional.EvaluateAt( vFwd ) + this->BwdFunctional.EvaluateAt( vBwd );
  }
  
  virtual typename Self::ReturnType Evaluate () 
  {
    return this->FwdFunctional.Evaluate() + this->BwdFunctional.Evaluate();
  }
  
  /// Get parameter stepping in milimeters.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    if ( idx < this->FwdFunctional.ParamVectorDim() )
      return this->FwdFunctional.GetParamStep( idx, mmStep );
    else
      return this->BwdFunctional.GetParamStep( idx - this->FwdFunctional.ParamVectorDim(), mmStep );
  }
  
  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const
  {
    return this->FwdFunctional.ParamVectorDim() + this->BwdFunctional.ParamVectorDim();
  }
  
  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const 
  {
    return this->FwdFunctional.VariableParamVectorDim() + this->BwdFunctional.VariableParamVectorDim();
  }
};

/// Constructor function.
SymmetricElasticFunctional*
CreateSymmetricElasticFunctional( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSymmetricElasticFunctional_h_included_
