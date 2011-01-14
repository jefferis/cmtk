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

#ifndef __cmtkFunctional_h_included_
#define __cmtkFunctional_h_included_

#include <cmtkconfig.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkConsole.h>

#include <Base/cmtkVector.h>
#include <Base/cmtkTypes.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class defining a real-valued functional on a multi-dimensional domain.
class Functional
{
public:
  /// This class.
  typedef Functional Self;

  /// Smart pointer to Functional.
  typedef SmartPointer<Self> SmartPtr;

  /** Functional return type.
   * We set this to float (rather than double) because with double we get
   * different results using different CPUs, even between Pentium4 Nocoma and
   * Core2/Xeon using the exact same compiled bainary. That is not good for regression
   * testing, nor is it for consistent results.
   */
  typedef Types::Combined<Types::Coordinate,Types::DataItem>::Type ReturnType;

  /// Parameter vector type.
  typedef Vector<Types::Coordinate> ParameterVectorType;

  /// Functional return type.
  typedef Types::Coordinate ParameterType;

  /// Set parameter vector.
  virtual void SetParamVector ( ParameterVectorType& )
  {
    StdErr << "ERROR: Functional::SetParamVector() was called but not implemented\n";
    exit( 1 );
  }

  /// Return parameter vector.
  virtual void GetParamVector ( ParameterVectorType& )
  {
    StdErr << "ERROR: Functional::GetParamVector() was called but not implemented\n";
    exit( 1 );
  }

  /// Evaluate functional.
  virtual Self::ReturnType Evaluate() { return 0; } // cannot make this abstract because we need to instantiate this for generic initialization

  /// Evaluate functional with new parameter vector.
  virtual Self::ReturnType EvaluateAt( ParameterVectorType& v )
  {
    this->SetParamVector( v );
    return this->Evaluate();
  }

#ifdef CMTK_BUILD_DEMO
  /// Create a snapshot (to disk) of current functional result.
  virtual void SnapshotAt( ParameterVectorType& ) {}
#endif

  /** Evaluate functional with new parameter vector along previously computed gradient direction.
   * By default this function simply calls Evaluate(), but derived classes can override it to
   * provide more computationally efficient, i.e., restricted, implementations.
   */
  virtual Self::ReturnType EvaluateAlongGradientAt( ParameterVectorType& v ) { return this->EvaluateAt( v ); }

  /// Evaluate functional and also return its gradient.
  virtual Self::ReturnType EvaluateWithGradient( ParameterVectorType& v, ParameterVectorType& g, const Types::Coordinate step = 1 );

  /// Return dimension of the parameter vector.
  virtual size_t ParamVectorDim() const = 0;

  /** Return dimension of the parmater vector's variable part.
   * For example, the rotation center of a rigid body transformation is not
   * considered variable. Therefore it should not be used for gradient
   * computation. By default this function returns the total parameter
   * vector length.
   */
  virtual size_t VariableParamVectorDim() const
  {
    return this->ParamVectorDim();
  }

  /// Virtual destructor.
  virtual ~Functional() {};

  /** Get stepping for one parameter.
   * This function serves to make sure that optimizers optimize all parameters with
   * approximately equal impact. If a parameter has smaller impact than another, it should
   * return a smaller value here.
   */
  virtual Types::Coordinate GetParamStep( const size_t, const Types::Coordinate mmStep = 1 ) const 
  {
    return mmStep;
  }

  /** Wiggle a little.
   * If the functional is not time-invariant, e.g., due to randomization,
   * then the optimizer can tell it to re-organize itself for example for
   * a repeated search. 
   *\return If this function returns true, then the functional did indeed
   * wiggle a little, i.e., potentially change. If the function returns 
   * false, then the functional does not support this operation, and the
   * optimizer can assume that further evaluations will produce exactly
   * the same values as before.
   */
  virtual bool Wiggle() { return false; }
};

//@}

} // namespace cmtk

#endif // #ifndef _FUNCTIONAL_H_
