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

#ifndef __cmtkTemplateMultiChannelRegistrationFunctional_h_included_
#define __cmtkTemplateMultiChannelRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkMultiChannelRMIRegistrationFunctional.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for transformation-templated multi-channel registration functional. */
template<class TXform, class TMetricFunctional>
class TemplateMultiChannelRegistrationFunctional :
  /* Inherit from multi-channel registration functional base class. */
  public TMetricFunctional
{
public:
  /** This class. */
  typedef TemplateMultiChannelRegistrationFunctional<TXform,TMetricFunctional> Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** This class. */
  typedef TMetricFunctional Superclass;

  /** The transformation type. */
  typedef TXform TransformationType;

  /** Get transformation. */
  TransformationType& GetTransformation() { return this->m_Transformation; }

  /** Get constant transformation. */
  const TransformationType& GetTransformation() const { return this->m_Transformation; }

  /** Set number of degrees of freedom for transformation. */
  void SetNumberDOFs( const int numberDOFs )
  {
    this->m_Transformation.SetNumberDOFs( numberDOFs );
  }

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v )  
  {
    this->m_Transformation.GetParamVector( v );
  };

  /// Return parameter vector.
  virtual void SetParamVector ( CoordinateVector& v )  
  {
    this->m_Transformation.SetParamVector( v );
  };

  /// Return parameter stepping.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    return this->m_Transformation.GetParamStep( idx, this->m_ReferenceSize, mmStep );
  }

  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const 
  {
    return this->m_Transformation.ParamVectorDim();
  }

  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const 
  {
    return this->m_Transformation.VariableParamVectorDim();
  }

protected:
  /** The templated coordinate transformation. */
  TransformationType m_Transformation;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTemplateMultiChannelRegistrationFunctional_h_included_
