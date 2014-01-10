/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012, 2014 SRI International
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

#ifndef __cmtkPolynomialXform_h_included_
#define __cmtkPolynomialXform_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkXform.h>
#include <Base/cmtkPolynomial.h>
#include <Base/cmtkMatrix3x3.h>

#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** 3D polynomial coordinate transformation.
 */
class PolynomialXform :
  /// Inherit from meta data information container.
  public Xform
{
public:
  /// This class.
  typedef PolynomialXform Self;

  /// This class.
  typedef Xform Superclass;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const for this class.
  typedef SmartConstPointer<Self> SmartConstPtr;
  
  /// Copy constructor.
  PolynomialXform( const PolynomialXform& other )
    : Xform( other ),
      m_Degree( other.m_Degree ),
      m_NumberOfMonomials( other.m_NumberOfMonomials )
  {
  }
  
  /// Default constructor.
  PolynomialXform( const byte degree = 0 /*!< Polynomial degree - 0 through 4 are supported. */ ) : m_Degree( degree )
  {
    // Sort out how many monomials a polynomial of the given degree has.
    this->m_NumberOfMonomials = PolynomialHelper::GetNumberOfMonomials( this->m_Degree );

    // Allocate one set per spatial dimension
    this->AllocateParameterVector( 3 * this->m_NumberOfMonomials );
  }

  /// Clone and return smart pointer.
  Self::SmartPtr Clone () const 
  {
    return Self::SmartPtr( this->CloneVirtual() );
  }

  /// Virtual destructor.
  virtual ~PolynomialXform() {}

  /// Get degree of the polynomial.
  byte Degree() const 
  { 
    return this->m_Degree;
  }
  
  /// Apply transformation to vector.
  virtual Self::SpaceVectorType Apply ( const Self::SpaceVectorType& v ) const
  {
    // initialize result vector
    Self::SpaceVectorType result = Self::SpaceVectorType( 0.0 );

    // now apply actual monomials
    size_t paramIdx = 0;
    for ( size_t monomialIdx = 0; monomialIdx < this->m_NumberOfMonomials; ++monomialIdx )
      {
      const Types::Coordinate monomialValue = Polynomial<4,Types::Coordinate>::EvaluateMonomialAt( monomialIdx, v[0], v[1], v[2] );
      for ( size_t dim = 0; dim < 3; ++dim, ++paramIdx )
	result[dim] += this->m_Parameters[paramIdx] * monomialValue;
      }

    return result;
  }

  /// Get monomial at coordinate.
  Types::Coordinate GetMonomialAt( const size_t idx /*!< Index of the parameter to gert the monomial for. */, const Self::SpaceVectorType& v ) const
  {
    // first three are constant (translational components)
    return Polynomial<4,Types::Coordinate>::EvaluateMonomialAt( (idx/3), v[0], v[1], v[2] );
  }

  /** Return inverse-transformed vector.
   */
  virtual bool ApplyInverse ( const Self::SpaceVectorType&, Self::SpaceVectorType&, const Types::Coordinate = 0.01  ) const { return false; }

  /// Get local Jacobian.
  virtual CoordinateMatrix3x3 GetJacobian( const Self::SpaceVectorType& v ) const;

  /// Compute Jacobian determinant at a certain location.
  virtual Types::Coordinate GetJacobianDeterminant ( const Self::SpaceVectorType& v ) const { return this->GetJacobian( v ).Determinant() ; }

protected:
  /// Polynomial degree.
  byte m_Degree;

  /// Number of monomials (per spatial dimension)
  size_t m_NumberOfMonomials;

  /// Actual virtual clone constructor function.
  virtual Self* CloneVirtual () const 
  {
    return new Self( *this );
  }
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkPolynomialXform_h_included_
