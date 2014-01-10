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
#include <System/cmtkException.h>

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
  
  /// Exception class thrown when unsupported degree is used.
  class DegreeUnsupported : public Exception {};

  /// Copy constructor.
  PolynomialXform( const PolynomialXform& other )
    : Xform( other ),
      m_Degree( other.m_Degree )
  {
  }
  
  /// Default constructor.
  PolynomialXform( const byte degree = 0 /*!< Polynomial degree - 0 through 4 are supported. */ ) : m_Degree( degree )
  {
    // Sort out how many monomials a polynomial of the given degree has.
    size_t numberOfMonomials = 0;
    switch ( this->m_Degree )
      {
      case 1: numberOfMonomials = Polynomial<1>::NumberOfMonomials; break;
      case 2: numberOfMonomials = Polynomial<2>::NumberOfMonomials; break;
      case 3: numberOfMonomials = Polynomial<3>::NumberOfMonomials; break;
      case 4: numberOfMonomials = Polynomial<4>::NumberOfMonomials; break;
      default: throw Self::DegreeUnsupported();
      }

    // add constant parameter (shift) to polynomial; allocate one set per spatial dimension
    this->AllocateParameterVector( 3 * (1 + numberOfMonomials) );
  }

  /// Clone and return smart pointer.
  Self::SmartPtr Clone () const 
  {
    return Self::SmartPtr( this->CloneVirtual() );
  }

  /// Virtual destructor.
  virtual ~PolynomialXform() {}
  
  /// Apply transformation to vector.
  virtual Self::SpaceVectorType Apply ( const Self::SpaceVectorType& v ) const
  {
    // first three parameters are constant
    Self::SpaceVectorType result;
    for ( size_t idx = 0; idx < 3; ++idx )
      result[idx] = this->m_Parameters[idx];

    // now apply actual monomials
    size_t monomialIdx = 0;
    for ( size_t idx = 3; idx < this->m_NumberOfParameters; ++monomialIdx )
      {
      for ( size_t dim = 0; dim < 3; ++dim, ++idx )
	result[dim] += this->m_Parameters[idx] * Polynomial<4,Types::Coordinate>::EvaluateMonomialAt( monomialIdx, v[0], v[1], v[2] );
      }

    return result;
  }

  /// Get monomial at coordinate.
  Types::Coordinate GetMonomialAt( const size_t idx /*!< Index of the parameter to gert the monomial for. */, const Self::SpaceVectorType& v ) const
  {
    // first three are constant (translational components)
    if ( idx < 3 )
      {
      return 1.0;
      }

    return Polynomial<4,Types::Coordinate>::EvaluateMonomialAt( (idx/3)-1, v[0], v[1], v[2] );
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

  /// Actual virtual clone constructor function.
  virtual Self* CloneVirtual () const 
  {
    return new Self( *this );
  }
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkPolynomialXform_h_included_
