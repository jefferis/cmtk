/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkXform_h_included_
#define __cmtkXform_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMetaInformationObject.h>
#include <Base/cmtkAnatomicalOrientationBase.h>
#include <Base/cmtkVector.h>
#include <Base/cmtkFixedVector.h>
#include <Base/cmtkBitVector.h>
#include <Base/cmtkLandmarkPairList.h>

#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** General 3D coordinate transformation.
 */
class Xform :
  /// Inherit from meta data information container.
  public MetaInformationObject
{
public:
  /// This class.
  typedef Xform Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const for this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Three-dimensional vector type.
  typedef FixedVector<3,Types::Coordinate> SpaceVectorType;

  /// Pointer to warp parameter array.
  Types::Coordinate *m_Parameters;

  /// Total number of parameters, ie. the values in Coefficients.
  size_t m_NumberOfParameters;

  /// Copy constructor.
  Xform( const Xform& other )
    : MetaInformationObject( other ),
      m_NumberOfParameters( other.m_NumberOfParameters ),
      m_ParameterVector( other.m_ParameterVector )
  {
    this->m_Parameters = this->m_ParameterVector->Elements;
    this->SetMetaInfo( cmtk::META_SPACE, cmtk::AnatomicalOrientationBase::ORIENTATION_STANDARD );
  }

  /// Default constructor.
  Xform()
    : m_Parameters( NULL ),
      m_NumberOfParameters( 0 ) 
  {
    this->SetMetaInfo( cmtk::META_SPACE, cmtk::AnatomicalOrientationBase::ORIENTATION_STANDARD );
  }

  /// Virtual destructor.
  virtual ~Xform() {}
  
  /// Check whether coordinate is in domain of transformation.
  virtual bool InDomain( const Self::SpaceVectorType& ) const { return true; }
  
  /// Get global scaling factor.
  virtual Types::Coordinate GetGlobalScaling() const { return 0.0; }

  /// Apply transformation to vector.
  virtual Self::SpaceVectorType Apply ( const Self::SpaceVectorType& ) const = 0;

  /// Apply transformation to vector in-place.
  virtual void ApplyInPlace ( Self::SpaceVectorType& ) const = 0;

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverse ( const Self::SpaceVectorType&, Self::SpaceVectorType&, const Types::Coordinate = 0.01  ) const = 0;

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverseInPlace( Self::SpaceVectorType&, const Types::Coordinate = 0.01  ) const = 0;

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverseInPlaceWithInitial( Self::SpaceVectorType& v, const Self::SpaceVectorType&, const Types::Coordinate error = 0.01 ) const
  {
    return this->ApplyInverseInPlace( v, error );
  }

  /// Clone and return smart pointer.
  Self::SmartPtr Clone () const 
  {
    return Self::SmartPtr( this->CloneVirtual() );
  }

  /// Return number of coefficients in parameter vector.
  virtual size_t ParamVectorDim () const 
  {
    return this->m_NumberOfParameters;
  }

  /** Get number of variable parameters in parameter vector.
   * The variable parameters are those that may be modified i.e. for an
   * optimization. They are located at the beginning of the complete parameter
   * vector.
   */
  virtual size_t VariableParamVectorDim () const
  {
    return this->ParamVectorDim();
  }

  /** Set Xform by parameter vector.
   * Be careful: This is NOT a one-way function. The Xform object may change
   * the parameter vector in order to ensure internal consistency 
   * (AffineXform) or to enhance efficiency.
   */
  virtual void SetParamVector ( CoordinateVector& v );

  /** Copy parameter vector from other transformation.
   * THERE ARE NO CHECKS WHETHER THE TWO TRANSFORMATIONS MATCH!!
   */
  virtual void CopyParamVector ( const Xform* other )
  {
    *(this->m_ParameterVector) = *(other->m_ParameterVector);
    this->m_Parameters = this->m_ParameterVector->Elements;
  }

  /// Set the parameter vector.
  virtual void SetParamVector ( const CoordinateVector& v );

  /// Set a single parameter value.
  virtual void SetParameter ( const size_t idx, const Types::Coordinate p )
  {
    this->m_Parameters[idx] = p;
  }

  /// Get a single parameter value.
  virtual Types::Coordinate GetParameter ( const size_t idx ) const
  {
    return this->m_Parameters[idx];
  }

  /// Copy parameter vector to existing vector object.
  virtual CoordinateVector& GetParamVector( CoordinateVector& v, const size_t targetOffset = 0 ) const;

  /// Get parameter step given a transformed volume size.
  virtual Types::Coordinate GetParamStep( const size_t, const Self::SpaceVectorType&, const Types::Coordinate step_mm = 1 ) const
  { 
    return step_mm;
  }
  
  /// Compute Jacobian determinant at a certain location.
  virtual Types::Coordinate GetJacobianDeterminant ( const Self::SpaceVectorType& ) const = 0;
  
  /** Return registration error for set of source/target landmarks.
   * What is actually returned is the mean squared distance of source
   * landmark after transformation and desired target landmark.
   */
  virtual Types::Coordinate GetLandmarksMSD( const LandmarkPairList& ll ) const;
  
protected:
  /** Encapsulated representation of the transformation parameters.
   * This vector object contains the parameter array pointed at by the public
   * member Coefficients. The latter is used for more efficient direct access
   * to the parameters where necessary.
   */
  CoordinateVector::SmartPtr m_ParameterVector;

  /** Allocate parameter vector.
   */
  void AllocateParameterVector( const size_t numberOfParameters );

  /// Actual virtual clone constructor function.
  virtual Self* CloneVirtual () const = 0;
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkXform_h_included_
