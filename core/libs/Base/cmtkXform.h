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

#ifndef __cmtkXform_h_included_
#define __cmtkXform_h_included_

#include <cmtkconfig.h>

#include <stdlib.h>

#include <cmtkInformationObject.h>

#include <cmtkVector3D.h>
#include <cmtkVector.h>
#include <cmtkBitVector.h>

#include <cmtkMatchedLandmarkList.h>

#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

// Forward declaration of "Volume" class.
// this is required to be able to declare the "Volume*" parameter of
// RegisterVolume
class UniformVolume;

/** General 3D coordinate transformation.
 */
class Xform :
  /// Inherit from meta data information container.
  public InformationObject
{
public:
  /// Smart pointer to Xform.
  typedef SmartPointer<Xform> SmartPtr;

  /// Pointer to warp parameter array.
  Types::Coordinate *m_Parameters;

  /// Total number of parameters, ie. the values in Coefficients.
  size_t m_NumberOfParameters;

  /// Copy constructor.
  Xform( const Xform& other )
    : InformationObject( other )
  {    
  }

  /// Default constructor.
  Xform()
  {
  }

  /// Virtual destructor.
  virtual ~Xform();
  
  /// Check whether coordinate is in domain of transformation.
  virtual bool InDomain( const Vector3D& ) const { return true; }
  
  /// Get global scaling factor.
  virtual Types::Coordinate GetGlobalScaling() const { return 0.0; }

  /// Apply transformation to vector.
  virtual Vector3D Apply ( const Vector3D& ) const = 0;

  /// Apply transformation to vector in-place.
  virtual void ApplyInPlace ( Vector3D& ) const = 0;

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverse ( const Vector3D&, Vector3D&, const Types::Coordinate = 0.01  ) const
  {
    throw Exception( "unimplemented function called" );
  }

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverseInPlace( Vector3D&, const Types::Coordinate = 0.01  ) const
  {
    throw Exception( "unimplemented function called" );
  }

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverseInPlaceWithInitial( Vector3D&, const Vector3D&, const Types::Coordinate = 0.01 ) const
  {
    throw Exception( "unimplemented function called" );
  }

  /// Clone transformation.
  virtual Xform* Clone () const { return NULL; }

  /// Make inverse transformation, or return NULL if this transformation does not support explicit inverses.
  virtual Xform* MakeInverse () const { return NULL; }

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
  virtual size_t VariableParamVectorDim () const = 0;

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

  virtual Types::Coordinate GetParamStep( const size_t, const Types::Coordinate*, const Types::Coordinate step_mm = 1 ) const
  { 
    return step_mm;
  }
  
  /// Compute Jacobian determinant at a certain location.
  virtual Types::Coordinate GetJacobianDeterminant ( const Vector3D& ) const 
  { 
    return 1;
  }
  
  /// Compute sequence of Jacobian determinants from given grid location.
  virtual void GetJacobianDeterminantSequence( double *const values, const int, const int, const int, const size_t numberOfPoints = 1 ) const
  {
    for ( size_t i = 0; i < numberOfPoints; ++i ) 
      values[i] = 1.0;
  }

  /** Return registration error for set of source/target landmarks.
   * What is actually returned is the mean squared distance of source
   * landmark after transformation and desired target landmark.
   */
  virtual Types::Coordinate GetLandmarksMSD( const MatchedLandmarkList* ll ) const;
  
  /** Return derivative of registration error with respect to one parameter.
   */
  virtual void GetDerivativeLandmarksMSD( double&, double&, const MatchedLandmarkList*, const unsigned int, const Types::Coordinate ) {}
  
  /// Get volume influenced by one parameter.
  virtual void GetVolumeOfInfluence( const size_t idx, const Vector3D&, const Vector3D&, Vector3D&, Vector3D&, const int = -1 ) const;
  
  /// Register the grid points of the deformed uniform or non-uniform volume.
  virtual void RegisterVolume ( const UniformVolume* );

  /// Unegister axes points, ie free all internal data structures.
  virtual void UnRegisterVolume () {};

  virtual void GetTransformedGrid( Vector3D&, const int, const int, const int ) const = 0;
  
  virtual void GetTransformedGridSequence( Vector3D *const, const int, const int, const int, const int ) const = 0;

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
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkXform_h_included_
