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

#ifndef __cmtkAffineXform_h_included_
#define __cmtkAffineXform_h_included_

#include <cmtkconfig.h>

#include <cmtkTypes.h>
#include <cmtkVector.h>
#include <cmtkVector3D.h>
#include <cmtkXform.h>
#include <cmtkMatrix4x4.h>

#include <cmtkSmartPtr.h>

#include <cmtkUnits.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** 3D affine transformation.
 * This transformation class allows translations, rotations, componentwise 
 * scalings, and shears. Transformation is done by vector-matrix-multiplication
 * with a homogeneous 4x4 matrix. The transformation can be defined by either
 * the matrix itself or a parameter vector. Both representations are 
 * permanently accessible and held mutually up-to-date whenever one of them
 * changes.
 *@author $Author$
 */
class AffineXform : 
  /// Inherit virtual interface from generic transformation.
  public Xform 
{
public:
  /// This class type.
  typedef AffineXform Self;

  /// Superclass type.
  typedef Xform Superclass;

  /// Smart pointer to AffineXform.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const AffineXform.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Transformation matrix.type.
  typedef Matrix4x4<Types::Coordinate> MatrixType;

  /// Create identity transformation.
  void MakeIdentityXform ();

protected:
  /** Compose this object's transformation matrix from parameter vector.
   */
  void ComposeMatrix();

  /** Decompose this object's transformation matrix into parameter vector.
   */
  bool DecomposeMatrix();

  /** Actual number of degrees of freedom.
   * This value should be one out of four choices: 6 (rigid transformation), 
   * 7 (rigid with global scaling), 9 (rigid with componentwise scaling), or
   * 12 (full affine). The number of elements of the parameter vector does NOT
   * change, even in the 7 DOF case. Then, all scaling factors are simply kept
   * equal.
   */
  int NumberDOFs;

public:
  /** Homogeneous transformation matrix.
   * Vectors are transformed by right-multiplication with this matrix, i.e.
   * v being a vector, the transformed vector is (v Matrix).
   */
  Self::MatrixType Matrix;

  /** Total number of parameters.
   * The values stored in the parameter vector are:
   *  Three translations (x, y, z), three rotation
   * angles (around x-, y-, and z-axis) in degrees, three scaling factors
   * (x-, y-, and z-direction), and three shears. The last three values
   * define the center of rotation, scaling, and shearing.
   */
  static const size_t TotalNumberOfParameters = 15;

  /** Apply rotation, scaling, and shearing but no translation to vector.
   * This function only used the upper left 3x3 sub-matrix of Matrix to
   * transform a given vector. This is needed for rotation-center changes.
   *\note It is safe to use the same address for source and destination of this
   * operation.
   *@param vM This three-component array holds the rotated, scaled, and
   * sheared vector after return from the function.
   *@param v The vector to be rotated, scaled, and sheared.
   */
  void RotateScaleShear ( Types::Coordinate[3], const Types::Coordinate[3] ) const;

  /** Create identity transformation.
   */
  AffineXform () :
    m_LogScaleFactors( false )
  {
    this->AllocateParameterVector( TotalNumberOfParameters );
    this->NumberDOFs = this->DefaultNumberOfDOFs();
    this->MakeIdentityXform();
  }

  /** Create transformation from parameter vector.
   *@param v The parameter vector defining the desired transformation. Refer
   * to 'Parameters' for a detailed description.
   *@see Parameters
   */
  AffineXform ( const CoordinateVector& v, const bool logScaleFactors = false ) :
      m_LogScaleFactors( logScaleFactors )
  {
    this->AllocateParameterVector( TotalNumberOfParameters );
    this->NumberDOFs = this->DefaultNumberOfDOFs();
    this->SetParamVector( v );
  }

  /** Create transformation from raw parameter array.
   *@param v The parameter vector defining the desired transformation. Refer
   * to 'Parameters' for a detailed description.
   *@see Parameters
   */
  AffineXform ( const Types::Coordinate v[15], const bool logScaleFactors = false ) :
      m_LogScaleFactors( logScaleFactors )
  {
    this->AllocateParameterVector( TotalNumberOfParameters );
    this->NumberDOFs = this->DefaultNumberOfDOFs();
    memcpy( this->m_Parameters, v, 15 * sizeof(Types::Coordinate) );
    this->ComposeMatrix();
    this->CanonicalRotationRange();
  }

  /** Create transformation from transformation matrix.
   *@param matrix The homogeneous 4x4 affine transformation matrix.
   *@param center If non-NULL, this parameter points to a three-coordinate
   * vector defining the rotation center of the constructed transformation.
   * This does only effect the computation of the translation vector when
   * decomposing the given matrix into its parameter representation.
   */
  AffineXform ( const Types::Coordinate matrix[4][4], const Types::Coordinate* center = NULL );

  /** Create transformation from transformation matrix.
   *@param matrix The homogeneous 4x4 affine transformation matrix.
   *@param center If non-NULL, this parameter points to a three-coordinate
   * vector defining the rotation center of the constructed transformation.
   * This does only effect the computation of the translation vector when
   * decomposing the given matrix into its parameter representation.
   */
  AffineXform ( const MatrixType& matrix, const Types::Coordinate* center = NULL );

  /** Create transformation from transformation matrix.
   * Translation vector and rotation center are given in addition to the
   * 4x4 transformation matrix..
   *@param matrix The homogeneous 4x4 affine transformation matrix. If a
   * translation is defined by this matrix, it is ignored.
   *@param xlate Translation of the constructed transformation.
   *@param center Rotation center of the constructed transformation.
   */
  AffineXform ( const Types::Coordinate matrix[4][4], const Types::Coordinate xlate[3], const Types::Coordinate center[3] );
  
  /// Copy transform by reference.
  AffineXform ( const AffineXform& other );
  
  /** Virtual destructor.
   * Frees the linked inverse transformation if one exists.
   */
  virtual ~AffineXform() 
  { 
    InverseXform = Self::SmartPtr::Null; 
  }
  
  /// Clone this object.
  virtual Self* Clone () const
  {
    return new AffineXform( *this );
  }
  
  /// Clone inverse of this transformation.
  virtual Self* MakeInverse () const;

  /** Compute difference to another affine transformation.
   *\todo It is probably not a good idea to return a dynamically allocated object here.
   */
  virtual Self::SmartPtr GetDifference( const AffineXform& other ) const;

  /// Get linked inverse of this transformation.
  const Self::SmartPtr GetInverse() const;

  /// Get global scaling factor.
  virtual Types::Coordinate GetGlobalScaling() const 
  { 
    if ( this->m_LogScaleFactors )
      {
      return exp( this->m_Parameters[6] + this->m_Parameters[7] + this->m_Parameters[8] );
      }
    else
      {
      return this->m_Parameters[6] * this->m_Parameters[7] * this->m_Parameters[8];
      }
  }

  /** Compute Jacobian determinant at a certain location.
   * For an affine transformation, the Jacobian determinant is the product
   * of the anisotropic scale factors at any location.
   */
  virtual Types::Coordinate GetJacobianDeterminant ( const Vector3D& ) const
  { 
    return this->GetGlobalScaling(); 
  }

  /// Concatenate this transformation with another.
  void Concat( const AffineXform& other );

  /// Insert another transformation before this one.
  void Insert( const AffineXform& other );

  /** Rotate around axis.
   *\param origin If this parameter is given, it defines the 3D coordinates of
   * a point on the rotation axis. If this parameter is not given, the
   * rotation axis goes through the rotation center of this transformation.
   *\param accumulate This parameter, if given, points to a 4x4 transformation
   * matrix that can accumulate successive rotations. If given, two things 
   * happen: a) the current rotation axis is rotated according to the 
   * accumulated transformation, and b) the new rotation is concatenated onto
   * the accumulated matrix by multiplication.
   */
  void RotateWXYZ( const Units::Radians angle, const Vector3D& direction, const Types::Coordinate* origin = NULL, Self::MatrixType *const accumulate = NULL );

  /// Change transformation coordinate system.
  void ChangeCoordinateSystem
  ( const Vector3D& newX, const Vector3D newY )
  {
    this->Matrix.ChangeCoordinateSystem( newX.XYZ, newY.XYZ );
    this->DecomposeMatrix();
  }

  /// Apply transformation to vector.
  virtual Vector3D Apply ( const Vector3D& vec ) const 
  {
    Vector3D Result( vec );
    this->ApplyInPlace( Result );
    return Result;
  }
  
  /// Apply transformation to existing vector.
  void ApplyInPlace( Vector3D& vec ) const 
  {
    this->Matrix.Multiply( vec.XYZ );
  }

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverse ( const Vector3D& v, Vector3D& u, const Types::Coordinate = 0.01  ) const
  {
    u = this->GetInverse()->Apply( v );
    return true;
  }

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverseInPlace( Vector3D& v, const Types::Coordinate = 0.01  ) const
  {
    this->GetInverse()->ApplyInPlace( v );
    return true;
  }

  /**@name Read-only parameter retrieval.
   */
  //@{
  /// Return pointer to translation parameters.
  const Types::Coordinate* RetXlate () const  { return this->m_Parameters; }
  /// Return pointer to rotation angles.
  const Types::Coordinate* RetAngles () const { return this->m_Parameters+3; }
  /// Return pointer to scaling factors.
  const Types::Coordinate* RetScales () const { return this->m_Parameters+6; }
  /// Return pointer to shear coefficients.
  const Types::Coordinate* RetShears () const { return this->m_Parameters+9; }
  /// Return pointer to center of rotation, scaling, and shearing.
  const Types::Coordinate* RetCenter () const { return this->m_Parameters+12; }
  //@}

  /**@name Modifyable parameter retrieval.
   */
  //@{
  /// Return pointer to translation parameters.
  Types::Coordinate* RetXlate ()  { return this->m_Parameters; }
  /// Return pointer to rotation angles.
  Types::Coordinate* RetAngles () { return this->m_Parameters+3; }
  /// Return pointer to scaling factors.
  Types::Coordinate* RetScales () { return this->m_Parameters+6; }
  /// Return pointer to shear coefficients.
  Types::Coordinate* RetShears () { return this->m_Parameters+9; }
  /// Return pointer to center of rotation, scaling, and shearing.
  Types::Coordinate* RetCenter () { return this->m_Parameters+12; }
  //@}

  /**@name Direct parameter modifications.
   */
  //@{
  /// Set transformation's translation vector.
  void SetXlate ( const Types::Coordinate* xlate ) 
  { 
    if ( this->RetXlate() != xlate )
      memcpy( this->RetXlate(), xlate, 3 * sizeof(Types::Coordinate) );
    this->ComposeMatrix();
  }

  /// Set transformation's translation vector.
  void SetXlate ( const Types::Coordinate dx, const Types::Coordinate dy, const Types::Coordinate dz ) 
  { 
    this->m_Parameters[0] = dx; this->m_Parameters[1] = dy; this->m_Parameters[2] = dz;
    this->ComposeMatrix(); 
  }

  /// Add to transformation's translation vector.
  void Translate( const Types::Coordinate dx, const Types::Coordinate dy, const Types::Coordinate dz ) 
  {
    this->m_Parameters[0] += dx; this->m_Parameters[1] += dy; this->m_Parameters[2] += dz; 
    this->ComposeMatrix(); 
  }
  
  /// Add to transformation's translation vector.
  void Translate( const Types::Coordinate delta[3] ) 
  { 
    for ( int dim = 0; dim < 3; ++dim ) 
      this->m_Parameters[dim] += delta[dim];
    this->ComposeMatrix();
  }

  /// Set transformation's rotation angles.
  void SetAngles ( const Types::Coordinate* angles ) 
  {
    if ( angles != this->RetAngles() )
      memcpy( this->RetAngles(), angles, 3 * sizeof(Types::Coordinate) );
    this->ComposeMatrix();
  }

  /// Get scale factors with implicit conversion of log scales.
  template<class T> void GetScales( T (&scales)[3] ) const
  { 
    if ( this->m_LogScaleFactors )
      {
      for ( size_t i = 0; i < 3; ++i )
	{
	scales[i] = exp( this->m_Parameters[6+i] );
	}
      }
    else
      {
      for ( size_t i = 0; i < 3; ++i )
	{
	scales[i] = this->m_Parameters[6+i];
	}
      }
  }

  /// Set transformation's scaling factors.
  void SetScales ( const Types::Coordinate* scales ) 
  { 
    if ( this->RetScales() != scales )
      memcpy( this->RetScales(), scales, 3 * sizeof(Types::Coordinate) );
    this->ComposeMatrix();
  }
  
  /// Set transformation's scaling factors.
  void SetScales ( const Types::Coordinate sx, const Types::Coordinate sy, const Types::Coordinate sz )
  { 
    this->m_Parameters[6] = sx; 
    this->m_Parameters[7] = sy; 
    this->m_Parameters[8] = sz; }

  /// Set transformation's shears.
  void SetShears ( const Types::Coordinate* shears ) 
  { 
    if ( this->RetShears() != shears )
      memcpy( this->RetShears(), shears, 3 * sizeof(Types::Coordinate) );
    this->ComposeMatrix();
  }

  /// Set transformation's rotation, scaling, and shearing center.
  void SetCenter ( const Types::Coordinate* center ) 
  {
    if ( this->RetCenter() != center )
      memcpy( RetCenter(), center, 3 * sizeof(Types::Coordinate) );
    this->ComposeMatrix();
  }
  //@}
  
  /**@name Matrix access.
   */
  //@{

  /// Return transformation matrix.
  const Types::Coordinate* RetMatrix () const { return &Matrix[0][0]; }
    
  /// Return transformation matrix.
  Types::Coordinate* RetMatrix () { return &Matrix[0][0]; }

  /// Set transformation matrix.
  void SetMatrix( const float matrix[4][4] );
  void SetMatrix( const double matrix[4][4] );
  void SetMatrix( const MatrixType& matrix );

  /// Get transformation matrix.
  template<class T> void GetMatrix( T (&matrix)[4][4] ) const;
  //@}

  /** Create equivalent transformation with different rotational center.
   * In fact, this function only computes a new translation vector reflecting
   * the changed rotation center. The transformation matrix itself is not
   * changed.
   */
  void ChangeCenter ( const Types::Coordinate* center );

  /// Get dimension of parameter vector.
  virtual size_t ParamVectorDim () const { return 15; }

  /** Get dimension of variable parameter vector.
   * The rotation center is not considered variable, therefore 6 is returned.
   */
  virtual size_t VariableParamVectorDim () const { return NumberDOFs; }

  /// Set the number of degrees of freedom for this object.
  virtual void SetNumberDOFs ( const int numberDOFs = 12 );

  /// Get the number of degrees of freedom for this object.
  virtual unsigned int GetNumberDOFs () const { return NumberDOFs; }

  /// Return flag for log scale factors.
  bool GetUseLogScaleFactors() const
  {
    return this->m_LogScaleFactors;
  }

  /// Switch between log and ordinary scale factors.
  void SetUseLogScaleFactors( const bool logScaleFactors );

  /// Set Xform by parameter vector.
  virtual void SetParamVector ( CoordinateVector& v );

  /** Set Xform by constant parameter vector.
   * Other than setting the transformation from a non-constant vector object,
   * this function does not update the given parameter vector to match the
   * internal representation of the igsAfffineXform object. It is therefore
   * not guaranteed that subsequent calls to GetParamVector() will yield the
   * same parameters.
   */
  virtual void SetParamVector ( const CoordinateVector& v );

  /// Set a single parameter value.
  virtual void SetParameter ( const size_t idx, const Types::Coordinate p );

  /// Get parameter stepping.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate* vol_size, const Types::Coordinate step_mm = 1 ) const;
  
  /// Assignment operator.
  AffineXform& operator=( const AffineXform& other );

  /// Print object.
  virtual void Print() const;

private:
  /// Flag for logarithmic vs. ordinary scale factors.
  bool m_LogScaleFactors;

  /// Return the default number of degrees of freedom.
  static int DefaultNumberOfDOFs() { return 12; }

  /// Link to an auto-updated inverse transformation.
  mutable Self::SmartPtr InverseXform;

  /// Update linked inverse transformation.
  void UpdateInverse() const;

  /** Set transformation matrix.
   * Other than the plain SetMatrix() function, this implementation DOES NOT
   * also update the associated inverse transformation. Therefore, it may
   * safely be called by all other SetXXX() functions in order to update
   * the inverse transformation without producing an infinite loop.
   */
  void SetMatrixDirect ( const Types::Coordinate* matrix ) 
  {
    Matrix.Set( matrix );
    this->DecomposeMatrix();
  }

  /// Correct rotation parameters to canonical range -180 to 180.
  void CanonicalRotationRange();
};

//@}

} // namespace cmtk

#endif // #define __cmtkAffineXform_h_included_
