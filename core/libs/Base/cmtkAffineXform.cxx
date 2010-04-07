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

#include <cmtkAffineXform.h>

#include <cmtkMathUtil.h>
#include <cmtkUniformVolume.h>

#include <cmtkConsole.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
AffineXform::MakeIdentityXform () 
{
  this->m_ParameterVector->Clear();
  Types::Coordinate* scales = this->RetScales();
  if ( ! m_LogScaleFactors )
    scales[0] = scales[1] = scales[2] = 1;
  this->ComposeMatrix();
}

AffineXform::AffineXform ( const AffineXform& other ) :
  Xform( other ),
  Matrix( NULL ), 
  m_LogScaleFactors( false ),
  InverseXform( NULL )
{ 
  this->AllocateParameterVector( TotalNumberOfParameters );
  (*this->m_ParameterVector) = (*other.m_ParameterVector);
  this->m_LogScaleFactors = other.m_LogScaleFactors;
  this->NumberDOFs = other.NumberDOFs;
  this->ComposeMatrix();
  memset( RegisteredVolumeAxes, 0, sizeof( RegisteredVolumeAxes ) );
}

AffineXform::AffineXform
( const Types::Coordinate matrix[4][4], const Types::Coordinate* center ) :
  Matrix( &matrix[0][0] ), 
  m_LogScaleFactors( false ),
  InverseXform( NULL )
{
  this->AllocateParameterVector( TotalNumberOfParameters );
  this->NumberDOFs = this->DefaultNumberOfDOFs();
  if ( center )
    memcpy( this->RetCenter(), center, 3 * sizeof( Types::Coordinate ) );
  else
    memset( this->RetCenter(), 0, 3 * sizeof( Types::Coordinate ) );
  this->DecomposeMatrix();
  memset( RegisteredVolumeAxes, 0, sizeof( RegisteredVolumeAxes ) );
}

AffineXform::AffineXform
( const MatrixType& matrix, const Types::Coordinate* center ) :
  Matrix( matrix ), 
  m_LogScaleFactors( false ),
  InverseXform( NULL )
{
  this->AllocateParameterVector( TotalNumberOfParameters );
  this->NumberDOFs = this->DefaultNumberOfDOFs();
  if ( center )
    memcpy( this->RetCenter(), center, 3 * sizeof( Types::Coordinate ) );
  else
    memset( this->RetCenter(), 0, 3 * sizeof( Types::Coordinate ) );
  this->DecomposeMatrix();
  memset( RegisteredVolumeAxes, 0, sizeof( RegisteredVolumeAxes ) );
}

AffineXform::AffineXform
( const Types::Coordinate matrix[4][4], const Types::Coordinate xlate[3],
  const Types::Coordinate center[3] ) :
  Matrix( &matrix[0][0] ), 
  m_LogScaleFactors( false ),
  InverseXform( NULL )
{
  this->AllocateParameterVector( TotalNumberOfParameters );
  this->NumberDOFs = this->DefaultNumberOfDOFs();
  Types::Coordinate cM[3] = 
    {
      center[0]*Matrix[0][0] + center[1]*Matrix[1][0] + center[2]*Matrix[2][0],
      center[0]*Matrix[0][1] + center[1]*Matrix[1][1] + center[2]*Matrix[2][1],
      center[0]*Matrix[0][2] + center[1]*Matrix[1][2] + center[2]*Matrix[2][2]
    };
  
  Matrix[3][0] = xlate[0] + center[0] - cM[0];
  Matrix[3][1] = xlate[1] + center[1] - cM[1];
  Matrix[3][2] = xlate[2] + center[2] - cM[2];

  this->Matrix.Decompose( this->m_Parameters, center, this->m_LogScaleFactors );
  memset( RegisteredVolumeAxes, 0, sizeof( RegisteredVolumeAxes ) );
}

AffineXform& 
AffineXform::operator=( const AffineXform& other )
{
  (*this->m_ParameterVector) = (*other.m_ParameterVector);
  this->NumberDOFs = other.NumberDOFs;
  this->m_LogScaleFactors = other.m_LogScaleFactors;
  this->ComposeMatrix();
  return *this;
}

void
AffineXform::SetNumberDOFs ( const int numberDOFs ) 
{
  this->NumberDOFs = numberDOFs; 
  if ( this->NumberDOFs == 7 )
    {
    this->m_Parameters[8] = (this->m_Parameters[7] = this->m_Parameters[6]);
    this->ComposeMatrix();
    }
}

void
AffineXform::SetUseLogScaleFactors( const bool logScaleFactors )
{
  if ( logScaleFactors != this->m_LogScaleFactors )
    {
    if ( logScaleFactors )
      {
      for ( int i = 6; i < 9; ++i )
	this->m_Parameters[i] = log( this->m_Parameters[i] );
      }
    else
      {
      for ( int i = 6; i < 9; ++i )
	this->m_Parameters[i] = exp( this->m_Parameters[i] );	  
      }
    this->m_LogScaleFactors = logScaleFactors;
    }
}

void
AffineXform::ComposeMatrix () 
{
  // For 7 parameter form (rigid plus global scaling) be sure to use equal
  // scalings for all coordinates.
  if ( this->NumberDOFs == 7 ) 
    this->m_Parameters[8] = (this->m_Parameters[7] = this->m_Parameters[6]);
    
  // Now build matrix.
  this->Matrix.Compose( this->m_Parameters, this->m_LogScaleFactors );
  this->UpdateInverse();
}

bool
AffineXform::DecomposeMatrix () 
{
  return this->Matrix.Decompose( this->m_Parameters, this->RetCenter(), this->m_LogScaleFactors );
}

void
AffineXform::RotateScaleShear
( Types::Coordinate vM[3], const Types::Coordinate v[3] ) const
{
  Types::Coordinate V[3];
  V[0] = v[0] * Matrix[0][0] + v[1] * Matrix[1][0] + v[2] * Matrix[2][0];
  V[1] = v[0] * Matrix[0][1] + v[1] * Matrix[1][1] + v[2] * Matrix[2][1];
  V[2] = v[0] * Matrix[0][2] + v[1] * Matrix[1][2] + v[2] * Matrix[2][2];
  memcpy( vM, V, 3 * sizeof( Types::Coordinate ) );
}

AffineXform* 
AffineXform::MakeInverse () const
{
  Self* inverseXform = new AffineXform();
  inverseXform->m_LogScaleFactors = this->m_LogScaleFactors;
  inverseXform->SetNumberDOFs( this->NumberDOFs );
  inverseXform->Matrix = this->Matrix;
  inverseXform->Matrix.Invert();
  inverseXform->DecomposeMatrix();

  Types::Coordinate newCenter[3];
  this->Matrix.Multiply( this->RetCenter(), newCenter );
  inverseXform->ChangeCenter( newCenter );

  if ( this->NumberDOFs == 7 ) 
    {
    inverseXform->m_Parameters[8] = (inverseXform->m_Parameters[7] = inverseXform->m_Parameters[6]);
    inverseXform->Matrix.Compose( inverseXform->m_Parameters, this->m_LogScaleFactors );
    }

  inverseXform->m_MetaInformation[META_SPACE] = this->m_MetaInformation[META_SPACE];

  return inverseXform;
}

void
AffineXform::ChangeCenter ( const Types::Coordinate* newCenter ) 
{
  Types::Coordinate* xlate = this->RetXlate();
  Types::Coordinate* center = this->RetCenter();
  Types::Coordinate deltaCenter[3] = { newCenter[0] - center[0], newCenter[1] - center[1], newCenter[2] - center[2] };
  xlate[0] -= deltaCenter[0];
  xlate[1] -= deltaCenter[1];
  xlate[2] -= deltaCenter[2];

  this->RotateScaleShear( deltaCenter, deltaCenter );
  
  xlate[0] += deltaCenter[0];
  xlate[1] += deltaCenter[1];
  xlate[2] += deltaCenter[2];

  memcpy( center, newCenter, 3*sizeof(Types::Coordinate) );
}

void
AffineXform::SetMatrix( const MatrixType& matrix ) 
{
  this->Matrix = matrix;
  this->DecomposeMatrix();
  this->UpdateInverse();
}

void
AffineXform::SetMatrix( const float matrix[4][4] ) 
{
  for ( unsigned int j = 0; j < 4; ++j )
    for ( unsigned int i = 0; i < 4; ++i )
      Matrix[j][i] = matrix[j][i];
  this->DecomposeMatrix();
  this->UpdateInverse();
}

void
AffineXform::SetMatrix( const double matrix[4][4] ) 
{
  for ( unsigned int j = 0; j < 4; ++j )
    for ( unsigned int i = 0; i < 4; ++i )
      Matrix[j][i] = matrix[j][i];
  this->DecomposeMatrix();
  this->UpdateInverse();
}

template<> 
void
AffineXform::GetMatrix( float (&matrix)[4][4] ) const
{
  for ( unsigned int j = 0; j < 4; ++j )
    for ( unsigned int i = 0; i < 4; ++i )
      matrix[j][i] = static_cast<float>( Matrix[j][i] );
}

template<> 
void 
AffineXform::GetMatrix( double (&matrix)[4][4] ) const
{
  for ( unsigned int j = 0; j < 4; ++j )
    for ( unsigned int i = 0; i < 4; ++i )
      matrix[j][i] = static_cast<double>( Matrix[j][i] );
}

void
AffineXform::SetParamVector ( CoordinateVector& v )
{
  Superclass::SetParamVector( v );
  this->CanonicalRotationRange();
  this->ComposeMatrix();
  v = (*this->m_ParameterVector);
}

void 
AffineXform::SetParamVector ( const CoordinateVector& v )
{
  Superclass::SetParamVector( v );
  this->CanonicalRotationRange();
  this->ComposeMatrix();
}

void
AffineXform::SetParameter ( const size_t idx, const Types::Coordinate p )
{
  Superclass::SetParameter( idx, p );
  this->CanonicalRotationRange();
  this->ComposeMatrix();
}


Types::Coordinate
AffineXform::GetParamStep
( const size_t idx, const Types::Coordinate* volSize, const Types::Coordinate mmStep ) 
  const
{
  if ( (int)idx >= this->NumberDOFs ) return 0.0;
  
  switch ( idx )
    {
    case 0: case 1: case 2:
      return mmStep;
    case 3:
      return mmStep * 180 / (M_PI * sqrt( MathUtil::Square( volSize[1] ) + MathUtil::Square( volSize[2] ) ) );
    case 4:
      return mmStep * 180 / (M_PI * sqrt( MathUtil::Square( volSize[0] ) + MathUtil::Square( volSize[2] ) ) );
    case 5:
      return mmStep * 180 / (M_PI * sqrt( MathUtil::Square( volSize[0] ) + MathUtil::Square( volSize[1] ) ) );
    case 6: case 7: case 8:
      if ( this->m_LogScaleFactors )
	return log( 1 + 0.5 * mmStep / MathUtil::Max( 3, volSize ) );
      else
	return 0.5 * mmStep / MathUtil::Max( 3, volSize );
    case 9: case 10: case 11:
      return 0.5 * mmStep / MathUtil::Max( 3, volSize );
    }
  return mmStep;
}

AffineXform::SmartPtr
AffineXform::GetDifference( const AffineXform& other ) const
{
  Self::SmartPtr result( this->MakeInverse() );
  result->Concat( other );
  return result;
}

void
AffineXform::Concat( const AffineXform& other )
{
  this->Matrix *= other.Matrix;
  this->DecomposeMatrix();
}

void
AffineXform::Insert( const AffineXform& other )
{
  Self::MatrixType composed( other.Matrix );
  composed *= this->Matrix;
  this->Matrix = composed;
  this->DecomposeMatrix();
}

void 
AffineXform::RotateWXYZ
( const Types::Coordinate angle, const Vector3D& direction,
  const Types::Coordinate* center, Self::MatrixType *const accumulate )
{
  Vector3D unit( direction );

  Vector3D center3D;
  if ( center ) 
    center3D.Set( center );
  else
    center3D.Set( this->RetCenter() );

  if ( accumulate ) 
    {
    unit += center3D;
    accumulate->Multiply( unit.XYZ );
    accumulate->Multiply( center3D.XYZ );
    unit -= center3D;
    }

  // translation into rotation center
  Self::MatrixType xlate;
  for ( int dim = 0; dim < 3; ++dim )
    xlate[3][dim] = -center3D[dim];

  if ( accumulate ) 
    {
    *accumulate *= xlate;
    }

  this->Matrix *= xlate;

  double x = unit.XYZ[0];
  double y = unit.XYZ[1];
  double z = unit.XYZ[2];

  // make a normalized quaternion
  const double w = cos(0.5*angle);
  const double f = sin(0.5*angle)/sqrt(x*x+y*y+z*z);
  x *= f;
  y *= f;
  z *= f;

  // convert the quaternion to a matrix
  Self::MatrixType matrix;

  const double ww = w*w;
  const double wx = w*x;
  const double wy = w*y;
  const double wz = w*z;

  const double xx = x*x;
  const double yy = y*y;
  const double zz = z*z;

  const double xy = x*y;
  const double xz = x*z;
  const double yz = y*z;

  const double s = ww - xx - yy - zz;

  matrix[0][0] = xx*2 + s;
  matrix[1][0] = (xy + wz)*2;
  matrix[2][0] = (xz - wy)*2;

  matrix[0][1] = (xy - wz)*2;
  matrix[1][1] = yy*2 + s;
  matrix[2][1] = (yz + wx)*2;

  matrix[0][2] = (xz + wy)*2;
  matrix[1][2] = (yz - wx)*2;
  matrix[2][2] = zz*2 + s;

  this->Matrix *= matrix;

  xlate.Invert();
  this->Matrix *= xlate;
  this->DecomposeMatrix();

  if ( accumulate ) 
    {
    *accumulate *= matrix;
    *accumulate *= xlate;
    }
}

AffineXform::SmartPtr
AffineXform::GetInverse() const
{
  if ( !InverseXform ) 
    {
    InverseXform = AffineXform::SmartPtr( this->MakeInverse() );
    } 
  else 
    {
    this->UpdateInverse();
    }
  
  return InverseXform;
}

void
AffineXform::UpdateInverse() const
{
  if ( InverseXform ) 
    {
    InverseXform->NumberDOFs = this->NumberDOFs;
    InverseXform->m_LogScaleFactors = this->m_LogScaleFactors;
    InverseXform->Matrix = this->Matrix;
    InverseXform->Matrix.Invert();
    InverseXform->DecomposeMatrix();
    }
}

void
AffineXform::CanonicalRotationRange()
{
  for ( int rotIdx = 0; rotIdx<3; ++rotIdx ) 
    {
    while ( this->m_Parameters[3+rotIdx] >  180 ) 
      this->m_Parameters[3+rotIdx] -= 360;
    while ( this->m_Parameters[3+rotIdx] < -180 ) 
      this->m_Parameters[3+rotIdx] += 360;
    }
}

void
AffineXform::RegisterVolume ( const UniformVolume* volume )
{
  this->UnRegisterVolume();

  // define volume corners
  Vector3D dX(1,0,0), dY(0,1,0), dZ(0,0,1);
  Vector3D V(0,0,0);
    
  // compute inverse transformation for subsequent (direct) use.
  this->GetInverse();

  InverseXform->ApplyInPlace(V);
  InverseXform->ApplyInPlace(dX);
  dX -= V;
  InverseXform->ApplyInPlace(dY);
  dY -= V;
  InverseXform->ApplyInPlace(dZ);
  dZ -= V;
  
  // Apply post-transformation scaling
  /**
     // This code allows for the integration of an implicit coordinate-to-
     // grid-index scaling. Potentially, this saves some additional computation
     // time, but it does not really fit into the unified handling of 
     // registered volumes for affine and non-rigid transformations. Therefore,
     // it has been excluded from the current code base.
  if ( scales ) {
    dX.Set( dX[0] * scales[0], dX[1] * scales[1], dX[2] * scales[2] );
    dY.Set( dY[0] * scales[0], dY[1] * scales[1], dY[2] * scales[2] );
    dZ.Set( dZ[0] * scales[0], dZ[1] * scales[1], dZ[2] * scales[2] );
    V.Set (  V[0] * scales[0],  V[1] * scales[1],  V[2] * scales[2] );
  }
  */

  const Types::Coordinate deltaX = volume->m_Delta[0];
  const Types::Coordinate deltaY = volume->m_Delta[1];
  const Types::Coordinate deltaZ = volume->m_Delta[2];

  const int *volumeDims = volume->GetDims();
  for ( int dim = 0; dim<3; ++dim ) 
    {
    RegisteredVolumeAxes[dim] = Memory::AllocateArray<Vector3D>( volumeDims[dim] );
    }
  
  int idx;
  for ( idx = 0; idx < volumeDims[0]; ++idx )
    RegisteredVolumeAxes[0][idx] = deltaX*idx*dX;
  for ( idx = 0; idx < volumeDims[1]; ++idx )
    RegisteredVolumeAxes[1][idx] = deltaY*idx*dY;
  for ( idx = 0; idx < volumeDims[2]; ++idx )
    (RegisteredVolumeAxes[2][idx] = deltaZ*idx*dZ) += V;
}

void
AffineXform::UnRegisterVolume ()
{
  for ( int idx=0; idx<3; ++idx ) 
    {
    if ( RegisteredVolumeAxes[idx] ) 
      {
      delete[] RegisteredVolumeAxes[idx];
      RegisteredVolumeAxes[idx] = NULL;
      }
    }
}

void
AffineXform::Print() const
{
  StdErr.printf( "AffineXform at %p:\n", this );
  StdErr.printf( "\tNumber DOFs: %d\n", NumberDOFs );
  StdErr.printf( "\tTranslation: [%f,%f,%f]\n", this->m_Parameters[0], this->m_Parameters[1], this->m_Parameters[2] );
  StdErr.printf( "\tRotation: [%f,%f,%f] around [%f,%f,%f]\n", 
		    this->m_Parameters[3], this->m_Parameters[4], this->m_Parameters[5], 
		    this->m_Parameters[12], this->m_Parameters[13], this->m_Parameters[14] );
  StdErr.printf( "\tScale: [%f,%f,%f]\n", this->m_Parameters[6], this->m_Parameters[7], this->m_Parameters[8] );
  StdErr.printf( "\tShear: [%f,%f,%f]\n", this->m_Parameters[9], this->m_Parameters[10], this->m_Parameters[11] );

  this->Matrix.Print( StdErr );
}

void 
AffineXform::ApplyToAll
( CoordinateVector& v, BitVector& valid, const bool inverse, const Types::Coordinate epsilon, const int* gridDims )
{
  if ( inverse )
    this->GetInverse()->ApplyToAll( v, valid, false, epsilon, gridDims );
  else
    {
    const size_t numberOfPoints = v.Dim / 3;
    
    Types::Coordinate* p = v.Elements;
    for ( size_t i = 0; i < numberOfPoints; ++i, p+=3 )
      if ( valid[i] )
	this->Matrix.Multiply( p );
    }
}

} // namespace
