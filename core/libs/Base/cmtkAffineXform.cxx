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

#include <Base/cmtkAffineXform.h>

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkUniformVolume.h>

#include <System/cmtkConsole.h>

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
  m_LogScaleFactors( false )
{ 
  this->AllocateParameterVector( TotalNumberOfParameters );
  (*this->m_ParameterVector) = (*other.m_ParameterVector);
  this->m_LogScaleFactors = other.m_LogScaleFactors;
  this->NumberDOFs = other.NumberDOFs;
  this->ComposeMatrix();
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

AffineXform::SpaceVectorType
AffineXform::RotateScaleShear
( const Self::SpaceVectorType& v ) const
{
  Self::SpaceVectorType Mv;
  for ( size_t i = 0; i<3; ++i )
    {
    Mv[i] = v[0] * Matrix[0][i] + v[1] * Matrix[1][i] + v[2] * Matrix[2][i];
    }
  return Mv;
}

AffineXform* 
AffineXform::MakeInverse () const
{
  Self* inverseXform = new AffineXform();
  inverseXform->m_LogScaleFactors = this->m_LogScaleFactors;
  inverseXform->SetNumberDOFs( this->NumberDOFs );
  inverseXform->Matrix = this->Matrix.GetInverse();
  inverseXform->DecomposeMatrix();

  const Self::SpaceVectorType newCenter = Self::SpaceVectorType::FromPointer( this->RetCenter() ) * this->Matrix;
  inverseXform->ChangeCenter( newCenter );
  
  if ( this->NumberDOFs == 7 ) 
    {
    inverseXform->m_Parameters[8] = (inverseXform->m_Parameters[7] = inverseXform->m_Parameters[6]);
    inverseXform->Matrix.Compose( inverseXform->m_Parameters, this->m_LogScaleFactors );
    }
  
  inverseXform->CopyMetaInfo( *this, META_SPACE );
  inverseXform->CopyMetaInfo( *this, META_XFORM_FIXED_IMAGE_PATH );
  inverseXform->CopyMetaInfo( *this, META_XFORM_MOVING_IMAGE_PATH );
  
  return inverseXform;
}

void
AffineXform::ChangeCenter ( const Self::SpaceVectorType& newCenter ) 
{
  Types::Coordinate *const xlate = this->RetXlate();
  Types::Coordinate *const center = this->RetCenter();
  Self::SpaceVectorType deltaCenter = newCenter - Self::SpaceVectorType::FromPointer( center );
  
  for ( size_t i = 0; i<3; ++i )
    xlate[i] -= deltaCenter[i];

  deltaCenter = this->RotateScaleShear( deltaCenter );

  for ( size_t i = 0; i<3; ++i )
    {
    xlate[i] += deltaCenter[i];
    center[i] = newCenter[i];
    }
}

void
AffineXform::SetMatrix( const MatrixType& matrix ) 
{
  this->Matrix = matrix;
  this->DecomposeMatrix();
  this->UpdateInverse();
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
( const size_t idx, const Self::SpaceVectorType& volSize, const Types::Coordinate mmStep ) 
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
    case 6: 
    case 7:
    case 8:
      if ( this->NumberDOFs == 603 )
	{ // special case: 6 DOFs rigid plus 3 DOFs shear, but no scale
	return 0;
	}
      else
	{
	if ( this->m_LogScaleFactors )
	  return log( 1 + 0.5 * mmStep / volSize.MaxValue() );
	else
	  return 0.5 * mmStep / volSize.MaxValue();
	}
    case 9: 
    case 10:
    case 11:
      return 0.5 * mmStep / volSize.MaxValue();
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
( const Units::Radians angle, const Self::SpaceVectorType& direction,
  const Types::Coordinate* center, Self::MatrixType *const accumulate )
{
  Self::SpaceVectorType unit( direction );

  Self::SpaceVectorType center3D;
  if ( center ) 
    center3D = Self::SpaceVectorType::FromPointer( center );
  else
    center3D = Self::SpaceVectorType::FromPointer( this->RetCenter() );

  if ( accumulate ) 
    {
    unit += center3D;
    unit *= *accumulate;
    center3D *= *accumulate;
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

  double x = unit[0];
  double y = unit[1];
  double z = unit[2];

  // make a normalized quaternion
  const double w = MathUtil::Cos(0.5*angle);
  const double f = MathUtil::Sin(0.5*angle)/sqrt(x*x+y*y+z*z);
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

  xlate = xlate.GetInverse();
  this->Matrix *= xlate;
  this->DecomposeMatrix();

  if ( accumulate ) 
    {
    *accumulate *= matrix;
    *accumulate *= xlate;
    }
}

AffineXform::SmartPtr&
AffineXform::GetInverse()
{
  if ( !InverseXform ) 
    {
    InverseXform = AffineXform::SmartPtr( this->MakeInverse() );
    } 
  else 
    {
    this->UpdateInverse();
    }
  
  return this->InverseXform;
}

const AffineXform::SmartPtr&
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
  
  return this->InverseXform;
}

void
AffineXform::UpdateInverse() const
{
  if ( InverseXform ) 
    {
    InverseXform->NumberDOFs = this->NumberDOFs;
    InverseXform->m_LogScaleFactors = this->m_LogScaleFactors;
    InverseXform->Matrix = this->Matrix.GetInverse();
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

} // namespace
