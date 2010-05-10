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

#include <cmtkMatrix4x4.h>

#include <cmtkConsole.h>
#include <cmtkMathUtil.h>
#include <cmtkMatrix.h>
#include <cmtkQRDecomposition.h>

#include <string.h>
#include <math.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T>
Matrix4x4<T>::Matrix4x4()
{
  memset( Matrix, 0, sizeof( Matrix ) );
  Matrix[0][0] = Matrix[1][1] = Matrix[2][2] = Matrix[3][3] = 1.0;
}

template<class T>
Matrix4x4<T>::Matrix4x4( const Self& other )
{
  memcpy( this->Matrix, other.Matrix, sizeof( this->Matrix ) );
}

template<class T>
Matrix4x4<T>::Matrix4x4( const Matrix3x3<T>& other )
{
  for ( int j=0; j<3; ++j ) 
    {
    for ( int i=0; i<3; ++i ) 
      {
      this->Matrix[i][j] = other[i][j];
      }
    }

  for ( int j=0; j<3; ++j ) 
    {
    this->Matrix[3][j] = this->Matrix[j][3] = 0.0;
    }
  this->Matrix[3][3] = 1.0;
}

template<class T>
Matrix4x4<T>&
Matrix4x4<T>::Set( const T *const values )
{
  memcpy( this->Matrix, values, sizeof( this->Matrix ) );
  return *this;
}

template<class T>
Matrix4x4<T>
Matrix4x4<T>::GetTranspose() const
{
  Self transpose;
  for ( int i = 0; i < 4; ++i ) 
    {
    for ( int j = 0; j < 4; ++j )
      transpose[i][j] = this->Matrix[j][i];
    }
  return transpose;
}
  
template<class T>
Matrix4x4<T>&
Matrix4x4<T>::Compose
( const Types::Coordinate params[15], const bool logScaleFactors )
{
  const Units::Radians alpha = Units::Degrees( params[3] );
  const Units::Radians theta = Units::Degrees( params[4] );
  const Units::Radians   phi = Units::Degrees( params[5] );

  const double cos0 = MathUtil::Cos(alpha), sin0 = MathUtil::Sin(alpha);
  const double cos1 = MathUtil::Cos(theta), sin1 = MathUtil::Sin(theta);
  const double cos2 = MathUtil::Cos(  phi), sin2 = MathUtil::Sin(  phi);

  const double sin0xsin1 = sin0 * sin1;
  const double cos0xsin1 = cos0 * sin1;

  const double scaleX = (logScaleFactors) ? exp( params[6] ) : params[6];
  const double scaleY = (logScaleFactors) ? exp( params[7] ) : params[7];
  const double scaleZ = (logScaleFactors) ? exp( params[8] ) : params[8];

  Matrix[0][0] = static_cast<T>( cos1*cos2 * scaleX );
  Matrix[0][1] = static_cast<T>( -cos1*sin2 * scaleX );                     
  Matrix[0][2] = static_cast<T>( -sin1 * scaleX );
  Matrix[0][3] = static_cast<T>( 0 );
  Matrix[1][0] = static_cast<T>(  (sin0xsin1*cos2 + cos0*sin2) * scaleY );
  Matrix[1][1] = static_cast<T>( (-sin0xsin1*sin2 + cos0*cos2) * scaleY ); 
  Matrix[1][2] = static_cast<T>(  sin0*cos1 * scaleY );
  Matrix[1][3] = static_cast<T>( 0 );
  Matrix[2][0] = static_cast<T>(  (cos0xsin1*cos2 - sin0*sin2) * scaleZ );
  Matrix[2][1] = static_cast<T>( (-cos0xsin1*sin2 - sin0*cos2) * scaleZ );
  Matrix[2][2] = static_cast<T>(  cos0*cos1 * scaleZ );
  Matrix[2][3] = static_cast<T>( 0 );

  Matrix[3][0] = Matrix[3][1] = Matrix[3][2] = static_cast<T>( 0 );
  Matrix[3][3] = static_cast<T>( 1.0 );

  // generate shears
  for ( int i = 2; i >= 0; --i )
    {
    Self shear;
    shear[i/2][(i/2)+(i%2)+1] = params[9+i];
    *this *= shear;
    }
  
  // transform rotation center
  const Types::Coordinate cM[3] = 
    {
      params[12]*Matrix[0][0] + params[13]*Matrix[1][0] + params[14]*Matrix[2][0],
      params[12]*Matrix[0][1] + params[13]*Matrix[1][1] + params[14]*Matrix[2][1],
      params[12]*Matrix[0][2] + params[13]*Matrix[1][2] + params[14]*Matrix[2][2]
    };
  
  // set translations
  Matrix[3][0] = params[0] - cM[0] + params[12];
  Matrix[3][1] = params[1] - cM[1] + params[13];
  Matrix[3][2] = params[2] - cM[2] + params[14];
  
  return *this;
}

template<class T>
bool
Matrix4x4<T>::Decompose
( Types::Coordinate params[12], const Types::Coordinate *center, const bool logScaleFactor ) const
{
  // make a working copy of the matrix for step-by-step decomposition
  Self matrix( *this );

  // translation entries
  params[0] = matrix[3][0];
  params[1] = matrix[3][1];
  params[2] = matrix[3][2];

  if ( center )
    {
    const Types::Coordinate cM[3] = 
      {
	center[0]*matrix[0][0] + center[1]*matrix[1][0] + center[2]*matrix[2][0],
	center[0]*matrix[0][1] + center[1]*matrix[1][1] + center[2]*matrix[2][1],
	center[0]*matrix[0][2] + center[1]*matrix[1][2] + center[2]*matrix[2][2],
      };
    
    params[0] += cM[0] - center[0];
    params[1] += cM[1] - center[1];
    params[2] += cM[2] - center[2];
    
    if ( center != params+12)
      {
      // sometimes we may get a pointer to our own parameters
      memcpy( params+12, center, 3*sizeof( Types::Coordinate ) );
      }
    } 
  else
    {
    memset( params+12, 0, 3*sizeof( Types::Coordinate ) );
    }

  // use QR decomposition to get shears (gets the scales, too, but we'll discard those 
  // and do them explicitly later.
  Matrix2D<T> matrix2d( 3, 3 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      matrix2d[i][j] = matrix[i][j];

  std::vector<T> R_diagonal( 3 );

  QRDecomposition<T> qr( matrix2d );

  Matrix2D<T> R = *(qr.GetR());

  for ( int i = 0; i < 3; i++ )
    R_diagonal[i] = R[i][i];
                                                                   
  // shear
  for ( int k = 0; k<3; ++k ) 
    {
    const int i = k / 2;           // i.e. i := { 0, 0, 1 }
    const int j = i + (k%2) + 1;   // i.e. j := { 0, 1, 2 } -- so i,j index the upper triangle of aMat, which is R from QR
    params[9+k] = R[i][j] / R_diagonal[i];

    // remove contribution from transformation matrix
    Self shear;
    shear[i][j] = params[9+k];
    matrix *= shear.Invert();
    }
  
/*=========================================================================

THE FOLLOWING CODE WAS ADOPTED AND MODIFIED FROM VTK, The Visualization
Toolkit.

  Program:   Visualization Toolkit
  Language:  C++
  Thanks:    Thanks to David G. Gobbi who developed this class.

Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#define VTK_AXIS_EPSILON 0.001
  double   d;
  double   d1;
  double   d2;
  double   dot;
  double   cosPhi, sinPhi;
  double   cosTheta, sinTheta;
  double   cosAlpha, sinAlpha;
  double   x2, y2, z2;
  double   x3, y3, z3;
  double   x3p, y3p;
    
  for ( int i=0; i<3; ++i ) 
    {
    // scale
    params[6 + i] = sqrt( MathUtil::Square( matrix[i][0] ) + MathUtil::Square( matrix[i][1] ) + MathUtil::Square( matrix[i][2] ) );
    
    // report error on singular matrices.
    if ( fabs(params[6+i]) < VTK_AXIS_EPSILON ) 
      {
      StdErr <<"Matrxi4x4::Decompose encountered singular matrix.";
      return false;
      }
    }

  // if negative determinant, negate x scale
  const Types::Coordinate determinant = 
    this->Matrix[0][0]*this->Matrix[1][1]*this->Matrix[2][2] + 
    this->Matrix[0][1]*this->Matrix[1][2]*this->Matrix[2][0] + 
    this->Matrix[0][2]*this->Matrix[1][0]*this->Matrix[2][1] - 
    this->Matrix[0][2]*this->Matrix[1][1]*this->Matrix[2][0] - 
    this->Matrix[0][0]*this->Matrix[1][2]*this->Matrix[2][1] - 
    this->Matrix[0][1]*this->Matrix[1][0]*this->Matrix[2][2];

  if ( determinant < 0 )
    params[6] = -params[6];
  
  // rotation
  // first rotate about y axis
  x2 = matrix[0][1] / params[6];
  y2 = matrix[0][2] / params[6];
  z2 = matrix[0][0] / params[6];
    
  x3 = matrix[2][1] / params[8];
  y3 = matrix[2][2] / params[8];
  z3 = matrix[2][0] / params[8];
    
  dot = x2 * x2 + z2 * z2;
  d1 = sqrt (dot);
    
  if (d1 < VTK_AXIS_EPSILON) 
    {
    cosTheta = 1.0;
    sinTheta = 0.0;
    } 
  else 
    {
    cosTheta = z2 / d1;
    sinTheta = x2 / d1;
    }
  
  params[5] = Units::Degrees( -MathUtil::ArcTan2( sinTheta, cosTheta ) ).Value(); // theta
    
    // now rotate about x axis
  dot = x2 * x2 + y2 * y2 + z2 * z2;
  d = sqrt (dot);
    
  if (d < VTK_AXIS_EPSILON) 
    {    
    sinPhi = 0.0;
    cosPhi = 1.0;
    } 
  else 
    if (d1 < VTK_AXIS_EPSILON) 
      {
      sinPhi = y2 / d;
      cosPhi = z2 / d;
      } 
    else 
      {
      sinPhi = y2 / d;
      cosPhi = ( x2 * x2 + z2 * z2) / (d1 * d);
      }
  
  params[4] = Units::Degrees( -MathUtil::ArcTan2( sinPhi, cosPhi ) ).Value(); // phi 
  
  // finally, rotate about z
  x3p = x3 * cosTheta - z3 * sinTheta;
  y3p = - sinPhi * sinTheta * x3 + cosPhi * y3 - sinPhi * cosTheta * z3;
  dot = x3p * x3p + y3p * y3p;
  
  d2 = sqrt (dot);
  if (d2 < VTK_AXIS_EPSILON) 
    {
    cosAlpha = 1.0;
    sinAlpha = 0.0;
    } 
  else
    {
    cosAlpha = y3p / d2;
    sinAlpha = x3p / d2;
    }
  
  params[3] = Units::Degrees( -MathUtil::ArcTan2( sinAlpha, cosAlpha ) ).Value(); // alpha
  
  if ( logScaleFactor )
    {
    for ( int i = 6; i < 9; ++i )
      params[i] = log( params[i] );
    }
  
  return true;
  
  /** END OF ADOPTED VTK CODE **/
}

template<class T>
Matrix4x4<T>&
Matrix4x4<T>::Invert()
{
  Self inverse;
  
  T rowBuff[4];
  for ( int col = 0; col<4; ++col ) 
    {    
    int pivIdx = col;
    T pivAbs = fabs( this->Matrix[col][col] );

    for ( int row = col+1; row<3; ++row )  // 3 to exclude last row!
      {
      T nextAbs = fabs( this->Matrix[row][col] );
      if (nextAbs > pivAbs ) 
	{
	pivIdx = row;
	pivAbs = nextAbs;
	}
      }
    
    if ( col != pivIdx )
      {
      memcpy( rowBuff, this->Matrix[col], sizeof(rowBuff) );
      memcpy( this->Matrix[col], this->Matrix[pivIdx], sizeof(rowBuff) );
      memcpy( this->Matrix[pivIdx], rowBuff, sizeof(rowBuff) );
      
      memcpy( rowBuff, inverse[col],sizeof(rowBuff));
      memcpy( inverse[col], inverse[pivIdx], sizeof(rowBuff) );
      memcpy( inverse[pivIdx], rowBuff, sizeof(rowBuff) );
      }
    
    for ( int c=0; c<4; ++c ) 
      {
      if (c>col )
	this->Matrix[col][c] /= this->Matrix[col][col];
      inverse[col][c] /= this->Matrix[col][col];
      }
    this->Matrix[col][col] = 1.0;
    
    for ( int row = 0; row<4; ++row ) 
      {
      if (row != col ) 
	{
	for ( int c=0; c<4; ++c ) 
	  {
	  if ( c>col ) 
	    this->Matrix[row][c] -= 
	      this->Matrix[row][col] * this->Matrix[col][c];
	  inverse[row][c] -= this->Matrix[row][col] * inverse[col][c];
	  }
	this->Matrix[row][col] = 0;
	}
      }
    }
  
  // finally copy inverse into this object.
  return (*this = inverse);
}

template<class T>
Matrix4x4<T>& 
Matrix4x4<T>::operator*=( const Self& other )
{
  return (*this = ((*this) * other));
}

template<class T>
const Matrix4x4<T>
Matrix4x4<T>::operator*
( const Self& other ) const
{
  Self result( NULL );

  for ( int j=0; j<4; ++j ) 
    {
    for ( int i=0; i<4; ++i ) 
      {
      result[i][j] = 0;
      for ( int k=0; k<4; ++k )
	result[i][j] += this->Matrix[i][k] * other.Matrix[k][j];
      }
    }

  return result;
}

template<class T>
Matrix4x4<T>& 
Matrix4x4<T>::operator=( const Self& other )
{
  memcpy( this->Matrix, other.Matrix, sizeof( this->Matrix ) );
  return *this;
}

template<class T>
Matrix4x4<T>& 
Matrix4x4<T>::operator=( const Matrix3x3<T>& other )
{
  for ( int j=0; j<3; ++j ) 
    {
    for ( int i=0; i<3; ++i ) 
      {
      this->Matrix[i][j] = other[i][j];
      }
    }

  for ( int j=0; j<3; ++j ) 
    {
    this->Matrix[3][j] = this->Matrix[j][3] = 0.0;
    }
  this->Matrix[3][3] = 1.0;
  
  return *this;
}

template<class T>
Matrix4x4<T>& 
Matrix4x4<T>::ChangeCoordinateSystem
( const FixedVector<3,T>& newX, const FixedVector<3,T>& newY )
{
  // rotate x axis to match new coordinate system
  Self rotate = RotateX( -atan2( newX[1], newX[2] ) );
  rotate *= RotateY( acos( newX[0] ) );

  // rotate previously rotated y axis further to match new coordinate system
  FixedVector<3,T> newYrot;
  rotate.Multiply( newY, newYrot );
  rotate *= RotateX( atan2( newYrot[2], newYrot[1] ) );

  // z axis now matches automatically

  // apply rotation matrix to previous transformation to apply change of
  // coordinate systems.
  *this *= rotate;
  rotate.Invert();
  *this = rotate * *this;

  return *this;
}


template<class T>
Matrix4x4<T>
Matrix4x4<T>::RotateX( const T angle )
{
  Self rot;
  rot[1][1] = rot[2][2] = cos( angle );
  rot[1][2] = -1.0 * (rot[2][1] = sin( angle ) );

  return rot;
}
  
template<class T>
Matrix4x4<T> 
Matrix4x4<T>::RotateY( const T angle )
{
  Self rot;
  rot[0][0] = rot[2][2] = cos( angle );
  rot[0][2] = -1.0 * (rot[2][0] = sin( angle ) );

  return rot;
}
  
template<class T>
Matrix4x4<T> 
Matrix4x4<T>::RotateZ( const T angle )
{
  Self rot;
  rot[0][0] = rot[1][1] = cos( angle );
  rot[0][1] = -1.0 * (rot[1][0] = sin( angle ) );

  return rot;
}

template<class T>
T
Matrix4x4<T>::FrobeniusNorm() const
{
  T norm = 0.0;
  for ( int i = 0; i < 4; ++i ) 
    {
    for ( int j = 0; j < 4; ++j )
      norm += MathUtil::Square( this->Matrix[i][j] );
    }
  return sqrt( norm );
}

template class Matrix4x4<Types::Coordinate>;

} // namespace cmtk
