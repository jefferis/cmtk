/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkMatrix4x4.h"

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkMatrix.h>
#include <Base/cmtkQRDecomposition.h>

#include <System/cmtkConsole.h>

#include <string.h>
#include <math.h>
#include <vector>

namespace
cmtk
{

template<class T>
Matrix4x4<T>::Matrix4x4( const Matrix3x3<T>& other )
{
  for ( int j=0; j<3; ++j ) 
    {
    for ( int i=0; i<3; ++i ) 
      {
      this->m_Matrix[i][j] = other[i][j];
      }
    }

  for ( int j=0; j<3; ++j ) 
    {
    this->m_Matrix[3][j] = this->m_Matrix[j][3] = 0.0;
    }
  this->m_Matrix[3][3] = 1.0;
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

  Self rotation = Self::Identity();
  rotation[0][0] = static_cast<T>( cos1*cos2 );
  rotation[0][1] = static_cast<T>( -cos1*sin2 );                     
  rotation[0][2] = static_cast<T>( -sin1 );
  rotation[1][0] = static_cast<T>(  (sin0xsin1*cos2 + cos0*sin2) );
  rotation[1][1] = static_cast<T>( (-sin0xsin1*sin2 + cos0*cos2) ); 
  rotation[1][2] = static_cast<T>(  sin0*cos1 );
  rotation[2][0] = static_cast<T>(  (cos0xsin1*cos2 - sin0*sin2) );
  rotation[2][1] = static_cast<T>( (-cos0xsin1*sin2 - sin0*cos2) );
  rotation[2][2] = static_cast<T>(  cos0*cos1 );

  // generate shears
  Self scaleShear = Self::Identity();
  for ( int i = 0; i < 3; ++i )
    {    
    scaleShear[i][i] = (logScaleFactors) ? exp( params[6+i] ) : params[6+i];
    scaleShear[(i/2)+(i%2)+1][i/2] = params[9+i];
    }
  *this = scaleShear * rotation;
  
  // transform rotation center
  const Types::Coordinate cM[3] = 
    {
      params[12]*this->m_Matrix[0][0] + params[13]*this->m_Matrix[1][0] + params[14]*this->m_Matrix[2][0],
      params[12]*this->m_Matrix[0][1] + params[13]*this->m_Matrix[1][1] + params[14]*this->m_Matrix[2][1],
      params[12]*this->m_Matrix[0][2] + params[13]*this->m_Matrix[1][2] + params[14]*this->m_Matrix[2][2]
    };
  
  // set translations
  this->m_Matrix[3][0] = params[0] - cM[0] + params[12];
  this->m_Matrix[3][1] = params[1] - cM[1] + params[13];
  this->m_Matrix[3][2] = params[2] - cM[2] + params[14];
  
  return *this;
}

template<class T>
bool
Matrix4x4<T>::Decompose
( Types::Coordinate params[15], const Types::Coordinate *center, const bool logScaleFactor ) const
{
  // translation entries
  params[0] = this->m_Matrix[3][0];
  params[1] = this->m_Matrix[3][1];
  params[2] = this->m_Matrix[3][2];

  if ( center )
    {
    const Types::Coordinate cM[3] = 
      {
	center[0]*this->m_Matrix[0][0] + center[1]*this->m_Matrix[1][0] + center[2]*this->m_Matrix[2][0],
	center[0]*this->m_Matrix[0][1] + center[1]*this->m_Matrix[1][1] + center[2]*this->m_Matrix[2][1],
	center[0]*this->m_Matrix[0][2] + center[1]*this->m_Matrix[1][2] + center[2]*this->m_Matrix[2][2],
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

  // use QR decomposition to separate rotation (Q) from shear and scales (R)
  Matrix2D<T> matrix2d( 3, 3 );
  for ( int i = 0; i < 3; ++i )
    {
    for ( int j = 0; j < 3; ++j )
      {
      // We use QR decomposition, but because we do left-multiplication of coordinate vector times transformation matrix,
      // the matrix is actually (RQ)^T, and we instead need to decompose the transformation matrix transpose.
      matrix2d[i][j] = this->m_Matrix[j][i];
      }
    }

  QRDecomposition<T> qr( matrix2d );
  const Matrix2D<T> R = qr.GetR();
  const Matrix2D<T> Q = qr.GetQ();

  for ( int k=0; k<3; ++k ) 
    {
    // if scale is negative, make positive and correct Q and R accordingly (we will figure out later if the overall transformation is a true rotation or has a negative determinant)
    if ( R[k][k] < 0 )
      {
      for ( int i=0; i<3; ++i ) 
	{
	R[k][i] = -R[k][i];
	Q[i][k] = -Q[i][k];
	}
      }

    // scale
    params[6 + k] = R[k][k];
    
    // report error on singular matrices.
    if ( params[6+k]  < std::numeric_limits<T>::epsilon() ) 
      {
      throw typename Self::SingularMatrixException();
      }

    // shear
    const int i = k / 2;           // i.e. i := { 0, 0, 1 }
    const int j = i + (k%2) + 1;   // i.e. j := { 1, 2, 2 } -- so i,j index the upper triangle of aMat, which is R from QR
    params[9+k] = R[i][j];
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
    
  const Types::Coordinate determinant = 
    this->m_Matrix[0][0]*this->m_Matrix[1][1]*this->m_Matrix[2][2] + 
    this->m_Matrix[0][1]*this->m_Matrix[1][2]*this->m_Matrix[2][0] + 
    this->m_Matrix[0][2]*this->m_Matrix[1][0]*this->m_Matrix[2][1] - 
    this->m_Matrix[0][2]*this->m_Matrix[1][1]*this->m_Matrix[2][0] - 
    this->m_Matrix[0][0]*this->m_Matrix[1][2]*this->m_Matrix[2][1] - 
    this->m_Matrix[0][1]*this->m_Matrix[1][0]*this->m_Matrix[2][2];

  // if negative determinant, this is not a true rotation
  if ( determinant < 0 )
    {
    // negate x scale
    params[6] = -params[6];
    // also negative shears related to x
    params[9] = -params[9];
    params[10] = -params[10];
    }
  
  // Now deal with the rotation Q:
  //  First rotate about y axis
  x2 = Q[1][0] / params[6];
  y2 = Q[2][0] / params[6];
  z2 = Q[0][0] / params[6];
    
  x3 = Q[1][2] / params[8];
  y3 = Q[2][2] / params[8];
  z3 = Q[0][2] / params[8];
    
  dot = x2 * x2 + z2 * z2;
  d1 = sqrt (dot);
    
//  if (d1 < std::numeric_limits<T>::epsilon()) 
  if (d1 < 1e-8) 
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
    
//  if (d < std::numeric_limits<T>::epsilon()) 
  if (d < 1e-8) 
    {    
    sinPhi = 0.0;
    cosPhi = 1.0;
    } 
  else 
    if (d1 < std::numeric_limits<T>::epsilon()) 
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
//  if (d2 < std::numeric_limits<T>::epsilon()) 
  if (d2 < 1e-8) 
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
Matrix4x4<T>::ChangeCoordinateSystem
( const FixedVector<3,T>& newX, const FixedVector<3,T>& newY )
{
  // rotate x axis to match new coordinate system
  Self rotate = RotateX( -atan2( newX[1], newX[2] ) );
  rotate *= RotateY( acos( newX[0] ) );

  // rotate previously rotated y axis further to match new coordinate system
  const FixedVector<3,T> newYrot = newY * rotate;
  rotate *= RotateX( atan2( newYrot[2], newYrot[1] ) );

  // z axis now matches automatically

  // apply rotation matrix to previous transformation to apply change of
  // coordinate systems.
  *this *= rotate;
  *this = rotate.GetInverse() * *this;

  return *this;
}


template<class T>
Matrix4x4<T>
Matrix4x4<T>::RotateX( const T angle )
{
  Self rot = Self::Identity();
  rot[1][1] = rot[2][2] = cos( angle );
  rot[1][2] = -1.0 * (rot[2][1] = sin( angle ) );

  return rot;
}
  
template<class T>
Matrix4x4<T> 
Matrix4x4<T>::RotateY( const T angle )
{
  Self rot = Self::Identity();
  rot[0][0] = rot[2][2] = cos( angle );
  rot[0][2] = -1.0 * (rot[2][0] = sin( angle ) );

  return rot;
}
  
template<class T>
Matrix4x4<T> 
Matrix4x4<T>::RotateZ( const T angle )
{
  Self rot = Self::Identity();
  rot[0][0] = rot[1][1] = cos( angle );
  rot[0][1] = -1.0 * (rot[1][0] = sin( angle ) );

  return rot;
}

template class Matrix4x4<Types::Coordinate>;

} // namespace cmtk
