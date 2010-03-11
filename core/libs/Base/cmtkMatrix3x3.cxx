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

#include <cmtkMatrix3x3.h>

#include <string.h>
#include <math.h>

#include <cmtkConsole.h>
#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T>
Matrix3x3<T>::Matrix3x3()
{
  memset( Matrix, 0, sizeof( Matrix ) );
  Matrix[0][0] = Matrix[1][1] = Matrix[2][2] = 1.0;
}

template<class T>
Matrix3x3<T>::Matrix3x3( const Self& other )
{
  memcpy( this->Matrix, other.Matrix, sizeof( this->Matrix ) );
}

template<class T>
Matrix3x3<T>&
Matrix3x3<T>::Set( const T *const values )
{
  memcpy( this->Matrix, values, sizeof( this->Matrix ) );
  return *this;
}

template<class T>
Matrix3x3<T>&
Matrix3x3<T>::Fill( const T value )
{
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      this->Matrix[i][j] = value;
  return *this;
}

template<class T>
Matrix3x3<T>&
Matrix3x3<T>::Compose( const Types::Coordinate params[8] )
{
  Types::Coordinate alpha = MathUtil::DegToRad(params[2]);

  this->Matrix[0][0] = static_cast<T>(  cos( alpha ) * params[3] );
  this->Matrix[0][1] = static_cast<T>( -sin( alpha ) * params[3] );
  this->Matrix[0][2] = static_cast<T>( 0.0 );
  this->Matrix[1][0] = static_cast<T>(  sin( alpha ) * params[4] );
  this->Matrix[1][1] = static_cast<T>(  cos( alpha ) * params[4] );
  this->Matrix[1][2] = static_cast<T>( 0.0 );
  this->Matrix[2][0] = static_cast<T>( 0.0 );
  this->Matrix[2][1] = static_cast<T>( 0.0 );
  this->Matrix[2][2] = static_cast<T>( 1.0 );

  // generate shears
  Self shearMatrix;
  shearMatrix[0][1] = static_cast<T>( params[5] );
  *this *= shearMatrix;

  // transform rotation center
  Types::Coordinate cM[2] = 
    {
      params[6]*this->Matrix[0][0] + params[7]*this->Matrix[1][0],
      params[6]*this->Matrix[0][1] + params[7]*this->Matrix[1][1],
    };
  
  // set translations
  this->Matrix[2][0] = static_cast<T>( params[0] - cM[0] + params[6] );
  this->Matrix[2][1] = static_cast<T>( params[1] - cM[1] + params[7] );

  return *this;
}

/**\todo Currently, we cannot correctly recover the parameters of a
 * transformation that includes a mirroring. Rotation also gets messed up
 * in this case. A possible solution may involve using the determinant to
 * detect the mirror and then come up with something that respects this.
 * Shear is not even supported yet.
 */
template<class T>
bool
Matrix3x3<T>::Decompose
( Types::Coordinate params[8], const Types::Coordinate *center ) const
{
  // make a working copy of the matrix for step-by-step decomposition
  Types::Coordinate matrix[3][3];
  memcpy( matrix, Matrix, sizeof( matrix ) );

  // translation entries
  params[0] = matrix[2][0];
  params[1] = matrix[2][1];

  if ( center ) 
    {
    Types::Coordinate cM[2] = { center[0]*matrix[0][0] + center[1]*matrix[1][0], center[0]*matrix[0][1] + center[1]*matrix[1][1] };
  
    params[0] += cM[0] - center[0];
    params[1] += cM[1] - center[1];

    memcpy( params+6, center, 2*sizeof( Types::Coordinate ) );
    }
  else
    {
    memset( params+6, 0, 2*sizeof( Types::Coordinate ) );
    }

#define IGS_EPSILON 0.001
  for ( int i=0; i<2; ++i ) 
    {
    // scale
    params[3+i] = sqrt( MathUtil::Square( matrix[i][0] ) + MathUtil::Square( matrix[i][1] ) );

    // report error on singular matrices.
    if ( fabs(params[3+i]) < IGS_EPSILON ) 
      {
      StdErr <<"igsMatrxi3x3::Decompose encountered singular matrix.";
      return false;
      }
    }

  // rotation
  // first rotate about y axis
  double x2 = matrix[0][0] / params[3];
  double y2 = -matrix[0][1] / params[3];
    
  double dot = x2 * x2 + y2 * y2;
  double d1 = sqrt (dot);

  double cosTheta, sinTheta;
  if (d1 < IGS_EPSILON) 
    {
    cosTheta = 1.0;
    sinTheta = 0.0;
    }
  else 
    {
    cosTheta = x2 / d1;
    sinTheta = y2 / d1;
    }
    
  params[2] = static_cast<T>( MathUtil::RadToDeg( atan2 (sinTheta, cosTheta) ) );
    
  return true;
#undef IGS_EPSILON
}

template<class T>
Matrix3x3<T>
Matrix3x3<T>::GetTranspose() const
{
  Self transpose;
  for ( int i = 0; i < 3; ++i ) 
    {
    for ( int j = 0; j < 3; ++j )
      transpose[i][j] = this->Matrix[j][i];
    }
  return transpose;
}
  
template<class T>
Matrix3x3<T>&
Matrix3x3<T>::Invert2x2()
{
  Self inverse;
  
  T rowBuff[3];
  for ( int col = 0; col<3; ++col ) 
    {    
    int pivIdx = col;
    T pivAbs = fabs( this->Matrix[col][col] );

    for ( int row = col+1; row<2; ++row ) // 2 to exclude last row!
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
    
    for ( int c=0; c<3; ++c ) 
      {
      if (c>col )
	this->Matrix[col][c] /= this->Matrix[col][col];
      inverse[col][c] /= this->Matrix[col][col];
      }
    this->Matrix[col][col] = 1.0;
    
    for ( int row = 0; row<3; ++row ) 
      {
      if (row != col ) 
	{
	for ( int c=0; c<3; ++c ) 
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
Matrix3x3<T>&
Matrix3x3<T>::Invert3x3()
{
  Self inverse;
  
  T rowBuff[3];
  for ( int col = 0; col<3; ++col ) {

    int pivIdx = col;
    T pivAbs = fabs( this->Matrix[col][col] );

    for ( int row = col+1; row<3; ++row ) {
      T nextAbs = fabs( this->Matrix[row][col] );
      if (nextAbs > pivAbs ) {
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
    
    for ( int c=0; c<3; ++c ) 
      {
      if (c>col )
	this->Matrix[col][c] /= this->Matrix[col][col];
      inverse[col][c] /= this->Matrix[col][col];
      }
    this->Matrix[col][col] = 1.0;
    
    for ( int row = 0; row<3; ++row ) 
      {
      if (row != col ) 
	{
	for ( int c=0; c<3; ++c ) 
	  {
	  if ( c>col ) 
	    this->Matrix[row][c] -= this->Matrix[row][col] * this->Matrix[col][c];
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
Matrix3x3<T>& 
Matrix3x3<T>::operator*=( const Self& other )
{
  return (*this = ((*this) * other));
}

template<class T>
Matrix3x3<T>& 
Matrix3x3<T>::operator*=( const T scalar )
{
  for ( int j=0; j<3; ++j ) 
    {
    for ( int i=0; i<3; ++i ) 
      {
      this->Matrix[i][j] *= scalar;
      }
    }
  return *this;
}

template<class T>
Matrix3x3<T>
Matrix3x3<T>::operator*
( const Self& other ) const
{
  Self result( NULL );

  for ( int j=0; j<3; ++j ) 
    {
    for ( int i=0; i<3; ++i ) 
      {
      result[i][j] = 0;
      for ( int k=0; k<3; ++k )
	result[i][j] += this->Matrix[i][k] * other.Matrix[k][j];
      }
    }
  
  return result;
}

template<class T>
Matrix3x3<T>& 
Matrix3x3<T>::operator=( const Self& other )
{
  memcpy( this->Matrix, other.Matrix, sizeof( this->Matrix ) );
  return *this;
}

template<class T>
T
Matrix3x3<T>::FrobeniusNorm() const
{
  T norm = 0.0;
  for ( int i = 0; i < 3; ++i ) 
    {
    for ( int j = 0; j < 3; ++j )
      norm += MathUtil::Square( this->Matrix[i][j] );
    }
  return sqrt( norm );
}

template<class T>
void
Matrix3x3<T>::ComputeEigenvalues( T (&lambda)[3] ) const
{
  const double M11 = this->Matrix[0][0];
  const double M12 = this->Matrix[0][1];
  const double M13 = this->Matrix[0][2];
  const double M22 = this->Matrix[1][1];
  const double M23 = this->Matrix[1][2];
  const double M33 = this->Matrix[2][2];

// <begin copyright notice>
// ---------------------------------------------------------------------------
//
//                Copyright (c) 2000-2003 TargetJr Consortium
//               GE Corporate Research and Development (GE CRD)
//                             1 Research Circle
//                            Niskayuna, NY 12309
//                            All Rights Reserved
//              Reproduction rights limited as described below.
//
//      Permission to use, copy, modify, distribute, and sell this software
//      and its documentation for any purpose is hereby granted without fee,
//      provided that (i) the above copyright notice and this permission
//      notice appear in all copies of the software and related documentation,
//      (ii) the name TargetJr Consortium (represented by GE CRD), may not be
//      used in any advertising or publicity relating to the software without
//      the specific, prior written permission of GE CRD, and (iii) any
//      modifications are clearly marked and summarized in a change history
//      log.
//
//      THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
//      EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
//      WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
//      IN NO EVENT SHALL THE TARGETJR CONSORTIUM BE LIABLE FOR ANY SPECIAL,
//      INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND OR ANY
//      DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
//      WHETHER OR NOT ADVISED OF THE POSSIBILITY OF SUCH DAMAGES, OR ON
//      ANY THEORY OF LIABILITY ARISING OUT OF OR IN CONNECTION WITH THE
//      USE OR PERFORMANCE OF THIS SOFTWARE.
//
// ---------------------------------------------------------------------------
// <end copyright notice>

  // Characteristic eqtn |M - xI| = 0
  // x^3 + b x^2 + c x + d = 0
  const double b = -M11-M22-M33;
  const double c =  M11*M22 +M11*M33 +M22*M33  -M12*M12 -M13*M13 -M23*M23;
  const double d = M11*M23*M23 +M12*M12*M33 +M13*M13*M22 -2.0*M12*M13*M23 -M11*M22*M33;
  // Using a numerically tweaked version of the real cubic solver http://www.1728.com/cubic2.htm
  const double b_3 = b/3.0;
  const double f = b_3*b_3 -  c/3.0 ;
  const double g = b*c/6.0 - b_3*b_3*b_3 - 0.5*d;
  
  if (f == 0.0 && g == 0.0)
    {
    lambda[0] = lambda[1] = lambda[2] = static_cast<T>( - b_3 );
    return;
    }
  
  const double f3 = f*f*f;
  const double g2 = g*g;
  const double sqrt_f = -sqrt(f);
       
   // deal explicitly with repeated root and treat
   // complex conjugate roots as numerically inaccurate repeated roots.
   
   // first check we are not too numerically innacurate
//  assert((g2 - f3) / vnl_math_sqr(b*b*b) < 1e-8);  
   
  if (g2 >= f3)
    {
    if (g < 0.0)
      {
      lambda[0] = static_cast<T>( 2.0 * sqrt_f  - b_3 );
      lambda[1] = lambda[2] = static_cast<T>( - sqrt_f - b_3 );
      }
    else
      {
      lambda[0] = lambda[1] = static_cast<T>( sqrt_f  - b_3 );
      lambda[2] = static_cast<T>( -2.0 * sqrt_f - b_3 );
      }
    return;
    }
   
 
  const double sqrt_f3 = sqrt_f * sqrt_f * sqrt_f;
  const double k = acos(g / sqrt_f3) / 3.0;
  const double j = 2.0 * sqrt_f;
  lambda[0] = static_cast<T>( j * cos(k) - b_3 );
  lambda[1] = static_cast<T>( j * cos(k + M_PI * 2.0 / 3.0) - b_3 );
  lambda[2] = static_cast<T>( j * cos(k - M_PI * 2.0 / 3.0) - b_3 );

  T tmp;
  if (lambda[1] < lambda[0]) 
    {
    tmp = lambda[1]; lambda[1] = lambda[0]; lambda[0] = tmp;
    }

  if (lambda[2] < lambda[1])
    {
    tmp = lambda[1]; lambda[1] = lambda[2]; lambda[2] = tmp;
    if (lambda[1] < lambda[0]) 
      {
      tmp = lambda[1]; lambda[1] = lambda[0]; lambda[0] = tmp;
      }
    }
}

template class Matrix3x3<float>;
template class Matrix3x3<double>;

} // namespace cmtk
