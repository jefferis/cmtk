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

#include "cmtkMatrix3x3.h"

#include <System/cmtkConsole.h>
#include <Base/cmtkMathUtil.h>

#include <string.h>
#include <math.h>
#include <limits>

namespace
cmtk
{

template<class T>
Matrix3x3<T>&
Matrix3x3<T>::Compose( const Types::Coordinate params[8] )
{
  const Units::Radians alpha = Units::Degrees( params[2] );

  (*this)[0][0] = static_cast<T>(  MathUtil::Cos( alpha ) * params[3] );
  (*this)[0][1] = static_cast<T>( -MathUtil::Sin( alpha ) * params[3] );
  (*this)[0][2] = static_cast<T>( 0.0 );
  (*this)[1][0] = static_cast<T>(  MathUtil::Sin( alpha ) * params[4] );
  (*this)[1][1] = static_cast<T>(  MathUtil::Cos( alpha ) * params[4] );
  (*this)[1][2] = static_cast<T>( 0.0 );
  (*this)[2][0] = static_cast<T>( 0.0 );
  (*this)[2][1] = static_cast<T>( 0.0 );
  (*this)[2][2] = static_cast<T>( 1.0 );

  // generate shears
  Self shearMatrix = Self::Identity();
  shearMatrix[0][1] = static_cast<T>( params[5] );
  *this *= shearMatrix;

  // transform rotation center
  Types::Coordinate cM[2] = 
    {
      params[6]*(*this)[0][0] + params[7]*(*this)[1][0],
      params[6]*(*this)[0][1] + params[7]*(*this)[1][1],
    };
  
  // set translations
  (*this)[2][0] = static_cast<T>( params[0] - cM[0] + params[6] );
  (*this)[2][1] = static_cast<T>( params[1] - cM[1] + params[7] );

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
  memcpy( matrix, this->m_Matrix, sizeof( matrix ) );

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

  for ( int i=0; i<2; ++i ) 
    {
    // scale
    params[3+i] = sqrt( MathUtil::Square( matrix[i][0] ) + MathUtil::Square( matrix[i][1] ) );

    // report error on singular matrices.
    if ( fabs(params[3+i]) < std::numeric_limits<T>::epsilon() ) 
      {
      throw typename Self::SingularMatrixException();
      }
    }

  // rotation
  // first rotate about y axis
  double x2 = matrix[0][0] / params[3];
  double y2 = -matrix[0][1] / params[3];
    
  double dot = x2 * x2 + y2 * y2;
  double d1 = sqrt (dot);

  double cosTheta, sinTheta;
//  if (d1 < std::numeric_limits<T>::epsilon() ) 
  if (d1 < 1e-8 ) 
    {
    cosTheta = 1.0;
    sinTheta = 0.0;
    }
  else 
    {
    cosTheta = x2 / d1;
    sinTheta = y2 / d1;
    }
    
  params[2] = static_cast<T>( Units::Degrees( MathUtil::ArcTan2 (sinTheta, cosTheta) ).Value() );
    
  return true;
}

template<class T>
void
Matrix3x3<T>::ComputeEigenvalues( T (&lambda)[3] ) const
{
  const double M11 = (*this)[0][0];
  const double M12 = (*this)[0][1];
  const double M13 = (*this)[0][2];
  const double M22 = (*this)[1][1];
  const double M23 = (*this)[1][2];
  const double M33 = (*this)[2][2];

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
