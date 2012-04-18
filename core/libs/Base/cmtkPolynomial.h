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

#ifndef __cmtkPolynomial_h_included_
#define __cmtkPolynomial_h_included_

#include <cmtkconfig.h>

#include <stdlib.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/** Generic class template for polynomials of arbitrary degrees.
 * This must be implemented by partial specialization for each degree that is
 * needed in the program.
 */
template<unsigned int NDegree,class TRealType=double>
class Polynomial
{
public:
  /// This class.
  typedef Polynomial<NDegree,TRealType> Self;

  /// Real value type.
  typedef TRealType RealValueType;

  /// Number of monomials in x, y, and z of degree up to NDegree.
  enum { NumberOfMonomials = 0 };

  /// Evaluate the idx'th monomial at (x,y,z).
  static TRealType EvaluateMonomialAt( const size_t idx, const TRealType x, const TRealType y, const TRealType z )
  {
    return 0.0;
  }

  /** Evaluate all monomials at one point.
   * This is more efficient than calling EvaluateMonomialAt() repeatedly, because the
   * computation can proceed incrementally and save most multiplications in the process.
   */
  static void EvaluateAllMonomials( TRealType *const mvec, const TRealType x, const TRealType y, const TRealType z )
  {
  }
};

/// Generic class template for polynomials of degree 1.
template<class TRealType>
class Polynomial<0,TRealType>
{
public:
  /// This class.
  typedef Polynomial<0,TRealType> Self;

  /// Real value type.
  typedef TRealType RealValueType;

  /// Number of monomials in x, y, and z of degree up to 1.
  enum { NumberOfMonomials = 0 };

  /// Evaluate the idx'th monomial at (x,y,z).
  static TRealType EvaluateMonomialAt( const size_t idx, const TRealType x, const TRealType y, const TRealType z )
  {
    return 0.0;
  }

  /** Evaluate all monomials at one point.
   * This is more efficient than calling EvaluateMonomialAt() repeatedly, because the
   * computation can proceed incrementally and save most multiplications in the process.
   */
  static void EvaluateAllMonomials( TRealType *const, const TRealType, const TRealType, const TRealType ) {}
};

/// Generic class template for polynomials of degree 1.
template<class TRealType>
class Polynomial<1,TRealType>
{
public:
  /// This class.
  typedef Polynomial<1,TRealType> Self;

  /// Real value type.
  typedef TRealType RealValueType;

  /// Number of monomials in x, y, and z of degree up to 1.
  enum { NumberOfMonomials = 3 };

  /// Evaluate the idx'th monomial at (x,y,z).
  static TRealType EvaluateMonomialAt( const size_t idx, const TRealType x, const TRealType y, const TRealType z )
  {
    switch ( idx )
      {
      case 0 : return  x;
      case 1 : return  y;
      case 2 : return  z;
      }
    return 0.0;
  }

  /** Evaluate all monomials at one point.
   * This is more efficient than calling EvaluateMonomialAt() repeatedly, because the
   * computation can proceed incrementally and save most multiplications in the process.
   */
  static void EvaluateAllMonomials( TRealType *const mvec, const TRealType x, const TRealType y, const TRealType z )
  {
    mvec[0] = x;
    mvec[1] = y;
    mvec[2] = z;
  }
};

/// Generic class template for polynomials of degree 2.
template<class TRealType>
class Polynomial<2,TRealType>
{
public:
  /// This class.
  typedef Polynomial<2,TRealType> Self;

  /// Real value type.
  typedef TRealType RealValueType;

  /// Number of monomials in x, y, and z of degree up to 2.
  enum { NumberOfMonomials = 9 };

  /// Evaluate the idx'th monomial at (x,y,z).
  static TRealType EvaluateMonomialAt( const size_t idx, const TRealType x, const TRealType y, const TRealType z )
  {
    switch ( idx )
      {
      case 0 : return  x;
      case 1 : return  y;
      case 2 : return  z;

      case 3 : return  x*x;
      case 4 : return  x*y;
      case 5 : return  x*z;
      case 6 : return  y*y;
      case 7 : return  y*z;
      case 8 : return  z*z;
      }
    return 0.0;
  }

  /** Evaluate all monomials at one point.
   * This is more efficient than calling EvaluateMonomialAt() repeatedly, because the
   * computation can proceed incrementally and save most multiplications in the process.
   */
  static void EvaluateAllMonomials( TRealType *const mvec, const TRealType x, const TRealType y, const TRealType z )
  {
    Polynomial<1,TRealType>::EvaluateAllMonomials( mvec, x, y, z );
    mvec[3] = mvec[0]*x;
    mvec[4] = mvec[0]*y;
    mvec[5] = mvec[0]*z;
    mvec[6] = mvec[1]*y;
    mvec[7] = mvec[1]*z;
    mvec[8] = mvec[2]*z;
  }
};

/// Generic class template for polynomials of degree 3.
template<class TRealType>
class Polynomial<3,TRealType>
{
public:
  /// This class.
  typedef Polynomial<3,TRealType> Self;

  /// Real value type.
  typedef TRealType RealValueType;

  /// Number of monomials in x, y, and z of degree up to 3.
  enum { NumberOfMonomials = 19 };

  /// Evaluate the idx'th monomial at (x,y,z).
  static TRealType EvaluateMonomialAt( const size_t idx, const TRealType x, const TRealType y, const TRealType z )
  {
    switch ( idx )
      {
      case  0 : return  x;
      case  1 : return  y;
      case  2 : return  z;

      case  3 : return  x*x;
      case  4 : return  x*y;
      case  5 : return  x*z;
      case  6 : return  y*y;
      case  7 : return  y*z;
      case  8 : return  z*z;

      case  9 : return  x*x*x;
      case 10 : return  x*x*y;
      case 11 : return  x*x*z;
      case 12 : return  x*y*y;
      case 13 : return  x*y*z;
      case 14 : return  x*z*z;
      case 15 : return  y*y*y;
      case 16 : return  y*y*z;
      case 17 : return  y*z*z;
      case 18 : return  z*z*z;
      }
    return 0.0;
  }

  /** Evaluate all monomials at one point.
   * This is more efficient than calling EvaluateMonomialAt() repeatedly, because the
   * computation can proceed incrementally and save most multiplications in the process.
   */
  static void EvaluateAllMonomials( TRealType *const mvec, const TRealType x, const TRealType y, const TRealType z )
  {
    Polynomial<2,TRealType>::EvaluateAllMonomials( mvec, x, y, z );

    mvec[9]  = mvec[3]*x;
    mvec[10] = mvec[3]*y;
    mvec[11] = mvec[3]*z;

    mvec[12] = mvec[4]*y;
    mvec[13] = mvec[4]*z;

    mvec[14] = mvec[5]*z;

    mvec[15] = mvec[6]*y;
    mvec[16] = mvec[6]*z;

    mvec[17] = mvec[7]*z;
    mvec[18] = mvec[8]*z;
  }
};

/// Generic class template for polynomials of degree 4.
template<class TRealType>
class Polynomial<4,TRealType>
{
public:
  /// This class.
  typedef Polynomial<4,TRealType> Self;

  /// Real value type.
  typedef TRealType RealValueType;

  /// Number of monomials in x, y, and z of degree up to 4.
  enum { NumberOfMonomials = 34 };

  /// Evaluate the idx'th monomial at (x,y,z).
  static TRealType EvaluateMonomialAt( const size_t idx, const TRealType x, const TRealType y, const TRealType z )
  {
    switch ( idx )
      {
      case  0 : return  x;
      case  1 : return  y;
      case  2 : return  z;

      case  3 : return  x*x;
      case  4 : return  x*y;
      case  5 : return  x*z;
      case  6 : return  y*y;
      case  7 : return  y*z;
      case  8 : return  z*z;

      case  9 : return  x*x*x;
      case 10 : return  x*x*y;
      case 11 : return  x*x*z;
      case 12 : return  x*y*y;
      case 13 : return  x*y*z;
      case 14 : return  x*z*z;
      case 15 : return  y*y*y;
      case 16 : return  y*y*z;
      case 17 : return  y*z*z;
      case 18 : return  z*z*z;

      case 19 : return  x*x*x*x;
      case 20 : return  x*x*x*y;
      case 21 : return  x*x*x*z;
      case 22 : return  x*x*y*y;
      case 23 : return  x*x*y*z;
      case 24 : return  x*x*z*z;
      case 25 : return  x*y*y*y;
      case 26 : return  x*y*y*z;
      case 27 : return  x*y*z*z;
      case 28 : return  x*z*z*z;

      case 29 : return  y*y*y*y;
      case 30 : return  y*y*y*z;
      case 31 : return  y*y*z*z;
      case 32 : return  y*z*z*z;
      case 33 : return  z*z*z*z;
      }
    return 0.0;
  }

  /** Evaluate all monomials at one point.
   * This is more efficient than calling EvaluateMonomialAt() repeatedly, because the
   * computation can proceed incrementally and save most multiplications in the process.
   */
  static void EvaluateAllMonomials( TRealType *const mvec, const TRealType x, const TRealType y, const TRealType z )
  {
    Polynomial<3,TRealType>::EvaluateAllMonomials( mvec, x, y, z );

    mvec[19] = mvec[9] * x;
    mvec[20] = mvec[9] * y;
    mvec[21] = mvec[9] * z;
    mvec[22] = mvec[10] * y;
    mvec[23] = mvec[10] * z;
    mvec[24] = mvec[11] * z;
    mvec[25] = mvec[12] * y;
    mvec[26] = mvec[12] * z;
    mvec[27] = mvec[13] * z;
    mvec[28] = mvec[14] * z;
    mvec[29] = mvec[15] * y;
    mvec[30] = mvec[15] * z;
    mvec[31] = mvec[16] * z;
    mvec[32] = mvec[17] * z;
    mvec[33] = mvec[18] * z;
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkPolynomial_h_included_

