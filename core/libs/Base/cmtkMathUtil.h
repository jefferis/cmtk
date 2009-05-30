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

#ifndef __cmtkMathUtil_h_included_
#define __cmtkMathUtil_h_included_

#include <cmtkconfig.h>

#include <cmtkArray.h>
#include <cmtkMatrix.h>

#include <math.h>
#include <stdlib.h>

#ifdef HAVE_STDINT_H
#  include <stdint.h>
#else
#  ifdef _MSC_VER
typedef unsigned int uint32_t;
#  endif // _MSC_VER
#endif

#include <algorithm>

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

#ifdef _MSC_VER
#include <float.h>

/* Some useful constants taken from SGI's math.h */
#define M_E             2.7182818284590452354
#define M_LOG2E         1.4426950408889634074
#define M_LOG10E        0.43429448190325182765
#define M_LN2           0.69314718055994530942
#define M_LN10          2.30258509299404568402
#define M_PI_2          1.57079632679489661923
#define M_PI_4          0.78539816339744830962
#define M_1_PI          0.31830988618379067154
#define M_2_PI          0.63661977236758134308
#define M_2_SQRTPI      1.12837916709551257390
#define M_SQRT2         1.41421356237309504880
#define M_SQRT1_2       0.70710678118654752440
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{


/** General-purpose mathematical functions and function templates.
 */
namespace MathUtil 
{

/// Union that is used for initializing FP values with special values.
typedef 
union FPInitializeUnion
{
  /// Bit pattern for special values.
  uint32_t i[2];
  /// Read out for single-precision special value.
  float f[2];
  /// Read out for double-precision special value.
  double d;
} FPInitializeUnion;

/// Constant to generate NaN values.
extern const FPInitializeUnion FPInitializeNaN;
/// Constant to generate Inf values.
extern const FPInitializeUnion FPInitializeInf;

/// Pointers to constants to generate double-precision NaN values.
extern const void* FPInitializeNaN_P;
/// Constant to generate Inf values.
extern const void* FPInitializeInf_P;

/// Pointers to constants to generate single-precision NaN values.
extern const void* FPInitializeNaN_fP;
/// Constant to generate Inf values.
extern const void* FPInitializeInf_fP;

/// Double-precision not-a-number (NaN) value.
#define CMTK_DOUBLE_NAN cmtk::MathUtil::FPInitializeNaN.d
#define CMTK_DOUBLE_NAN_P cmtk::MathUtil::FPInitializeNaN_P
#define CMTK_FLOAT_NAN_P cmtk::MathUtil::FPInitializeNaN_fP

/// Single-precision not-a-number (NaN) value.
#if WORDS_BIGENDIAN
#  define CMTK_FLOAT_NAN cmtk::MathUtil::FPInitializeNaN.f[0]
#else
#  define CMTK_FLOAT_NAN cmtk::MathUtil::FPInitializeNaN.f[1]
#endif

/// Double-precision infinite (Inf) value.
#define CMTK_DOUBLE_INF cmtk::MathUtil::FPInitializeInf.d
#define CMTK_DOUBLE_INF_P cmtk::MathUtil::FPInitializeInf_P
#define CMTK_FLOAT_INF_P cmtk::MathUtil::FPInitializeInf_fP

/// Single-precision infinite (Inf) value.
#if WORDS_BIGENDIAN
#  define CMTK_FLOAT_INF cmtk::MathUtil::FPInitializeInf.f[0]
#else
#  define CMTK_FLOAT_INF cmtk::MathUtil::FPInitializeInf.f[1]
#endif

#ifndef TINY
#define TINY 1e-8
#endif

/// Return square of a float value.
template<class T> T Square ( const T a) { return a*a; }

/** Return minimum of an array of ordered values.
   */
template<class T> T Min ( const int nValues, const T* Values )
{
  T Result = Values[0];
  for ( int idx=1; idx<nValues; ++idx )
    Result = std::min( Result, Values[idx] );
  
  return Result;
}

/** Return minimum of an array of ordered values.
   */
template<class T> T Max ( const int nValues, const T* Values )
{
  T Result = Values[0];
  for ( int idx=1; idx<nValues; ++idx )
    Result = std::max( Result, Values[idx] );
  
  return Result;
}

/// Swap two primitive values.
template<class T> void Swap( T& a, T& b ) 
{
  T temp = a;
  a = b;
  b = temp;
}

/// Return length of intersection of two intervals.
template<class T> T Intersect ( const T aMin, const T aMax, const T bMin, const T bMax )
{
  return ( std::min( aMax, bMax ) - std::max( aMin, bMin ) );
}

/// Round float value to the nearest integer.
template<class T> int Round ( const T a ) { return (int)floor(a+0.5); }

/** Return sign of float value.
   *@return -1 if a<0, 1 if a>0, 0 if a==0.
   */
template<class T> int Sign ( const T a ) { return (a<0)?-1:((a==0)?0:1); }

/** Check if some float value is within a range.
   *@return 0 if value is in range, -1 if value is below minumum, 1 if value
   * is above maximum.
   */ 
template<class T> int CheckRange ( const T value, const T a, const T b ) 
{
  int sigA = Sign(value-a);
  
  return (sigA == Sign(value-b))?sigA:0;
}

template<class T> int TruncAndFrac
( const T value, const T ofs, const T delta, T& weight )
{ 
  T x = (value-ofs) / delta; 
  int t = static_cast<int>( floor( x ) );
  weight = 1.0 - (x-t); 
  return t;
}

/// Compute p*log(p) for a single value.
template<class T> double plogp( const T p ) { return (p>0) ? p * log( p ) : 0.0; }

/** Computes average of an array of float values.
   */
template<class T>
T Mean
( const unsigned int nValues, const T* values );

/** Computes variance of an array of float values.
    *\param unbiased If this flag is set (default: unset), then the variance
    * will be computed over nValues-1; otherwise over nValues.
   */
template<class T>
T Variance
( const unsigned int nValues, const T* values, T mean, const bool unbiased = false );


/** Normalized correlation coefficient between two float vectors.
   */
template<class T>
T Correlation( const size_t n, const T* x, const T* y );

/// Compute t-statistic from coefficient of correlation.
double TStatFromCorrelation
( const double r, //!< Coefficient of correlation as computed by MathUtil::Correlation function.
  const size_t df ); //!< Number of degrees of freedom

/// Compute probability from T-statistic.
double ProbabilityFromTStat
( const double t, //!< T-statistic as returned for example from MathUtil::TStatFromCorrelation function.
  const size_t df ); //!< Number of degrees of freedom.

/** Performs t-test on two distributions.
   */
template<class T>
T TTest ( const unsigned int nValuesX, const T* valuesX, const unsigned int nValuesY, const T* valuesY, T& t );

/** Performs t-test on two distributions.
   * Also return average value for each distribution.
   */
template<class T>
T TTest ( const unsigned int nValuesX, const T* valuesX, const unsigned int nValuesY, const T* valuesY, T& t, T& avgX, T& avgY );

/** Performs t-test on two distributions.
   * Also return average value for each distribution.
   */
template<class T>
T PairedTTest ( const unsigned int nValues, const T* valuesX, const T* valuesY, T& t, T& avgX, T& avgY );

/** Performs one-sample t-test on distribution to test for zero mean.
   * Also return average value for each distribution.
   */
template<class T>
T TTest ( const unsigned int nValuesX, const T* valuesX, T& t, T& avgX );

/// Beta-i function.
double Betai( const double a, const double b, const double x );

/// Beta-Cf function.
double BetaCf( const double a, const double b, const double x );

/// GammaLn function.
double GammaLn( const double xx );

/// Singular Value Decomposition
void SVD( double **U, size_t m, size_t n, double *W, double **V );

/** Convert degrees to radians.
   * No range checking is done. If the given angle is outside the range 
   * [0..360], the returned value will be outside the range [0..2Pi]. However,
   * the result will still be correct in the sense that it will denote the
   * same angle.
   *@param deg Angle in degrees (0 through 360).
   *@return Equivalent angle in radians (0 through 2 Pi).
   */
inline double DegToRad ( const double deg ) 
{ return static_cast<double>( deg * (M_PI/180.0) ); }

/** Convert radians to degrees.
   * No range checking is done. If the given angle is outside the range 
   * [0..2 Pi], the returned value will be outside the range [0..2Pi]. However,
   * the result will still be correct in the sense that it will denote the
   * same angle.
   *@param rad Angle in radians (0 through 2 Pi).
   *@return Equivalent angle in degrees (0 through 360).
   */
inline double RadToDeg ( const double rad ) 
{ return static_cast<float>( rad * (180.0/M_PI) ); }

/// Function that compares two floats; to be used in qsort().
inline int CompareFloat( const void *a, const void *b ) 
{
  const float* A = static_cast<const float*>( a );
  const float* B = static_cast<const float*>( b );
  if ( *A > *B ) return +1;
  if ( *A < *B ) return -11;
  return 0;
}

/// Function that compares two doubles; to be used in qsort().
inline int CompareDouble( const void *a, const void *b ) 
{
  const double* A = static_cast<const double*>( a );
  const double* B = static_cast<const double*>( b );
  if ( *A > *B ) return +1;
  if ( *A < *B ) return -11;
  return 0;
}

/** Generate normally distributed random numbers.
   * This function uses the Box-Muller method to transform a pair of uniformly
   * distributed random numbers into a pair of normally (ie., Gaussian) 
   * distributed random numbers. One of the two generated numbers is returned
   * while the other is stored so that, when this function is called the next 
   * time, the previously computed value can be returned without further 
   * computational expense.
   *@param sigma Standard deviation of the resulting distribution.
   */
inline double NormalRandom( const double sigma ) {
  static bool secondNumberReady = false;
  static double secondNumber = 0;

  if ( secondNumberReady ) {
  secondNumberReady = false;
  return secondNumber;
  }

  double x1, x2, w;
  do {
  x1 = 2.0 * (random()&0xffffff)/0x1000000 - 1.0;
  x2 = 2.0 * (random()&0xffffff)/0x1000000 - 1.0;
  w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
    
  w = sqrt( (-2.0 * log( w ) ) / w );
  secondNumber = x1 * w * sigma;
  return x2 * w * sigma;
}

/** Generate normally distributed random numbers with explicit seed.
   *@param seed Random seed given to srandom() function.
   */
inline double NormalRandom( const double sigma, const unsigned int seed ) {
  srandom( seed );
  return NormalRandom( sigma );
}

/// Uniform random number generator.
double UniformRandom();

/// Compute eigensystem and eigenvalues for square real matrix using Jacobi rotation.
template<class T> void ComputeEigensystem( const Matrix2D<T>& matrix, Matrix2D<T>& eigensystem, Array<T>& eigenvalues );

/// Compute eigenvalues for square real matrix using Jacobi rotation.
template<class T> void ComputeEigenvalues( const Matrix2D<T>& matrix, Array<T>& eigenvalues );

/// Determinant of an n x n square matrix.
template<class T> T CholeskyDeterminant( const Matrix2D<T>& matrix, int n);

} // namespace MathUtil

//@}

} // namespace cmtk

#include <cmtkMathUtilStatistics.txx>

#endif // #ifndef __cmtkMathUtil_h_included_
