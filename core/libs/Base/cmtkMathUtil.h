/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <Base/cmtkUnits.h>
#include <Base/cmtkMatrix.h>

#ifdef HAVE_STDINT_H
#  include <stdint.h>
#else
#  ifdef _MSC_VER
typedef unsigned int uint32_t;
#  endif // _MSC_VER
#endif

#include <algorithm>
#include <cfloat>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

#ifdef _MSC_VER
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
class MathUtil 
{
private:
  /// This class.
  typedef MathUtil Self;

public:
  /// Portable test for "not a number" values.
  template<class T>
  static bool IsNaN( const T value )
  {
    return isnan( value );
  }

  /// Get double-precision not-a-number (NaN) value.
  static double GetDoubleNaN()
  {
    return Self::GetInternalNaN().m_Union.d;
  }
  
  /// Get single-precision not-a-number (NaN) value.
  static float GetFloatNaN()
  {
#if WORDS_BIGENDIAN
    return Self::GetInternalNaN().m_Union.f[0];
#else
    return Self::GetInternalNaN().m_Union.f[1];
#endif
  }
  
  /// Get double-precision infinite (Inf) value.
  static double GetDoubleInf()
  {
    return Self::GetInternalInf().m_Union.d;
  }

  /// Get single-precision infinite (Inf) value.
  static float GetFloatInf()
  {
#if WORDS_BIGENDIAN
    return Self::GetInternalInf().m_Union.f[0];
#else
    return Self::GetInternalInf().m_Union.f[1];
#endif
  }

  /// Unit-safe sin() function.
  static double Sin( const Units::Radians& radians )
  {
    return sin( radians.Value() );
  }
  
  /// Unit-safe cos() function.
  static double Cos( const Units::Radians& radians )
  {
    return cos( radians.Value() );
  }
  
  /// Unit-safe tan() function.
  static double Tan( const Units::Radians& radians )
  {
    return tan( radians.Value() );
  }
  
  /// Unit-safe asin() function.
  static const Units::Radians ArcSin( const double value )
  {
    return Units::Radians( asin( value ) );
  }

  /// Unit-safe acos() function.
  static const Units::Radians ArcCos( const double value )
  {
    return Units::Radians( acos( value ) );
  }
  
  /// Unit-safe aten() function.
  static const Units::Radians ArcTan( const double value )
  {
    return Units::Radians( atan( value ) );
  }
  
  /// Unit-safe atan2() function.
  static const Units::Radians ArcTan2( const double y, const double x )
  {
    return Units::Radians( atan2( y, x ) );
  }
  
  /// Return square of a float value.
  template<class T> static T Square ( const T a) { return a*a; }

  /** Return minimum of an array of ordered values.
   */
  template<class T> static T Min ( const int nValues, const T* Values )
  {
    T Result = Values[0];
    for ( int idx=1; idx<nValues; ++idx )
      Result = std::min( Result, Values[idx] );
    
    return Result;
  }
  
  /** Return minimum of an array of ordered values.
   */
  template<class T> static T Max ( const int nValues, const T* Values )
  {
    T Result = Values[0];
    for ( int idx=1; idx<nValues; ++idx )
      Result = std::max( Result, Values[idx] );
    
    return Result;
  }
  
  /// Return length of intersection of two intervals.
  template<class T> static T Intersect ( const T aMin, const T aMax, const T bMin, const T bMax )
  {
    return ( std::min( aMax, bMax ) - std::max( aMin, bMin ) );
  }

  /// Round float value to the nearest integer.
  template<class T> static int Round ( const T a ) { return (int)floor(a+0.5); }
  
  /** Return sign of float value.
   *@return -1 if a<0, 1 if a>0, 0 if a==0.
   */
  template<class T> static int Sign ( const T a ) { return (a<0)?-1:((a==0)?0:1); }

  /** Check if some float value is within a range.
   *@return 0 if value is in range, -1 if value is below minumum, 1 if value
   * is above maximum.
   */ 
  template<class T> static int CheckRange ( const T value, const T a, const T b ) 
  {
    const int sigA = Self::Sign(value-a);    
    return (sigA == Self::Sign(value-b))?sigA:0;
  }
  
  /// Compute p*log(p) for a single value.
  template<class T> static double plogp( const T p ) { return (p>0) ? p * log( p ) : 0.0; }
  
  /** Computes average of an array of float values.
   */
  template<class T> static
  T Mean
  ( const unsigned int nValues, const T* values );
  
  /** Computes average of a vector of float values.
   */
  template<class T> static
  T Mean
  ( const std::vector<T>& values );
  
  /** Computes variance of an array of float values.
    *\param unbiased If this flag is set (default: unset), then the variance
    * will be computed over nValues-1; otherwise over nValues.
   */
  template<class T> static
  T Variance
  ( const unsigned int nValues, const T* values, T mean, const bool unbiased = false );
  
  /** Computes variance of a vector of float values.
    *\param values Vector of values to compute variance from.
    *\param mean Previously computed mean of vector values.
    *\param unbiased If this flag is set (default: unset), then the variance
    * will be computed over nValues-1; otherwise over nValues.
   */
  template<class T> static
  T Variance
  ( const std::vector<T>& values, T mean, const bool unbiased = false );
  
  /** Normalized correlation coefficient between two float vectors.
   */
  template<class T> static
  T Correlation( const std::vector<T>& x, const std::vector<T>& y );
  
  /// Compute t-statistic from coefficient of correlation.
  static double TStatFromCorrelation
  ( const double r, //!< Coefficient of correlation as computed by MathUtil::Correlation function.
    const size_t df //!< Number of degrees of freedom
    );
  
  /// Compute probability from T-statistic.
  static double ProbabilityFromTStat
  ( const double t, //!< T-statistic as returned for example from MathUtil::TStatFromCorrelation function.
    const size_t df //!< Number of degrees of freedom.
    );
  
  /** Performs two-tailed unpaired t-test on two distributions.
   */
  template<class T> static
  T TTest ( const std::vector<T>& valuesX, const std::vector<T>& valuesY, T& t );
  
  /** Performs two-tailed unpaired t-test on two distributions.
   * Also return average value for each distribution.
   */
  template<class T> static
  T TTest ( const std::vector<T>& valuesX, const std::vector<T>& valuesY, T& t, T& avgX, T& avgY );
  
  /** Performs two-tailed paired t-test on two distributions.
   * Also return average value for each distribution.
   */
  template<class T> static
  T PairedTTest ( const std::vector<T>& valuesX, const std::vector<T>& valuesY, T& t, T& avgX, T& avgY );
  
  /** Performs one-sample t-test on distribution to test for zero mean.
   * Also return average value for each distribution.
   */
  template<class T> static
  T TTest ( const std::vector<T>& valuesX, T& t, T& avgX );
  
  /// Beta-i function.
  static double Betai( const double a, const double b, const double x );
  
  /// Beta-Cf function.
  static double BetaCf( const double a, const double b, const double x );

  /// GammaLn function.
  static double GammaLn( const double xx );

  /// Singular Value Decomposition
  static void SVD( Matrix2D<double> *U, size_t m, size_t n, std::vector<double> *W, Matrix2D<double> *V );

  /// Linear Regression using SVD results
  static void
  SVDLinearRegression( Matrix2D<double> *U, size_t m, size_t n, std::vector<double> *W, Matrix2D<double> *V, double *b, std::vector<double>& lm_params );
  
  /// Function that compares two floats; to be used in qsort().
  static inline int CompareFloat( const void *a, const void *b ) 
  {
    const float* A = static_cast<const float*>( a );
    const float* B = static_cast<const float*>( b );
    if ( *A > *B ) return +1;
    if ( *A < *B ) return -11;
    return 0;
  }
  
  /// Function that compares two doubles; to be used in qsort().
  static inline int CompareDouble( const void *a, const void *b ) 
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
  static inline double NormalRandom( const double sigma ) 
  {
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
   *@param sigma Standard deviation of the resulting distribution.
   *@param seed Random seed given to srandom() function.
   */
  static inline double NormalRandom( const double sigma, const unsigned int seed ) 
  {
    srandom( seed );
    return NormalRandom( sigma );
  }
  
  /// Uniform random number generator.
  static double UniformRandom();
  
  /// Compute eigensystem and eigenvalues for square real matrix using Jacobi rotation.
  template<class T> static void ComputeEigensystem( const Matrix2D<T>& matrix, Matrix2D<T>& eigensystem, std::vector<T>& eigenvalues );

  /// Compute eigenvalues for square real matrix using Jacobi rotation.
  template<class T> static void ComputeEigenvalues( const Matrix2D<T>& matrix, std::vector<T>& eigenvalues );

  /// Determinant of an n x n square matrix.
  template<class T> static T CholeskyDeterminant( const Matrix2D<T>& matrix, int n);

private:
  /// Helper class to initialize a constant 64bit field.
  class BitInitializer
  {
  public:
    /// Union that is used for initializing variables of different types with bitfield values.
    typedef union
    {
      /// Write bit pattern for special values here.
      uint32_t i[2];
      /// Read out for single-precision special value.
      float f[2];
      /// Read out for double-precision special value.
      double d;
    } InitializeUnion;
    
    /// The actual initializer memory.
    InitializeUnion m_Union;
    
  public:
    /// Initialize 64bit field with two given 32bit integers.
    BitInitializer( const uint32_t i0, const uint32_t i1 )
    {
      this->m_Union.i[0] = i0;
      this->m_Union.i[1] = i1;
    }
  };
  
  /// Get reference to internal representation of NaN.
  static const Self::BitInitializer& GetInternalNaN()
  {
#if WORDS_BIGENDIAN
    static const Self::BitInitializer bits( 0x7fffffff, 0xffffffff );
#else
    static const Self::BitInitializer bits( 0xffffffff, 0x7fffffff );
#endif
    return bits;
  }
  
  /// Get reference to internal representation of Inf.
  static const Self::BitInitializer& GetInternalInf()
  {
#if WORDS_BIGENDIAN
    static const Self::BitInitializer bits( 0x7f800000, 0x00000000 );
#else
    static const Self::BitInitializer bits( 0x00000000, 0x7f800000 );
#endif
    return bits;
  }
};

//@}

} // namespace cmtk

#include "cmtkMathUtil_Statistics.txx"

#endif // #ifndef __cmtkMathUtil_h_included_
