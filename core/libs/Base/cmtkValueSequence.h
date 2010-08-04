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

#ifndef __cmtkValueSequence_h_included_
#define __cmtkValueSequence_h_included_

#include <cmtkconfig.h>

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for computing characteristic values of number sequences.
 * Instances of this class take a sequence of (real) numbers via calls to the
 * Proceed() member function. After last last number in the sequence, the
 * object can be queried for minimum and maximum values (both absolute and
 * standard), variance, and average value. Further calls to Proceed() are
 * allowed, thus enabling incremental computations.
 *\note This class performs a lot of potentially unnecessary computations if,
 * for example, only the sum of values is to be computed. If you need efficient
 * evaluation of more elementary expressions, you are encouraged to implement
 * these directly; this present class is merely a convenience tool.
 */
template<class T=float> 
class ValueSequence 
{
public:
  /// This class.
  typedef ValueSequence<T> Self;

  /// Default constructor.
  ValueSequence() 
  { 
    this->Reset(); 
  }

  /// Reset all computations.
  void Reset() 
  {
    NValues = 0; Sum = 0; SumAbs = 0; SumOfSquares = 0; 
    Minimum = Maximum = MinimumAbs = MaximumAbs = 0; 
  }
  
  /// Proceed with the next number in the sequence.
  void Proceed( const T v ) 
  {
    if ( ! NValues ) 
      {
      Minimum = Maximum = v;
      MinimumAbs = MaximumAbs = fabs( v );
      } 
    else
      {
      if ( v < Minimum ) Minimum = v;
      if ( v > Maximum ) Maximum = v;
      if ( fabs( v ) < MinimumAbs ) MinimumAbs = fabs( v );
      if ( fabs( v ) > MaximumAbs ) MaximumAbs = fabs( v );
      }
    ++NValues; Sum += v; SumAbs += fabs( v ); SumOfSquares += v*v; 
  }

  /// Return minimum of all values.
  double GetMinimum() const { return Minimum; }

  /// Return maximum of all values.
  double GetMaximum() const { return Maximum; }

  /// Return minimum of all absolute values.
  double GetMinimumAbs() const { return MinimumAbs; }

  /// Return maximum of all absolute values.
  double GetMaximumAbs() const { return MaximumAbs; }

  /// Return total number of values.
  int GetNValues() const { return NValues; }

  /// Return variance of all values.
  double GetVariance( const bool unbiased = true ) const
  { 
    const double mu = this->GetAverage();
    return ( NValues * mu * mu - 2 * mu * Sum + SumOfSquares ) / ( unbiased ? (NValues-1) : NValues );
  }

  /// Return sum of all values.
  double GetSum() const { return static_cast<double>( Sum ); }

  /// Return sum of squres of all values.
  double GetSumOfSquares() const { return static_cast<double>( SumOfSquares ); }

  /// Return average value.
  double GetAverage() const { return static_cast<double>( Sum / NValues ); }

  /// Return average value.
  double GetAverageAbs() const { return static_cast<double>( SumAbs / NValues ); }

  /// Assignment operator.
  ValueSequence<T>& operator=( const ValueSequence<T>& other );

private:
  /// Sum of all values.
  T Sum;

  /// Sum of all absolute values.
  T SumAbs;

  /// Sum of the squares of all values.
  T SumOfSquares;

  /// Number of values.
  int NValues;

  /// Minimum value so far.
  T Minimum;

  /// Maximum value so far.
  T Maximum;

  /// Minimum absolute value so far.
  T MinimumAbs;

  /// Maximum absolute value so far.
  T MaximumAbs;

  /// Allow addition operator direct access.
  template<class TT>
  friend ValueSequence<TT> operator+( const ValueSequence<TT>& a, const ValueSequence<TT>& b );
};

//@}

} // namespace cmtk

#include "cmtkValueSequence.txx"

#endif // #ifndef __cmtkValueSequence_h_included_
