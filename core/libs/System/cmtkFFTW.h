/*
//
//  Copyright 2012 SRI International
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

#ifndef __cmtkFFTW_h_included_
#define __cmtkFFTW_h_included_

#include <cmtkconfig.h>

#include <fftw3.h>
#include <math.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Interface class to FFTW library for fast Fourier transforms.
 * Among other things, this class handles the initialization and destruction of global FFTW properties, such as threads support.
 */
class FFTW
{
public:
  /// This class.
  typedef FFTW Self;

  /// Get static instance.
  static Self& GetStatic()
  {
    static Self staticInstance;
    return staticInstance;
  }

  /// In-place complex multiplication.
  static void MultiplyInPlace( fftw_complex& lhs, const fftw_complex& rhs )
  {
    fftw_complex product;
    product[0] = lhs[0]*rhs[0] - lhs[1]*rhs[1];
    product[1] = lhs[1]*rhs[0] + lhs[0]*rhs[1];
    
    lhs[0] = product[0];
    lhs[1] = product[1];
  }

  /// Return sum of squares of real and imaginary parts of a complex number.
  static double SumOfSquares( const fftw_complex& c )
  {
    return c[0]*c[0] + c[1]*c[1];
  }

  /// Return magnitude of complex number.
  static double Magnitude( const fftw_complex& c )
  {
    return sqrt( Self::SumOfSquares( c ) );
  }

  /** Set number of threads for FFTW plans.
   * The number of threads applies only to subsequently created plans. Previously created plans
   * retain their original number of threads when executed.
   */
  void SetNumberOfThreads( const int nThreads )
  {
#ifdef CMTK_USE_SMP
    fftw_plan_with_nthreads( nThreads );
#endif
  }

protected:
  /// Constructor.
  FFTW()
  {
#ifdef CMTK_USE_SMP
    fftw_init_threads();
#endif
  }

  /// Destructor.
  ~FFTW()
  {
#ifdef CMTK_USE_SMP
    fftw_cleanup_threads();
#endif
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkFFTW_h_included_
