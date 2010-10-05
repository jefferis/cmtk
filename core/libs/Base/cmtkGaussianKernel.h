/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkGaussianKernel_h_included_
#define __cmtkGaussianKernel_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUnits.h>
#include <Base/cmtkMathUtil.h>

#include <vector>

/** \addtogroup Base */
//@{

namespace
cmtk
{

/// Utility class for generating Gaussian kernels.
template<class TFloat=double>
class GaussianKernel
{
public:
  /// This class.
  typedef GaussianKernel<TFloat> Self;

  /// Create symmetric kernel.
  static std::vector<TFloat> GetSymmetricKernel( const Units::GaussianSigma& sigma /*!< Sigma parameter (standard deviation) of the kernel */, 
						 const TFloat maxError = 1e-5 /*!< Maximum approximation error: the kernel radius is computed so that truncated elements are below this value */ )
  {
    const double normFactor = 1.0/(sqrt(2*M_PI) * sigma.Value());
    const size_t radius = Self::GetRadius( sigma, normFactor, maxError );
    
    std::vector<TFloat> kernel( 2 * radius + 1 );
    for ( size_t i = 0; i <= radius; ++i )
      {
      kernel[radius-i] = kernel[radius+i] = normFactor * exp( -MathUtil::Square( 1.0 * i / sigma.Value() ) / 2 );
      }

    return kernel;
  }

  /// Create half kernel, starting with center element.
  static std::vector<TFloat> GetHalfKernel( const Units::GaussianSigma& sigma /*!< Sigma parameter (standard deviation) of the kernel */, 
					    const TFloat maxError = 1e-5 /*!< Maximum approximation error: the kernel radius is computed so that truncated elements are below this value */ )
  {
    const double normFactor = 1.0/(sqrt(2*M_PI) * sigma.Value());
    const size_t radius = Self::GetRadius( sigma, normFactor, maxError );
    
    std::vector<TFloat> kernel( radius + 1 );
    for ( size_t i = 0; i <= radius; ++i )
      {
      kernel[i] = normFactor * exp( -MathUtil::Square( 1.0 * i / sigma.Value() ) / 2 );
      }
    
    return kernel;
  }

private:
  /// Compute kernel radius based on maximum approximation error and normalization factor.
  static TFloat GetRadius(  const Units::GaussianSigma& sigma, const TFloat normFactor, const TFloat maxError )
  {
    if ( maxError >= normFactor ) // if normFactor is less than max error, then we really need no kernel at all
      return 0;
    else
      return static_cast<size_t>( sqrt( -2.0 * log( maxError / normFactor ) ) * sigma.Value() );
  }
};

} // namespace cmtk

//@}

#endif // #ifndef __cmtkGaussianKernel_h_included_
