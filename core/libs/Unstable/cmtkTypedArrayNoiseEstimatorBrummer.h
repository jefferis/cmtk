/*
//
//  Copyright 2008-2010 SRI International
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

#include <cmtkconfig.h>

#include "Base/cmtkTypedArray.h"
#include "Base/cmtkTypedArrayNoiseEstimatorNaiveGaussian.h"
#include "Base/cmtkHistogram.h"

namespace 
cmtk
{

/** Estimate noise level in data stored in a TypedArray.
 * Estimate Rician noise variance using Brummer's method.
 *\author Mike Hasak
 */
class
TypedArrayNoiseEstimatorBrummer :
  /// Inherit interface from naive estimator class.
  protected TypedArrayNoiseEstimatorNaiveGaussian
{
public:
  /// This class.
  typedef TypedArrayNoiseEstimatorBrummer Self;

  /// Base class.
  typedef TypedArrayNoiseEstimatorNaiveGaussian Superclass;

  /// Constructor.
  TypedArrayNoiseEstimatorBrummer( const TypedArray* data, const size_t histogramBins = 255 );
  
protected:
  /// Default constructor; should not be invoked by user code.
  TypedArrayNoiseEstimatorBrummer() {}

  /** Compute bias in an ML noise estimate.
   *  Eq. 26 from Sijbers et al, 2007 
   */
  static double SijbersBiasHat ( const Histogram<unsigned int>::SmartPtr histogram, const double sigmaHat, const int numBinsToUse );

  /** Compute the log-likelihood for the ML noise estimate.
   *  Eq. 26 from Sijbers et al, 2007 
   */
  static double SijbersLogLikelihood ( const Histogram<unsigned int>::SmartPtr histogram, const double sigma, const int numBinsToUse );

  /** Compute variance of an ML noise estimate.
   *  Eq. 21 from Sijbers et al, 2007 
   */
  static double SijbersVarHat ( const Histogram<unsigned int>::SmartPtr histogram, const double sigmaHat, const int numBinsToUse );

  /** Figure out how many bins to use in ML noise estimate.
   *  Section 2.3.1 from Sijbers et al, 2007
   */
  static double EstimateNumBinsToUse ( const TypedArray* data, const Histogram<unsigned int>::SmartPtr histogram, const double sigmaHat );
};

} // namespace cmtk
