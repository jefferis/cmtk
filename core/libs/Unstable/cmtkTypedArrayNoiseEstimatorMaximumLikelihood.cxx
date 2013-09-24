/*
//
//  Copyright 2008-2011, 2013 SRI International
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

#include "cmtkTypedArrayNoiseEstimatorMaximumLikelihood.h"

namespace
cmtk
{

/*   Maximum-likelihood estimation of the noise variance.
 */
TypedArrayNoiseEstimatorMaximumLikelihood::TypedArrayNoiseEstimatorMaximumLikelihood
( const TypedArray* data, const size_t histogramBins ) 
{
//  Histogram<unsigned int>::SmartPtr histogram( data->GetHistogram( histogramBins ) );
  UNUSED(histogramBins);

  size_t overriddenHistogramBins = 255;
  Histogram<unsigned int>::SmartPtr histogram( data->GetHistogram( overriddenHistogramBins ) );

  // Find first maximum
  size_t idx = 0;
  while ( (idx < (overriddenHistogramBins-1)) && ( (*histogram)[idx] <= (*histogram)[idx+1] ) )
    {
    ++idx;
    }
  // Now find following minimum
  while ( (idx < (overriddenHistogramBins-1)) && ( (*histogram)[idx] > (*histogram)[idx+1] ) )
    {
    ++idx;
    }
  
   /*  Part of estimating the noise variance by this method is
    *  deciding how many histogram bins to use in the computations
    *  below (i.e. determining the value of K in Sijbers).  
    *  K is called numBinsToUse here, and we'll initialize it
    *  to the bin with the minimum following the first maximum.
    *  It will then be refined using the EstimateNumBinsToUse method.
    */
  int numBinsToUse = idx;

  float sigmaUpperBound = 100;
  
  // Step size in search for sigma
  double stepSize = 0.1;

  // Threshold for how close two iterations of the below should
  // be before deciding we've found a good sigma
  double sigmaThresh = 0.1; 

  // maxLikelySigma corresponds to sigma-hat in the Sijbers paper
  double maxLikelySigma = 23.0;
  double prevMaxLikelySigma = 0.0;

  while ( sqrt( pow( maxLikelySigma - prevMaxLikelySigma, 2 ) ) > sigmaThresh )
    {
    
    std::cout << "K: " << numBinsToUse << "\tprevMaxLikelySigma: " << prevMaxLikelySigma << "\tmaxLikelySigma: " << maxLikelySigma << std::endl;
    prevMaxLikelySigma = maxLikelySigma;

    double minSigmaMinimizer = std::numeric_limits<double>::max(); 

    const double bin0Squared = pow( histogram->BinToValue( 0 ), 2 );
    const double binKSquared = pow( histogram->BinToValue(numBinsToUse-1), 2 );

     /*  This for loop iterates through values of sigma, assigning
      *  the one that produces the least value of minSigmaMinimizer
      *  to maxLikelySigma.
      */
    for ( double sigma = 1; sigma < sigmaUpperBound; sigma += stepSize )
      {
      const double twoSigmaSquared = 2 * pow( sigma, 2 );
      double sumForSearch = 0.0; 
  
      for ( int j = 1; j <= numBinsToUse; j++ )
        {
        const double curAddend = (*histogram)[j] 
	  * ( exp( -1 * ( pow( (double)(histogram->BinToValue(j-1)), 2 ) / twoSigmaSquared ) )
	      - exp( -1 * ( pow( (double)(histogram->BinToValue( j )), 2 ) / twoSigmaSquared ) ) );  
        sumForSearch += curAddend;
        }
      
      const double curSigmaMinimizer = ( exp( -bin0Squared / twoSigmaSquared ) - exp( -binKSquared / twoSigmaSquared ) ) - sumForSearch;
      
      if ( curSigmaMinimizer < minSigmaMinimizer ) 
        {
        minSigmaMinimizer = curSigmaMinimizer;
        maxLikelySigma = sigma;
        }
      }
    numBinsToUse = static_cast<int>( Superclass::EstimateNumBinsToUse( data, histogram, maxLikelySigma ) );
    }
    
    std::cout << maxLikelySigma << std::endl;

  //return static_cast<Types::DataItem>( maxLikelySigma );
    this->m_NoiseLevelSigma = static_cast<Types::DataItem>( 0 );
}

} // namespace cmtk
