/*
//
//  Copyright 2008-2011 SRI International
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

#include "cmtkTypedArrayNoiseEstimatorBrummer.h"

namespace
cmtk
{

/*   Brummer estimation of the noise variance.
 *   (As described in the Sijbers paper.)
 */
TypedArrayNoiseEstimatorBrummer::TypedArrayNoiseEstimatorBrummer
( const TypedArray* data, const size_t histogramBins )
{
  UNUSED(histogramBins);

  const size_t histogramBinsActual = 255;
  Histogram<unsigned int>::SmartPtr histogram( data->GetHistogram( histogramBinsActual ) );

//  int numBinsToUse;

  double sampleMean = 0;
  for ( size_t i = 0; i < histogramBinsActual; i++ )
    {
      sampleMean += (*histogram)[i] * histogram->BinToValue(i);
    }
  sampleMean /= histogram->SampleCount();
  std::cout << " sampleMean: " << sampleMean << std::endl;
  std::cout << " sampleCount: " << histogram->SampleCount() << std::endl;
  std::cout << " 1 / sampleCount: " << (1.0 / histogram->SampleCount()) << std::endl;
  
  double sumSquaredSampleDiffs = 0;
  for ( size_t i = 0; i < histogramBinsActual; i++ )
    {
      sumSquaredSampleDiffs += (*histogram)[i] * pow( ( histogram->BinToValue(i) - sampleMean ), 2 );
    }
  std::cout << " sumSquaredSampleDiffs: " << sumSquaredSampleDiffs << std::endl;
  
  const double sampleStdDev = sqrt( ( 1.0 / histogram->SampleCount() ) * sumSquaredSampleDiffs );
  std::cout << " sampleStdDev: " << sampleStdDev << std::endl;
  std::cout << " npow: " << pow( (double)histogram->SampleCount(), 1/5 ) << std::endl;

  /* Find an initial estimate of sigma using the Chang method.
   * (Chang L-C et al 2005 "An automatic method for estimating noise-induced signal 
   *  variance in magnitude-reconstructed magnetic resonance images")
   */
  const double changH = 1.0 / ( 1.06 * sampleStdDev * pow( (double)histogram->SampleCount(), 1/5 ) );
  double changSigma = 0;
  double sigmaUpperBound = 100;
  double valueToMaximize = std::numeric_limits<double>::min(); 
  for ( double maybeChangSigma = 1.0; maybeChangSigma < sigmaUpperBound; maybeChangSigma += 1.5 )
    {

    double changSum = 0;
    for ( size_t i = 0; i < histogramBinsActual; i++ )
      {
      const double valToAdd = ( (*histogram)[i]  
				* ( 1.0 / sqrt(2.0 * M_PI) )
				* exp( ( -1.0 * log( pow( (double)( maybeChangSigma - histogram->BinToValue(i) ) / changH, 2 )
                                                       / 2.0) ) ) );
        changSum += valToAdd;
        // std::cout << " valToAdd: " << valToAdd << std::endl;
      }

    std::cout << " changSum: " << changSum << std::endl;

    double curValueToMaximize = ( 1.0 / ( histogram->SampleCount() * changH ) ) * changSum;
    
    std::cout << " curValueToMaximize: " << curValueToMaximize << std::endl;

    if ( curValueToMaximize > valueToMaximize )
      {
      valueToMaximize = curValueToMaximize;
      changSigma = maybeChangSigma;
      std::cout << "New max: " << curValueToMaximize << " New changSigma: " << changSigma << std::endl;
      }

    }

  /*  Using the Chang sigma as a starting estimate, use
   *  the Brummer method to get a final sigma estimate.
   */
  valueToMaximize = std::numeric_limits<double>::min();
  double sigma0 = changSigma;
  //double sigma0 = 27.0;
  double maxAmplitude = 2000;
  double brummerSigma = changSigma;
  for ( double maybeSigma = 0; maybeSigma < histogramBinsActual; maybeSigma += 1.0 )
    {
    for ( double maybeAmplitude = 0; maybeAmplitude < maxAmplitude; maybeAmplitude += 10.0 )
      {
      double brummerSum = 0;
      for ( int i = 0; i <= 2 * sigma0; i++ )
        {
	brummerSum += log( pow( (double)(*histogram)[i] 
				- maybeAmplitude 
				* ( histogram->BinToValue(i) / pow( maybeSigma, 2) )
				* exp( -1 * ( pow( (double)histogram->BinToValue(i),2) / (2 * pow( maybeSigma, 2 ) ) ) ) 
                                , 2 ) );
        }
      std::cout << " sigma: " << maybeSigma << " K: " << maybeAmplitude << " brummerSum: " << brummerSum << std::endl;
      if ( brummerSum > valueToMaximize )
        {
          valueToMaximize = brummerSum;
          brummerSigma = maybeSigma;
          std::cout << "new valueToMaximize: " << valueToMaximize << std::endl;
          std::cout << "new brummerSigma: " << brummerSigma << std::endl;
        }
      }
    }
  std::cerr << brummerSigma << std::endl;
  //return static_cast<Types::DataItem>( maxLikelySigma );
  
  this->m_NoiseLevelSigma = static_cast<Types::DataItem>( 0 );
}

double
TypedArrayNoiseEstimatorBrummer::SijbersBiasHat
( const Histogram<unsigned int>::SmartPtr histogram, const double sigmaHat, const int numBinsToUse )
{

  double K = numBinsToUse;

  /* 
   *  Count the number of samples in the partial histogram
   *  (bins 0 thru K-1 of histogram)
   */
  int partialHistCount = 0;
  for ( int i = 0; i < K; i++ )
    {
    partialHistCount += (*histogram)[i];
    }

  /* 
   *  Pre-calculate some recurring quantities
   */
  const double twoSigmaSquared = 2 * sigmaHat * sigmaHat;
  const double denom = 
      exp( -1 * ( log( (double)pow( histogram->BinToValue( 0 ), 2 ) / twoSigmaSquared ) ) )
    - exp( -1 * ( log( (double)pow( histogram->BinToValue( K ), 2 ) / twoSigmaSquared ) ) ) ; 
  // This will be a denominator below, needs to be != 0             
  assert ( denom != 0 );

  /* 
   * lambdaStar is composed of two sums: one using
   * bins 0 thru K-1 of histogram, and
   * one using bins K through the end of the
   * histogram.  This loop computes the first sum.
   */
  double lambdaStar = 0.0;
  for ( int i = 1; i < K; i++ )
    {
    const double quotient = 
      ( exp( -1 * ( log( (double)pow( histogram->BinToValue(i-1), 2 ) ) ) / twoSigmaSquared ) 
	- exp( -1 * ( log( pow( (double)histogram->BinToValue( i ), 2 ) ) ) / twoSigmaSquared ) )
      / denom;
    const double f = partialHistCount * quotient;
    assert ( f != 0 );
    lambdaStar += pow( f - (*histogram)[i], 2 ) / f;
    }  

  /* 
   *  Second sum for lambdaStar
   */
  for ( size_t i = static_cast<size_t>( K ); i < histogram->GetNumBins(); i++ )
    {
    const double quotient = 
      ( exp( -1 * ( log( pow( (double)histogram->BinToValue(i-1), 2 ) ) ) / twoSigmaSquared ) 
	- exp( -1 * ( log( pow( (double)histogram->BinToValue( i ), 2 ) ) ) / twoSigmaSquared ) )
      / denom;
    const double f = partialHistCount * quotient;
    double numerator = ( f - (*histogram)[i] );
    if ( numerator < 0 ) numerator = 0 ;
    assert ( f != 0 );
    lambdaStar += pow( numerator, 2 ) / f;
                      
    }

  /* 
   *  Compute M (variable name same as in Sijbers)
   */
  int M = 0;
  for ( size_t i = static_cast<size_t>( K ); i < histogram->GetNumBins(); i++ )
    {
    const double quotient = 
      ( exp( -1 * ( log( pow( (double)histogram->BinToValue(i-1), 2 ) ) / twoSigmaSquared ) ) 
	- exp( -1 * ( log( pow( (double)histogram->BinToValue( i ), 2 ) ) / twoSigmaSquared ) ) )
      / denom;
    const double f = partialHistCount * quotient;
    const double applyStepFilterTo = ( f - (*histogram)[i] );
    M += ( applyStepFilterTo >= 0 ) ? 1 : 0;
    }
//  std::cout << "M: " << M << std::endl;

  /* 
   *  biasHat is a simple formula using the results 
   *  from above.
   */
  const double biasHat = ( lambdaStar - ( K - 2 + M ) )
                   / sqrt( K - 2 + M );
  std::cout << "biasHat: " << biasHat << std::endl;

  //return biasHat;
  return 0;
}

double
TypedArrayNoiseEstimatorBrummer::SijbersLogLikelihood
( const Histogram<unsigned int>::SmartPtr histogram, const double sigma, const int numBinsToUse )
{

  double s2 = 2 * sigma * sigma;
  int K = numBinsToUse;

  double logLikelihood = 0.0;
  for ( int i = 0; i < K; i++ )
    {
    const double curCount = (*histogram)[i];
    const double quotient =
      ( exp( -1 * ( log( pow( (double)histogram->BinToValue(i-1), 2 ) ) / s2 ) )
	- exp( -1 * ( log( pow( (double)histogram->BinToValue( i ), 2 ) ) / s2 ) ) )
      /
      ( exp( -1 * ( log( pow( (double)histogram->BinToValue( 0 ), 2 ) ) / s2 ) )
	- exp( -1 * ( log( pow( (double)histogram->BinToValue( K ), 2 ) ) / s2 ) ) );
    
    std::cout << exp( -1 * ( log( pow( histogram->BinToValue(i-1), 2 ) ) / s2 ) ) << std::endl;
    std::cout << exp( -1 * ( log( pow( histogram->BinToValue( i ), 2 ) ) / s2 ) ) << std::endl;
    std::cout << exp( -1 * ( log( pow( histogram->BinToValue( 0 ), 2 ) ) / s2 ) ) << std::endl;
    std::cout << exp( -1 * ( log( pow( histogram->BinToValue( K ), 2 ) ) / s2 ) ) << std::endl;
    
    logLikelihood += curCount * quotient;
    }
  return logLikelihood;
}

double
TypedArrayNoiseEstimatorBrummer::SijbersVarHat
( const Histogram<unsigned int>::SmartPtr histogram, const double sigmaHat, const int numBinsToUse )
{

  double h = 0.01;
  
  const double fXminusH = SijbersLogLikelihood( histogram, sigmaHat - h, numBinsToUse ); 
  const double fX = SijbersLogLikelihood( histogram, sigmaHat, numBinsToUse ); 
  const double fXplusH = SijbersLogLikelihood( histogram, sigmaHat + h, numBinsToUse ); 

  std::cout << "fXminusH: " << fXminusH << std::endl;
  std::cout << "fX: " << fX << std::endl;
  std::cout << "fXplusH: " << fXplusH << std::endl;
  
  const double secondDerivOfLogLikelihood = ( fXminusH - 2 * fX + fXplusH ) / ( h * h );

  const double varHat = -1.0 * 1 / secondDerivOfLogLikelihood;
  std::cout << "varHat: " << varHat << std::endl;
  //return varHat;
  return 0;
}

double
TypedArrayNoiseEstimatorBrummer::EstimateNumBinsToUse
( const TypedArray* data, const Histogram<unsigned int>::SmartPtr histogram, const double sigmaHat )
{
  UNUSED(data);

  double bestK = 0;
  double curKMinimizer;
  double minKMinimizer = std::numeric_limits<double>::max(); 
 
  double maxK = histogram->GetNumBins() / 2;
  
  for ( int curK = 1; curK < maxK; curK++ )
    {
    curKMinimizer = SijbersBiasHat( histogram, sigmaHat, curK )+ SijbersVarHat( histogram, sigmaHat, curK );
//std::cout << "minKMinimizer: " << minKMinimizer << std::endl;
//std::cout << "curKMinimizer: " << curKMinimizer << std::endl;
    if ( curKMinimizer < minKMinimizer ) 
      {
      minKMinimizer = curKMinimizer;
      bestK = curK;
      }
    }

//std::cout << std::endl << "bestK: " << bestK << std::endl;

  return bestK;
}

} // namespace cmtk
