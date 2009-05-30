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

#include <cmtkGroupwiseRegistrationRMIFunctional.h>

#include <cmtkMathUtil.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TXform>
GroupwiseRegistrationRMIFunctional<TXform>
::GroupwiseRegistrationRMIFunctional()
  : m_MutualInformation( false )
{
  this->SetNumberOfHistogramBins( 255 );
}

template<class TXform>
GroupwiseRegistrationRMIFunctional<TXform>
::~GroupwiseRegistrationRMIFunctional()
{
}

template<class TXform>
void
GroupwiseRegistrationRMIFunctional<TXform>
::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid,
  const int downsample, 
  const bool useTemplateData )
{  
  this->Superclass::SetTemplateGrid( templateGrid, downsample, useTemplateData );
}

template<class TXform>
typename GroupwiseRegistrationRMIFunctional<TXform>::ReturnType
GroupwiseRegistrationRMIFunctional<TXform>
::Evaluate()
{
  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  const size_t numberOfImages = this->m_ImageVector.size();

  this->m_CovarianceMatrix.Resize( numberOfImages, numberOfImages ); // needs no reset
  this->m_TotalNumberOfSamples = 0;

  this->m_SumOfProductsMatrix.resize( numberOfImages * (1+numberOfImages) / 2 );
  std::fill( this->m_SumOfProductsMatrix.begin(), this->m_SumOfProductsMatrix.end(), 0 );

  this->m_SumsVector.resize( numberOfImages );
  std::fill( this->m_SumsVector.begin(), this->m_SumsVector.end(), 0 );

  ThreadParameterArray<Self,EvaluateThreadParameters> params( this, numberOfThreads );
  this->m_ThreadSumOfProductsMatrix.resize( numberOfThreads );
  this->m_ThreadSumsVector.resize( numberOfThreads );
  
  if ( this->m_ProbabilisticSamples.size() )
    params.RunInParallel( &EvaluateProbabilisticThread );
  else
    params.RunInParallel( &EvaluateThread );

#ifdef CMTK_BUILD_MPI
  SumsAndProductsVectorType tmpVector( this->m_SumOfProductsMatrix.size() );
  MPI::COMM_WORLD.Allreduce( &(this->m_SumOfProductsMatrix[0]), &(tmpVector[0]), this->m_SumOfProductsMatrix.size(), MPI::LONG, MPI::SUM );
  std::copy( tmpVector.begin(), tmpVector.end(), this->m_SumOfProductsMatrix.begin() );
  
  tmpVector.resize( this->m_SumsVector.size() );
  MPI::COMM_WORLD.Allreduce( &(this->m_SumsVector[0]), &(tmpVector[0]), this->m_SumsVector.size(), MPI::LONG, MPI::SUM );
  std::copy( tmpVector.begin(), tmpVector.end(), this->m_SumsVector.begin() );

  unsigned int totalNumberOfSamples = this->m_TotalNumberOfSamples;
  MPI::COMM_WORLD.Allreduce( &totalNumberOfSamples, &this->m_TotalNumberOfSamples, 1, MPI::UNSIGNED, MPI::SUM );
#endif

  const float result = this->GetMetric( this->m_SumOfProductsMatrix, this->m_SumsVector, this->m_TotalNumberOfSamples, this->m_CovarianceMatrix );

  return result;
}

template<class TXform>
typename GroupwiseRegistrationRMIFunctional<TXform>::ReturnType
GroupwiseRegistrationRMIFunctional<TXform>
::GetMetric
( const SumsAndProductsVectorType& sumOfProductsMatrix, 
  const SumsAndProductsVectorType& sumsVector,
  const unsigned int totalNumberOfSamples,
  typename GroupwiseRegistrationRMIFunctional<TXform>::CovarianceMatrixType& covarianceMatrix ) const
{
  const size_t imagesFrom = this->m_ActiveImagesFrom;
  const size_t imagesTo = this->m_ActiveImagesTo;
  const size_t numberOfImages = imagesTo - imagesFrom;

  size_t midx = 0;
  for ( size_t j = 0; j < numberOfImages; ++j )
    {
    for ( size_t i = 0; i <= j; ++i, ++midx )
      {
      covarianceMatrix[i][j] = covarianceMatrix[j][i] = 
	(sumOfProductsMatrix[midx] - ((1.0 * sumsVector[i] * sumsVector[j]) / totalNumberOfSamples)) / totalNumberOfSamples;
      }
    }
  
  Array<typename Self::ReturnType> eigenvalues( numberOfImages );
  MathUtil::ComputeEigenvalues<typename Self::ReturnType>( covarianceMatrix, eigenvalues );

  const typename Self::ReturnType EIGENVALUE_THRESHOLD = 1e-6;
  typename Self::ReturnType determinant = 1.0;
  for ( size_t i = 0; i < numberOfImages; ++i )
    {
    if ( eigenvalues[i] > EIGENVALUE_THRESHOLD )
      determinant *= eigenvalues[i];
    }
  
  if ( determinant > 0 )
    {
    const static double alpha = 1.41893853320467;
    typename Self::ReturnType metric = numberOfImages*alpha + .5*log( determinant );
    if ( this->m_MutualInformation )
      {
      // subtract marginal entropy estimates (diagonal entries of the covariance matrix)
      for ( size_t i = 0; i < numberOfImages; ++i )
	if ( eigenvalues[i] > EIGENVALUE_THRESHOLD )
	  metric -= (alpha + .5*log( eigenvalues[i] ));
      }
    
    return -metric;
    }
  else
    {
    return -FLT_MAX;
    }
}

template<class TXform>
CMTK_THREAD_RETURN_TYPE
GroupwiseRegistrationRMIFunctional<TXform>
::EvaluateThread
( void *args )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  
  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const size_t numberOfImages = imagesTo - imagesFrom;
  
  const byte paddingValue = ThisConst->m_PaddingValue;

  SumsAndProductsVectorType& sumOfProductsMatrix = This->m_ThreadSumOfProductsMatrix[threadID];
  sumOfProductsMatrix.resize( numberOfImages * (1+numberOfImages) / 2 );
  std::fill( sumOfProductsMatrix.begin(), sumOfProductsMatrix.end(), 0 );
  
  SumsAndProductsVectorType& sumsVector = This->m_ThreadSumsVector[threadID];
  sumsVector.resize( numberOfImages );
  std::fill( sumsVector.begin(), sumsVector.end(), 0 );

  size_t totalNumberOfSamples = 0;

  const size_t numberOfPixels = ThisConst->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI  
  const size_t pixelsPerThread = numberOfPixels / ( numberOfThreads * ThisConst->m_SizeMPI );
  const size_t pixelFrom = ( threadID + ThisConst->m_RankMPI * numberOfThreads ) * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerThread );
#else
  const size_t pixelsPerThread = (numberOfPixels / numberOfThreads);
  const size_t pixelFrom = threadID * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerThread );
#endif

  for ( size_t ofs = pixelFrom; ofs < pixelTo; ++ofs )
    {
    bool allValid = This->m_Data[imagesFrom][ofs] != paddingValue;
    for ( size_t j = imagesFrom+1; allValid && (j < imagesTo); ++j )
      {
      allValid = This->m_Data[j][ofs] != paddingValue;
      }

    if ( allValid )
      {
      size_t midx = 0;
      for ( size_t j = imagesFrom; j < imagesTo; ++j )
	{
	const int dataJ = This->m_Data[j][ofs];
	sumsVector[j-imagesFrom] += dataJ;
	  
	for ( size_t i = imagesFrom; i <= j; ++i, ++midx )
	  {
	  const int dataI = This->m_Data[i][ofs];
	  sumOfProductsMatrix[midx] += dataI * dataJ;
	  ++totalNumberOfSamples;
	  }
	}
      }
    }

  // add our contribution to total results
  This->m_MutexLock.Lock();
  size_t midx = 0;
  for ( size_t j = imagesFrom; j < imagesTo; ++j )
    {
    This->m_SumsVector[j-imagesFrom] += sumsVector[j-imagesFrom];
    
    for ( size_t i = imagesFrom; i <= j; ++i, ++midx )
      {
      This->m_SumOfProductsMatrix[midx] += sumOfProductsMatrix[midx];
      }
    }
  This->m_TotalNumberOfSamples += totalNumberOfSamples;
  This->m_MutexLock.Unlock();

  return CMTK_THREAD_RETURN_VALUE;
}

template<class TXform>
CMTK_THREAD_RETURN_TYPE
GroupwiseRegistrationRMIFunctional<TXform>
::EvaluateProbabilisticThread
( void *args )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  
  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const size_t numberOfImages = imagesTo - imagesFrom;
  
  const byte paddingValue = ThisConst->m_PaddingValue;

  SumsAndProductsVectorType& sumOfProductsMatrix = This->m_ThreadSumOfProductsMatrix[threadID];
  sumOfProductsMatrix.resize( numberOfImages * (1+numberOfImages) / 2 );
  std::fill( sumOfProductsMatrix.begin(), sumOfProductsMatrix.end(), 0 );
  
  SumsAndProductsVectorType& sumsVector = This->m_ThreadSumsVector[threadID];
  sumsVector.resize( numberOfImages );
  std::fill( sumsVector.begin(), sumsVector.end(), 0 );

  const size_t numberOfSamples = ThisConst->m_ProbabilisticSamples.size();
#ifdef CMTK_BUILD_MPI  
  const size_t samplesPerThread = numberOfSamples / ( numberOfThreads * ThisConst->m_SizeMPI );
  const size_t sampleFrom = ( threadID + ThisConst->m_RankMPI * numberOfThreads ) * samplesPerThread;
  const size_t sampleTo = std::min( numberOfSamples, sampleFrom + samplesPerThread );
#else
  const size_t samplesPerThread = numberOfSamples / numberOfThreads;
  const size_t sampleFrom = threadID * samplesPerThread;
  const size_t sampleTo = std::min( numberOfSamples, sampleFrom + samplesPerThread );
#endif

  size_t totalNumberOfSamples = 0;
  for ( size_t sample = sampleFrom; sample < sampleTo; ++sample )
    {
    bool allValid = This->m_Data[imagesFrom][sample] != paddingValue;
    for ( size_t j = imagesFrom+1; allValid && (j < imagesTo); ++j )
      {
      allValid = (This->m_Data[j][sample] != paddingValue);
      }
    
    if ( allValid )
      {
      size_t midx = 0;
      for ( size_t j = imagesFrom; j < imagesTo; ++j )
	{
	const int dataJ = This->m_Data[j][sample];
	sumsVector[j-imagesFrom] += dataJ;
	  
	for ( size_t i = imagesFrom; i <= j; ++i, ++midx )
	  {
	  const int dataI = This->m_Data[i][sample];
	  sumOfProductsMatrix[midx] += dataI * dataJ;
	  ++totalNumberOfSamples;
	  }
	}
      }
    }
  
  // add our contribution to total results
  This->m_MutexLock.Lock();
  size_t midx = 0;
  for ( size_t j = imagesFrom; j < imagesTo; ++j )
    {
    This->m_SumsVector[j-imagesFrom] += sumsVector[j-imagesFrom];
    for ( size_t i = imagesFrom; i <= j; ++i, ++midx )
      {
      This->m_SumOfProductsMatrix[midx] += sumOfProductsMatrix[midx];    
      }
    }
  This->m_TotalNumberOfSamples += totalNumberOfSamples;
  This->m_MutexLock.Unlock();
  
  return CMTK_THREAD_RETURN_VALUE;
}

template<class TXform>
typename GroupwiseRegistrationRMIFunctional<TXform>::ReturnType
GroupwiseRegistrationRMIFunctional<TXform>
::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
  const typename Self::ReturnType baseValue = this->EvaluateAt( v );

  for ( size_t param = 0; param < this->ParamVectorDim(); ++param )
    {
    g[param] = 0.0;

    const size_t imageIndex = param / this->m_ParametersPerXform;
    const size_t paramIndex = param % this->m_ParametersPerXform;

    const Types::Coordinate pStep = this->GetParamStep( param, step );
    if ( pStep > 0 )
      {
      byte* tmp = this->m_Data[imageIndex];
      this->m_Data[imageIndex] = &(this->m_TempData[0]);

      const Types::Coordinate p0 = v[param];

      this->SetParameter( imageIndex, paramIndex, p0 + pStep );
      this->InterpolateImage( imageIndex, this->m_Data[imageIndex] );
      const  typename Self::ReturnType upper = this->Evaluate();

      this->SetParameter( imageIndex, paramIndex, p0 - pStep );
      this->InterpolateImage( imageIndex, this->m_Data[imageIndex] );
      const  typename Self::ReturnType lower = this->Evaluate();

      this->m_Data[imageIndex] = tmp;
      this->SetParameter( imageIndex, paramIndex, p0 );

      if ( (upper > baseValue) || (lower > baseValue) )
	{
	g[param] = (upper - lower);
	}
      }
    }

  if ( this->m_ForceZeroSum )
    {
    this->ForceZeroSumGradient( g );
    }

  return baseValue;
}

template<class TXform>
bool
GroupwiseRegistrationRMIFunctional<TXform>
::Wiggle()
{
  bool wiggle = this->Superclass::Wiggle();

  if ( wiggle )
    {
    }
  
  return wiggle;
}

//@}

} // namespace cmtk

#include <cmtkAffineXform.h>
#include <cmtkSplineWarpXform.h>

template class cmtk::GroupwiseRegistrationRMIFunctional<cmtk::AffineXform>;
template class cmtk::GroupwiseRegistrationRMIFunctional<cmtk::SplineWarpXform>;
