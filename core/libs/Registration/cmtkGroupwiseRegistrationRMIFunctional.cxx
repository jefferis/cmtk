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

#include <Registration/cmtkGroupwiseRegistrationRMIFunctional.h>

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkEigenValuesSymmetricMatrix.h>

#include <System/cmtkThreadPool.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TXform>
GroupwiseRegistrationRMIFunctional<TXform>
::GroupwiseRegistrationRMIFunctional()
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
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfImages = this->m_ImageVector.size();

  this->m_CovarianceMatrix.Resize( numberOfImages ); // needs no reset
  this->m_TotalNumberOfSamples = 0;

  this->m_SumOfProductsMatrix.resize( numberOfImages * (1+numberOfImages) / 2 );
  std::fill( this->m_SumOfProductsMatrix.begin(), this->m_SumOfProductsMatrix.end(), 0 );

  this->m_SumsVector.resize( numberOfImages );
  std::fill( this->m_SumsVector.begin(), this->m_SumsVector.end(), 0 );

  this->m_ThreadSumOfProductsMatrix.resize( this->m_NumberOfThreads );
  this->m_ThreadSumsVector.resize( this->m_NumberOfThreads );

  std::vector<EvaluateThreadParameters> params( this->m_NumberOfTasks );
  for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfTasks; ++taskIdx )
    {
    params[taskIdx].thisObject = this;
    }
  
  if ( this->m_ProbabilisticSamples.size() )
    threadPool.Run( EvaluateProbabilisticThread, params );
  else
    threadPool.Run( EvaluateThread, params );
  
  return this->GetMetric( this->m_SumOfProductsMatrix, this->m_SumsVector, this->m_TotalNumberOfSamples, this->m_CovarianceMatrix );
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
      covarianceMatrix(i,j) = (sumOfProductsMatrix[midx] - ((1.0 * sumsVector[i] * sumsVector[j]) / totalNumberOfSamples)) / totalNumberOfSamples;
      }
    }
  
  const std::vector<typename Self::ReturnType> eigenvalues = EigenValuesSymmetricMatrix<typename Self::ReturnType>( covarianceMatrix ).GetEigenvalues();
  
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
    return -metric;
    }
  else
    {
    return -FLT_MAX;
    }
}

template<class TXform>
void
GroupwiseRegistrationRMIFunctional<TXform>
::EvaluateThread
( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  
  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const size_t numberOfImages = imagesTo - imagesFrom;
  
  const byte paddingValue = ThisConst->m_PaddingValue;

  SumsAndProductsVectorType& sumOfProductsMatrix = This->m_ThreadSumOfProductsMatrix[threadIdx];
  sumOfProductsMatrix.resize( numberOfImages * (1+numberOfImages) / 2 );
  std::fill( sumOfProductsMatrix.begin(), sumOfProductsMatrix.end(), 0 );
  
  SumsAndProductsVectorType& sumsVector = This->m_ThreadSumsVector[threadIdx];
  sumsVector.resize( numberOfImages );
  std::fill( sumsVector.begin(), sumsVector.end(), 0 );

  size_t totalNumberOfSamples = 0;

  const size_t numberOfPixels = ThisConst->m_TemplateNumberOfPixels;
  const size_t pixelsPerTask = 1+(numberOfPixels / taskCnt);
  const size_t pixelFrom = taskIdx * pixelsPerTask;
  const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerTask );

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
}

template<class TXform>
void
GroupwiseRegistrationRMIFunctional<TXform>
::EvaluateProbabilisticThread
( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  
  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const size_t numberOfImages = imagesTo - imagesFrom;
  
  const byte paddingValue = ThisConst->m_PaddingValue;

  SumsAndProductsVectorType& sumOfProductsMatrix = This->m_ThreadSumOfProductsMatrix[threadIdx];
  sumOfProductsMatrix.resize( numberOfImages * (1+numberOfImages) / 2 );
  std::fill( sumOfProductsMatrix.begin(), sumOfProductsMatrix.end(), 0 );
  
  SumsAndProductsVectorType& sumsVector = This->m_ThreadSumsVector[threadIdx];
  sumsVector.resize( numberOfImages );
  std::fill( sumsVector.begin(), sumsVector.end(), 0 );

  const size_t numberOfSamples = ThisConst->m_ProbabilisticSamples.size();
  const size_t samplesPerTask = 1+(numberOfSamples / taskCnt);
  const size_t sampleFrom = taskIdx * samplesPerTask;
  const size_t sampleTo = std::min( numberOfSamples, sampleFrom + samplesPerTask );

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

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkSplineWarpXform.h>

template class cmtk::GroupwiseRegistrationRMIFunctional<cmtk::AffineXform>;
template class cmtk::GroupwiseRegistrationRMIFunctional<cmtk::SplineWarpXform>;
