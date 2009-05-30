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

#include <cmtkCongealingFunctional.h>

#include <cmtkMathUtil.h>
#include <cmtkScaleHistogramValueTrait.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TXform,class THistogramBinType>
CongealingFunctional<TXform,THistogramBinType>::CongealingFunctional() 
  : m_NeedsUpdateStandardDeviationByPixel( true )
{
  this->SetNumberOfHistogramBins( this->m_HistogramBins );
}

template<class TXform,class THistogramBinType>
CongealingFunctional<TXform,THistogramBinType>::~CongealingFunctional()
{
  for ( size_t idx = 0; idx < this->m_HistogramKernel.size(); ++idx )
    if ( this->m_HistogramKernel[idx] )
      delete[] this->m_HistogramKernel[idx];
  this->m_HistogramKernel.clear();
}

template<class TXform,class THistogramBinType>
void
CongealingFunctional<TXform,THistogramBinType>
::SetNumberOfHistogramBins( const size_t numberOfHistogramBins )
{
  this->m_HistogramBins = numberOfHistogramBins;
  this->m_HistogramKernelRadiusMax = this->m_HistogramBins / 2;
  this->CreateGaussianKernels();

  this->Superclass::SetNumberOfHistogramBins( numberOfHistogramBins );
}

template<class TXform,class THistogramBinType>
void
CongealingFunctional<TXform,THistogramBinType>::CreateGaussianKernels()
{
  for ( size_t idx = 0; idx < this->m_HistogramKernel.size(); ++idx )
    if ( this->m_HistogramKernel[idx] )
      delete[] this->m_HistogramKernel[idx];

  this->m_HistogramKernel.resize( this->m_HistogramKernelRadiusMax+1 );
  this->m_HistogramKernelRadius.resize( this->m_HistogramKernelRadiusMax+1 );
  for ( size_t idx = 0; idx <= this->m_HistogramKernelRadiusMax; ++idx )
    {
    const size_t radius = idx + 1;
    const double sigma = idx;
    
    this->m_HistogramKernelRadius[idx] = radius;
    this->m_HistogramKernel[idx] = Memory::AllocateArray<HistogramBinType>( radius );
    
    if ( idx < 1.0 )
      {
      this->m_HistogramKernel[idx][0] = cmtk::ScaleHistogramValueTrait<HistogramBinType>::Scale( 1.0 );
      for ( size_t i = 1; i < radius; ++i )
	this->m_HistogramKernel[idx][i] = cmtk::ScaleHistogramValueTrait<HistogramBinType>::Scale( 0.0 );
      }
    else
      {
      const double normFactor = 1.0/(sqrt(2*M_PI) * sigma);
      for ( size_t i = 0; i < radius; ++i )
	{
	this->m_HistogramKernel[idx][i] = cmtk::ScaleHistogramValueTrait<HistogramBinType>::Scale( normFactor * exp( -MathUtil::Square( 1.0 * i / sigma ) / 2 ) );
	}
      }
    }
}

template<class TXform,class THistogramBinType>
void
CongealingFunctional<TXform,THistogramBinType>::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid,
  const int downsample,
  const bool useTemplateData )
{  
  this->Superclass::SetTemplateGrid( templateGrid, downsample, useTemplateData );
  this->m_NeedsUpdateStandardDeviationByPixel = true;
}

template<class TXform,class THistogramBinType>
typename CongealingFunctional<TXform,THistogramBinType>::ReturnType
CongealingFunctional<TXform,THistogramBinType>::Evaluate()
{
  if ( this->m_NeedsUpdateStandardDeviationByPixel )
    this->UpdateStandardDeviationByPixel();
  
  double entropy = 0;
  unsigned int count = 0;

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  ThreadParameterArray<Self,EvaluateThreadParameters> params( this, numberOfThreads );
  this->m_ThreadHistograms.resize( numberOfThreads );
  
  if ( this->m_ProbabilisticSamples.size() )
    params.RunInParallel( &EvaluateProbabilisticThread );
  else
    params.RunInParallel( &EvaluateThread );
  
  // gather partial entropies from threads
  for ( size_t thread = 0; thread < numberOfThreads; ++thread )
    {
    entropy += params[thread].m_Entropy;
    count += params[thread].m_Count;
    }

#ifdef CMTK_BUILD_MPI
  double partialEntropy = entropy;
  MPI::COMM_WORLD.Allreduce( &partialEntropy, &entropy, 1, MPI::DOUBLE, MPI::SUM );
  
  unsigned int partialCount = count;
  MPI::COMM_WORLD.Allreduce( &partialCount, &count, 1, MPI::UNSIGNED, MPI::SUM );
#endif
  
  if ( count )
    return static_cast<typename Self::ReturnType>( entropy / count );
  else
    return -FLT_MAX;
}

template<class TXform,class THistogramBinType>
CMTK_THREAD_RETURN_TYPE
CongealingFunctional<TXform,THistogramBinType>::EvaluateThread
( void *args )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  
  HistogramType& histogram = This->m_ThreadHistograms[threadID];
  histogram.Resize( ThisConst->m_HistogramBins + 2 * ThisConst->m_HistogramKernelRadiusMax, false /*reset*/ );
  
  double entropy = 0;
  unsigned int count = 0;

  const size_t numberOfPixels = ThisConst->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI  
  const size_t pixelsPerThread = 1+numberOfPixels / ( numberOfThreads * ThisConst->m_SizeMPI );
  const size_t pixelFrom = ( threadID + ThisConst->m_RankMPI * numberOfThreads ) * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerThread );
#else
  const size_t pixelsPerThread = 1+(numberOfPixels / numberOfThreads);
  const size_t pixelFrom = threadID * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerThread );
#endif

  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const byte paddingValue = ThisConst->m_PaddingValue;
  
  for ( size_t ofs = pixelFrom; ofs < pixelTo; ++ofs )
    {
    histogram.Reset();
    const size_t kernelIdx = ThisConst->m_StandardDeviationByPixel[ofs];
    const size_t kernelRadius = ThisConst->m_HistogramKernelRadius[kernelIdx];
    const HistogramBinType* kernel = ThisConst->m_HistogramKernel[kernelIdx];

    bool fullCount = true;
    if ( ThisConst->m_UseTemplateData )
      {
      const byte templateValue = ThisConst->m_TemplateData[ofs];
      if ( (fullCount = (templateValue != paddingValue )) )
	{
	histogram.AddWeightedSymmetricKernel( templateValue, kernelRadius, kernel );
	}
      }
    
    for ( size_t idx = imagesFrom; (idx < imagesTo) && fullCount; ++idx )
      {
      const byte value = ThisConst->m_Data[idx][ofs];
      if ( value != paddingValue )
	{
	histogram.AddWeightedSymmetricKernel( value, kernelRadius, kernel );
	}
      else
	{
	fullCount = false;
	}
      }

    if ( fullCount )
      {
      entropy -= histogram.GetEntropy();
      ++count;
      }
    }
  
  threadParameters->m_Entropy = entropy;
  threadParameters->m_Count = count;

  return CMTK_THREAD_RETURN_VALUE;
}

template<class TXform,class THistogramBinType>
CMTK_THREAD_RETURN_TYPE
CongealingFunctional<TXform,THistogramBinType>::EvaluateProbabilisticThread
( void *args )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  
  HistogramType& histogram = This->m_ThreadHistograms[threadID];
  histogram.Resize( ThisConst->m_HistogramBins + 2 * ThisConst->m_HistogramKernelRadiusMax, false /*reset*/ );

  double entropy = 0;
  unsigned int count = 0;

  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const byte paddingValue = ThisConst->m_PaddingValue;

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

  for ( size_t sample = sampleFrom; sample < sampleTo; ++sample )
    {
    histogram.Reset();
    bool fullCount = true;
    
    const size_t kernelIdx = ThisConst->m_StandardDeviationByPixel[sample];
    const size_t kernelRadius = ThisConst->m_HistogramKernelRadius[kernelIdx];
    const HistogramBinType* kernel = ThisConst->m_HistogramKernel[kernelIdx];

    if ( ThisConst->m_UseTemplateData )
      {
      const byte templateValue = ThisConst->m_TemplateData[sample];
      if ( (fullCount = (templateValue != paddingValue)) )
	{
	histogram.AddWeightedSymmetricKernel( templateValue, kernelRadius, kernel );
	}
      }
    
    for ( size_t idx = imagesFrom; (idx < imagesTo) && fullCount; ++idx )
      {
      const byte value = ThisConst->m_Data[idx][sample];
      if ( value != paddingValue )
	{
	histogram.AddWeightedSymmetricKernel( value, kernelRadius, kernel );
	}
      else
	{
	fullCount = false;
	}
      }
    
    if ( fullCount )
      {
      entropy -= histogram.GetEntropy();
      ++count;
      }	
    else
      {
      }
    }
    
  threadParameters->m_Entropy = entropy;
  threadParameters->m_Count = count;

  return CMTK_THREAD_RETURN_VALUE;
}

template<class TXform,class THistogramBinType>
bool
CongealingFunctional<TXform,THistogramBinType>
::Wiggle()
{
  bool wiggle = this->Superclass::Wiggle();

  if ( wiggle )
    {
    this->m_NeedsUpdateStandardDeviationByPixel = true;
    }
  
  return wiggle;
}

template<class TXform,class THistogramBinType>
void
CongealingFunctional<TXform,THistogramBinType>
::UpdateStandardDeviationByPixel()
{
  if ( this->m_ProbabilisticSamples.size() )
    {
    const size_t numberOfSamples = this->m_ProbabilisticSamples.size();
#ifdef CMTK_BUILD_MPI
    const size_t samplesPerNode = (numberOfSamples+this->m_SizeMPI-1) / this->m_SizeMPI;
    this->m_StandardDeviationByPixel.resize( samplesPerNode * this->m_SizeMPI );
    this->m_StandardDeviationByPixelMPI.resize( samplesPerNode );
#else
    this->m_StandardDeviationByPixel.resize( numberOfSamples );
#endif
    }
  else
    {
    const size_t numberOfPixels = this->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI
    const size_t pixelsPerNode = (numberOfPixels+this->m_SizeMPI-1) / this->m_SizeMPI;
    this->m_StandardDeviationByPixel.resize( pixelsPerNode * this->m_SizeMPI );
    this->m_StandardDeviationByPixelMPI.resize( pixelsPerNode );
#else
    this->m_StandardDeviationByPixel.resize( numberOfPixels );
#endif
    }

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  ThreadParameterArray<Self,ThreadParametersType> params( this, numberOfThreads );
  params.RunInParallel( &UpdateStandardDeviationByPixelThreadFunc );

#ifdef CMTK_BUILD_MPI
  MPI::COMM_WORLD.Allgather
    ( &this->m_StandardDeviationByPixelMPI[0], this->m_StandardDeviationByPixelMPI.size(), MPI::CHAR,
      &this->m_StandardDeviationByPixel[0], this->m_StandardDeviationByPixelMPI.size(), MPI::CHAR );
#endif

  this->m_NeedsUpdateStandardDeviationByPixel = false;
}


template<class TXform,class THistogramBinType>
CMTK_THREAD_RETURN_TYPE
CongealingFunctional<TXform,THistogramBinType>
::UpdateStandardDeviationByPixelThreadFunc( void *args )
{
  ThreadParametersType* threadParameters = static_cast<ThreadParametersType*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;

  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const byte paddingValue = ThisConst->m_PaddingValue;
  
  if ( ThisConst->m_ProbabilisticSamples.size() )
    {
    const size_t numberOfSamples = ThisConst->m_ProbabilisticSamples.size();
#ifdef CMTK_BUILD_MPI
    const size_t samplesPerNode = (numberOfSamples+ThisConst->m_SizeMPI-1) / ThisConst->m_SizeMPI;
    const size_t samplesPerThread = 1 + (samplesPerNode / numberOfThreads);
    const size_t sampleFromNode = ThisConst->m_RankMPI * samplesPerNode;
    const size_t sampleFrom = sampleFromNode + threadID * samplesPerThread;
    const size_t sampleTo = std::min( numberOfSamples, std::min( sampleFromNode + samplesPerNode, 
								 sampleFrom + samplesPerThread ) );
    size_t mpiSmpl = threadID * samplesPerThread;
#else
    const size_t samplesPerThread = 1 + (numberOfSamples / numberOfThreads );
    const size_t sampleFrom = threadID * samplesPerThread;
    const size_t sampleTo = std::min( numberOfSamples, sampleFrom + samplesPerThread );
#endif

    for ( size_t smpl = sampleFrom; smpl < sampleTo; ++smpl )
      {
      float sum = 0, sumsq = 0;
      unsigned int count = 0;

      if ( ThisConst->m_UseTemplateData )
	{
	const byte templateValue = ThisConst->m_TemplateData[smpl];
	if ( templateValue != paddingValue )
	  {
	  sum += templateValue;
	  sumsq += templateValue * templateValue;
	  ++count;
	  }
	}

      for ( size_t idx = imagesFrom; idx < imagesTo; ++idx )
	{
	const byte value = ThisConst->m_Data[idx][smpl];
	if ( value != paddingValue )
	  {
	  const double data = static_cast<double>( value );
	  sum += data;
	  sumsq += data * data;
	  ++count;
	  }
	}
      
      if ( count )
	{
	const double mu = sum / count;
	const byte sdev = std::min<byte>( ThisConst->m_HistogramKernelRadiusMax, (byte)(sqrt(( count * mu * mu - 2 * mu * sum + sumsq ) / (count-1)) ) );

#ifdef CMTK_BUILD_MPI
	This->m_StandardDeviationByPixelMPI[mpiSmpl++] = sdev;
#else
	This->m_StandardDeviationByPixel[smpl] = sdev;
#endif
	}
      else
	{
	This->m_StandardDeviationByPixel[smpl] = 0;
	}
      }
    }
  else
    {
    const size_t numberOfPixels = ThisConst->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI
    const size_t pixelsPerNode = (numberOfPixels+ThisConst->m_SizeMPI-1) / ThisConst->m_SizeMPI;
    const size_t pixelsPerThread = 1 + (pixelsPerNode / numberOfThreads);
    const size_t pixelFromNode = ThisConst->m_RankMPI * pixelsPerNode;
    const size_t pixelFrom = pixelFromNode + threadID * pixelsPerThread;
    const size_t pixelTo = std::min( numberOfPixels, std::min( pixelFromNode + pixelsPerNode, 
								 pixelFrom + pixelsPerThread ) );
    size_t mpiPx = threadID * pixelsPerThread;
#else
    const size_t pixelsPerThread = 1 + (numberOfPixels / numberOfThreads);
    const size_t pixelFrom = threadID * pixelsPerThread;
    const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerThread );
#endif
    for ( size_t px = pixelFrom; px < pixelTo; ++px )
      {
      double sum = 0, sumsq = 0;
      unsigned int count = 0;

      if ( ThisConst->m_UseTemplateData )
	{
	const byte templateValue = ThisConst->m_TemplateData[px];
	if ( templateValue != paddingValue )
	  {
	  sum += templateValue;
	  sumsq += templateValue * templateValue;
	  ++count;
	  }
	}

      for ( size_t idx = imagesFrom; idx < imagesTo; ++idx )
	{
	const byte value = ThisConst->m_Data[idx][px];
	if ( value != paddingValue )
	  {
	  const double data = static_cast<double>( value );
	  sum += data;
	  sumsq += data * data;
	  ++count;
	  }
	}
      
      if ( count )
	{
	const double mu = sum / count;
	const byte sdev = std::min<byte>( ThisConst->m_HistogramKernelRadiusMax, (byte)(sqrt(( count * mu * mu - 2 * mu * sum + sumsq ) / (count-1)) ) );
#ifdef CMTK_BUILD_MPI
	This->m_StandardDeviationByPixelMPI[mpiPx++] = sdev;
#else
	This->m_StandardDeviationByPixel[px] = sdev;
#endif
	}
      else
	{
	This->m_StandardDeviationByPixel[px] = 0;
	}
      }
    }

  return CMTK_THREAD_RETURN_VALUE;
}

//@}

} // namespace cmtk

#include <cmtkAffineXform.h>
#include <cmtkSplineWarpXform.h>

template class cmtk::CongealingFunctional<cmtk::AffineXform>;
template class cmtk::CongealingFunctional<cmtk::SplineWarpXform>;
