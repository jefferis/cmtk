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

#include "cmtkCongealingFunctional.h"

#include "Base/cmtkMathUtil.h"
#include "Registration/cmtkScaleHistogramValueTrait.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TXform>
CongealingFunctional<TXform>::CongealingFunctional() 
  : m_NeedsUpdateStandardDeviationByPixel( true )
{
  this->SetNumberOfHistogramBins( this->m_HistogramBins );
}

template<class TXform>
CongealingFunctional<TXform>::~CongealingFunctional()
{
  for ( size_t idx = 0; idx < this->m_HistogramKernel.size(); ++idx )
    if ( this->m_HistogramKernel[idx] )
      delete[] this->m_HistogramKernel[idx];
  this->m_HistogramKernel.clear();
}

template<class TXform>
void
CongealingFunctional<TXform>
::SetNumberOfHistogramBins( const size_t numberOfHistogramBins )
{
  this->m_HistogramBins = numberOfHistogramBins;
  this->m_HistogramKernelRadiusMax = this->m_HistogramBins / 2;
  this->CreateGaussianKernels();

  this->Superclass::SetNumberOfHistogramBins( numberOfHistogramBins );
}

template<class TXform>
void
CongealingFunctional<TXform>::CreateGaussianKernels()
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

template<class TXform>
void
CongealingFunctional<TXform>::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid,
  const int downsample,
  const bool useTemplateData )
{  
  this->Superclass::SetTemplateGrid( templateGrid, downsample, useTemplateData );
  this->m_NeedsUpdateStandardDeviationByPixel = true;
}

template<class TXform>
typename CongealingFunctional<TXform>::ReturnType
CongealingFunctional<TXform>::Evaluate()
{
  if ( this->m_NeedsUpdateStandardDeviationByPixel )
    this->UpdateStandardDeviationByPixel();
  
  double entropy = 0;
  unsigned int count = 0;

  this->m_ThreadHistograms.resize( this->m_NumberOfThreads );

  std::vector<EvaluateThreadParameters> params( this->m_NumberOfTasks );
  for ( size_t idx = 0; idx < this->m_NumberOfTasks; ++idx )
    params[idx].thisObject = this;
  
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  if ( this->m_ProbabilisticSamples.size() )
    threadPool.Run( Self::EvaluateProbabilisticThread, params );
  else
    threadPool.Run( Self::EvaluateThread, params );
  
  // gather partial entropies from tasks
  for ( size_t task = 0; task < this->m_NumberOfTasks; ++task )
    {
    entropy += params[task].m_Entropy;
    count += params[task].m_Count;
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

template<class TXform>
void
CongealingFunctional<TXform>::EvaluateThread
( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  
  HistogramType& histogram = This->m_ThreadHistograms[threadIdx];
  histogram.Resize( ThisConst->m_HistogramBins + 2 * ThisConst->m_HistogramKernelRadiusMax, false /*reset*/ );
  
  double entropy = 0;
  unsigned int count = 0;

  const size_t numberOfPixels = ThisConst->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI  
  const size_t pixelsPerThread = 1+numberOfPixels / ( taskCnt * ThisConst->m_SizeMPI );
  const size_t pixelFrom = ( taskIdx + ThisConst->m_RankMPI * taskCnt ) * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerThread );
#else
  const size_t pixelsPerThread = 1+(numberOfPixels / taskCnt);
  const size_t pixelFrom = taskIdx * pixelsPerThread;
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
}

template<class TXform>
void
CongealingFunctional<TXform>::EvaluateProbabilisticThread
( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  
  HistogramType& histogram = This->m_ThreadHistograms[threadIdx];
  histogram.Resize( ThisConst->m_HistogramBins + 2 * ThisConst->m_HistogramKernelRadiusMax, false /*reset*/ );

  double entropy = 0;
  unsigned int count = 0;

  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const byte paddingValue = ThisConst->m_PaddingValue;

  const size_t numberOfSamples = ThisConst->m_ProbabilisticSamples.size();
#ifdef CMTK_BUILD_MPI  
  const size_t samplesPerThread = numberOfSamples / ( taskCnt * ThisConst->m_SizeMPI );
  const size_t sampleFrom = ( taskIdx + ThisConst->m_RankMPI * taskCnt ) * samplesPerThread;
  const size_t sampleTo = std::min( numberOfSamples, sampleFrom + samplesPerThread );
#else
  const size_t samplesPerThread = numberOfSamples / taskCnt;
  const size_t sampleFrom = taskIdx * samplesPerThread;
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
}

template<class TXform>
bool
CongealingFunctional<TXform>
::Wiggle()
{
  bool wiggle = this->Superclass::Wiggle();

  if ( wiggle )
    {
    this->m_NeedsUpdateStandardDeviationByPixel = true;
    }
  
  return wiggle;
}

template<class TXform>
void
CongealingFunctional<TXform>
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

  std::vector<ThreadParametersType> params( this->m_NumberOfTasks );
  for ( size_t idx = 0; idx < this->m_NumberOfTasks; ++idx )
    params[idx].thisObject = this;

  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  threadPool.Run( Self::UpdateStandardDeviationByPixelThreadFunc, params );
  
#ifdef CMTK_BUILD_MPI
  MPI::COMM_WORLD.Allgather( &this->m_StandardDeviationByPixelMPI[0], this->m_StandardDeviationByPixelMPI.size(), 
			     MPI::CHAR, &this->m_StandardDeviationByPixel[0], this->m_StandardDeviationByPixelMPI.size(), MPI::CHAR );
#endif
  
  this->m_NeedsUpdateStandardDeviationByPixel = false;
}

template<class TXform>
void
CongealingFunctional<TXform>
::UpdateStandardDeviationByPixelThreadFunc
( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  ThreadParametersType* taskParameters = static_cast<ThreadParametersType*>( args );
  
  Self* This = taskParameters->thisObject;
  const Self* ThisConst = taskParameters->thisObject;

  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const byte paddingValue = ThisConst->m_PaddingValue;
  
  if ( ThisConst->m_ProbabilisticSamples.size() )
    {
    const size_t numberOfSamples = ThisConst->m_ProbabilisticSamples.size();
#ifdef CMTK_BUILD_MPI
    const size_t samplesPerNode = (numberOfSamples+ThisConst->m_SizeMPI-1) / ThisConst->m_SizeMPI;
    const size_t samplesPerTask = 1 + (samplesPerNode / taskCnt);
    const size_t sampleFromNode = ThisConst->m_RankMPI * samplesPerNode;
    const size_t sampleFrom = sampleFromNode + taskIdx * samplesPerTask;
    const size_t sampleTo = std::min( numberOfSamples, std::min( sampleFromNode + samplesPerNode, sampleFrom + samplesPerTask ) );
    size_t mpiSmpl = taskIdx * samplesPerTask;
#else
    const size_t samplesPerTask = 1 + (numberOfSamples / taskCnt );
    const size_t sampleFrom = taskIdx * samplesPerTask;
    const size_t sampleTo = std::min( numberOfSamples, sampleFrom + samplesPerTask );
#endif

    for ( size_t smpl = sampleFrom; smpl < sampleTo; ++smpl )
      {
      double sum = 0, sumsq = 0;
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
    const size_t pixelsPerTask = 1 + (pixelsPerNode / taskCnt);
    const size_t pixelFromNode = ThisConst->m_RankMPI * pixelsPerNode;
    const size_t pixelFrom = pixelFromNode + taskIdx * pixelsPerTask;
    const size_t pixelTo = std::min( numberOfPixels, std::min( pixelFromNode + pixelsPerNode, pixelFrom + pixelsPerTask ) );
    size_t mpiPx = taskIdx * pixelsPerTask;
#else
    const size_t pixelsPerTask = 1 + (numberOfPixels / taskCnt);
    const size_t pixelFrom = taskIdx * pixelsPerTask;
    const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerTask );
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
}

//@}

} // namespace cmtk

#include "Base/cmtkAffineXform.h"
#include "Base/cmtkSplineWarpXform.h"

template class cmtk::CongealingFunctional<cmtk::AffineXform>;
template class cmtk::CongealingFunctional<cmtk::SplineWarpXform>;
