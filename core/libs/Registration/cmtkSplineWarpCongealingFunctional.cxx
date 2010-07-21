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

#include "cmtkSplineWarpCongealingFunctional.h"

#include "Base/cmtkMathUtil.h"
#include "Base/cmtkMatrix.h"

#include "System/cmtkThreadPool.h"
#include "System/cmtkThreadParameterArray.h"

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

void
SplineWarpCongealingFunctional
::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid, 
  const int downsample,
  const bool useTemplateData )
{
  this->Superclass::SetTemplateGrid( templateGrid, downsample, useTemplateData );
  // clear thread storage because we need to re-initialize these.
  this->m_StaticThreadStorage.resize(0);
}

void
SplineWarpCongealingFunctional
::InitializeXformsFromAffine
( const Types::Coordinate gridSpacing, std::vector<AffineXform::SmartPtr> initialAffineXformsVector, const bool exactSpacing )
{
  this->Superclass::InitializeXformsFromAffine( gridSpacing, initialAffineXformsVector, exactSpacing );

  // clear thread storage because we need to re-initialize these.
  this->m_StaticThreadStorage.resize(0);
}

void
SplineWarpCongealingFunctional
::RefineTransformationGrids()
{
  this->Superclass::RefineTransformationGrids();
  // clear thread storage because we need to re-initialize these.
  this->m_StaticThreadStorage.resize(0);
}

void
SplineWarpCongealingFunctional
::UpdateStandardDeviationByPixel()
{
  this->Superclass::UpdateStandardDeviationByPixel();
  this->UpdateActiveControlPoints();
  this->UpdateParamStepArray();
}

void
SplineWarpCongealingFunctional
::UpdateActiveControlPoints()
{
  if ( this->m_DeactivateUninformativeMode )
    {
    const size_t numberOfControlPoints = this->m_VolumeOfInfluenceArray.size();
    
    if ( numberOfControlPoints )
      {
      this->m_ActiveControlPointFlags.resize( numberOfControlPoints );
      this->m_NumberOfActiveControlPoints = 0;
      
      const Vector3D templateFrom( this->m_TemplateGrid->m_Offset );
      const Vector3D templateTo( this->m_TemplateGrid->m_Offset + this->m_TemplateGrid->Size );
      Vector3D fromVOI, toVOI;
      
      std::vector<DataGrid::RegionType>::const_iterator voi = this->m_VolumeOfInfluenceArray.begin();
      for ( size_t cp = 0; cp < numberOfControlPoints; ++cp, ++voi )
	{
	bool active = false;
	for ( int z = voi->From()[2]; (z < voi->To()[2]) && !active; ++z ) 
	  {
	  for ( int y = voi->From()[1]; (y < voi->To()[1]) && !active; ++y )
	    {
	    size_t ofs = this->m_TemplateGrid->GetOffsetFromIndex( voi->From()[0], y, z );
	    for ( int x = voi->From()[0]; (x < voi->To()[0])  && !active; ++x, ++ofs )
	      {
	      if ( this->m_StandardDeviationByPixel[ofs] > 0 )
		{
		active = true;
		}
	      }
	    }
	  }
	this->m_ActiveControlPointFlags[cp] = active;
	if ( active ) ++this->m_NumberOfActiveControlPoints;
	}

      StdErr << "Enabled " << this->m_NumberOfActiveControlPoints << "/" << this->m_ParametersPerXform / 3 << " control points.\n";
      }
    }
  else
    {
    this->m_NumberOfActiveControlPoints = this->m_VolumeOfInfluenceArray.size();
    }
}

SplineWarpCongealingFunctional::ReturnType
SplineWarpCongealingFunctional
::Evaluate()
{
  if ( this->m_NeedsUpdateStandardDeviationByPixel )
    this->UpdateStandardDeviationByPixel();

  const size_t numberOfPixels = this->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI
  const size_t pixelsPerNode = (numberOfPixels+this->m_SizeMPI-1) / this->m_SizeMPI;
  this->m_EntropyByPixel.resize( pixelsPerNode * this->m_SizeMPI );
  this->m_EntropyByPixelMPI.resize( pixelsPerNode );
#else
  this->m_EntropyByPixel.resize( numberOfPixels );
#endif

  double entropy = 0;
  unsigned int count = 0;
  
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfThreads = threadPool.GetNumberOfThreads();
  this->m_ThreadHistograms.resize( numberOfThreads );

  std::vector< EvaluateThreadParameters> params( numberOfThreads );
  for ( size_t taskIdx = 0; taskIdx < numberOfThreads; ++taskIdx )
    {
    params[taskIdx].thisObject = this;
    }
  threadPool.Run( EvaluateThread, params );
  
  // gather partial entropies from threads
  for ( size_t taskIdx = 0; taskIdx < numberOfThreads; ++taskIdx )
    {
    entropy += params[taskIdx].m_Entropy;
    count += params[taskIdx].m_Count;
    }
  
#ifdef CMTK_BUILD_MPI
  double partialEntropy = entropy;
  MPI::COMM_WORLD.Allreduce( &partialEntropy, &entropy, 1, MPI::DOUBLE, MPI::SUM );
  
  unsigned int partialCount = count;
  MPI::COMM_WORLD.Allreduce( &partialCount, &count, 1, MPI::UNSIGNED, MPI::SUM );

  const size_t totalSize = sizeof( this->m_EntropyByPixelMPI[0] ) * this->m_EntropyByPixelMPI.size();
  MPI::COMM_WORLD.Allgather( &this->m_EntropyByPixelMPI[0], totalSize, MPI::CHAR, &this->m_EntropyByPixel[0], totalSize, MPI::CHAR );  
#endif
  
  if ( count )
    {
    const double result = entropy / count;
    double constraint = 0;
    if ( this->m_JacobianConstraintWeight > 0 )
      {
      for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
	{
	constraint += dynamic_cast<const SplineWarpXform*>( this->m_XformVector[i].GetPtr() )->GetJacobianConstraint();
	}
      }
    return result - this->m_JacobianConstraintWeight * constraint;
    }
  else
    return -FLT_MAX;
}

void
SplineWarpCongealingFunctional
::EvaluateThread
( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  
  HistogramType& histogram = This->m_ThreadHistograms[threadIdx];
  histogram.Resize( ThisConst->m_HistogramBins + 2 * ThisConst->m_HistogramKernelRadiusMax, false /*reset*/ );

  double totalEntropy = 0;
  size_t count = 0;

  const size_t numberOfPixels = ThisConst->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI  
  const size_t pixelsPerNode = (numberOfPixels+ThisConst->m_SizeMPI-1) / ThisConst->m_SizeMPI;
  const size_t pixelsPerThread = 1 + (pixelsPerNode / taskCnt);
  const size_t pixelFromNode = ThisConst->m_RankMPI * pixelsPerNode;
  const size_t pixelFrom = pixelFromNode + taskIdx * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, std::min( pixelFromNode + pixelsPerNode, pixelFrom + pixelsPerThread ) );
  size_t mpiOfs = taskIdx * pixelsPerThread;
#else
  const size_t pixelsPerThread = (numberOfPixels / taskCnt);
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
      if ( (fullCount = (templateValue != paddingValue)) )
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
      const double entropy = histogram.GetEntropy();
#ifdef CMTK_BUILD_MPI
      This->m_EntropyByPixelMPI[mpiOfs++] = entropy;
#else
      This->m_EntropyByPixel[ofs] = entropy;
#endif
      totalEntropy -= entropy;
      ++count;
      }
    else
      {
#ifdef CMTK_BUILD_MPI
      This->m_EntropyByPixelMPI[mpiOfs++] = 0;
#else
      This->m_EntropyByPixel[ofs] = 0;
#endif
      }
    }
  
  threadParameters->m_Entropy = totalEntropy;
  threadParameters->m_Count = count;
}

//@}

} // namespace cmtk

#ifdef CMTK_BUILD_MPI
#  include "cmtkSplineWarpCongealingFunctionalMPI.txx"
#else
#  include "cmtkSplineWarpCongealingFunctional.txx"
#endif
