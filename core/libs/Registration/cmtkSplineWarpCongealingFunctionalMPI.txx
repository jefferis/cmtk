/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

namespace
cmtk
{

/** \addtogroup Registration */
//@{

SplineWarpCongealingFunctional::ReturnType
SplineWarpCongealingFunctional
::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
  const Self::ReturnType baseValue = this->EvaluateAt( v );

  const size_t numberOfXforms = this->m_XformVector.size();
  const size_t numberOfControlPoints = this->m_ParametersPerXform / 3;
  // use more chunks than nodes to get more fine-grained parallelism
  const size_t controlPointsPerNode = 1 + static_cast<size_t>( numberOfControlPoints / (4 * this->m_SizeMPI) );
  // assign less work per chunk to root node to leave time to check for incoming results and
  // keep worker nodes maximally busy
  const size_t controlPointsRootNode = 1 + controlPointsPerNode / 4;

  MPI::Status status;

  size_t resultsToReceive = 0;
  std::vector<size_t> taskByNode( this->m_SizeMPI, -1 );

  const size_t msgBufferSize = 3 * sizeof( Types::Coordinate ) * controlPointsPerNode * numberOfXforms;
  std::vector<char> msgBuffer( msgBufferSize );
  char* msgBufferPtr = &msgBuffer[0];

  // Allocate thread storage
  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  if ( this->m_StaticThreadStorage.size() != numberOfThreads )
    {
    this->m_StaticThreadStorage.resize( numberOfThreads );
    for ( size_t thread = 0; thread < numberOfThreads; ++thread )
      {
      this->m_StaticThreadStorage[thread].Initialize( this );
      }
    }
  else
    {
    for ( size_t thread = 0; thread < numberOfThreads; ++thread )
      {
      this->m_StaticThreadStorage[thread].m_NeedToCopyXformParameters = true;
      }
    }

  // Create thread parameter array and start threads asynchronously
  ThreadParameterArray<Self,EvaluateLocalGradientThreadParameters> threadParams( this, numberOfThreads );
  threadParams.RunInParallelAsynchronous( EvaluateLocalGradientThreadFunc );
  for ( size_t thread = 0; thread < numberOfThreads; ++thread )
    {
    threadParams[thread].m_Step = step;
    threadParams[thread].m_Gradient = (Types::Coordinate*)msgBufferPtr;
    }
  
  bool allJobsSent = false;
  size_t firstIndexToCompute = 0;
  while ( (firstIndexToCompute < numberOfControlPoints) && ! allJobsSent )
    {
    if ( this->m_RankMPI == 0 )
      {      
      // send out job instructions or termination signal
      for ( int dest = 1; (dest < this->m_SizeMPI) && !allJobsSent; ++dest )
	{
	if ( taskByNode[dest] == (size_t)-1 )
	  {
	  if ( firstIndexToCompute >= numberOfControlPoints )
	    {
	    allJobsSent = true;
	    }
	  else
	    {
	    taskByNode[dest] = firstIndexToCompute;
	    MPI::COMM_WORLD.Ssend( &firstIndexToCompute, sizeof( firstIndexToCompute ), MPI::CHAR, dest, MESSAGE_TAG_COMPUTE );
	    ++resultsToReceive;

	    firstIndexToCompute += controlPointsPerNode;
	    }
	  }
	}
      }
    else
      {
      // receive job instructions      
      MPI::COMM_WORLD.Recv( msgBufferPtr, msgBufferSize, MPI::CHAR, 0 /*root*/ , MPI::ANY_TAG, status );
      if ( status.Get_tag() == MESSAGE_TAG_FINISHED )
	{
	allJobsSent = true;
	}
      else
	{
	int position = 0;
	MPI::CHAR.Unpack( msgBufferPtr, msgBufferSize, &firstIndexToCompute, sizeof( firstIndexToCompute ), position, MPI::COMM_WORLD );
	}
      }
    
    if ( (firstIndexToCompute < numberOfControlPoints ) && ! allJobsSent )
      {
      const size_t cpsToCompute = std::min( this->m_RankMPI ? controlPointsPerNode : controlPointsRootNode, (int)numberOfControlPoints - firstIndexToCompute );
      // determine last control point processed
      this->m_ControlPointIndexNext = firstIndexToCompute;
      this->m_ControlPointIndexLast = firstIndexToCompute + cpsToCompute;

      for ( size_t thread = 0; thread < numberOfThreads; ++thread )
	{
	threadParams[thread].m_FirstIndexToCompute = firstIndexToCompute;
	// Post "work ready" semaphore for each thread
	this->m_ThreadWorkSemaphore.Post();
	}

      // Wait for "work done" semaphore from each thread
      for ( size_t nThreads = 0; nThreads < numberOfThreads; ++nThreads )
	{
	this->m_ThreadReadySemaphore.Wait();
	}
      
      // on root node, simply sort results into gradient fields
      if ( this->m_RankMPI == 0 )
	{
	// sort computed values into gradient vector
	this->ReorderGradientComponents( g.Elements, (Types::Coordinate*)msgBufferPtr, firstIndexToCompute, this->m_ControlPointIndexLast );
	firstIndexToCompute += cpsToCompute;
	}
      else
	{
	// send results
	MPI::COMM_WORLD.Ssend( msgBufferPtr, 3 * numberOfXforms * cpsToCompute * sizeof( Types::Coordinate ), MPI::CHAR, 0, MESSAGE_TAG_RESULTS );
	}
      }
    
    // are we still waiting for results from worker nodes?
    // are the worker nodes still trying to send data?
    if ( this->m_RankMPI == 0 )
      {
      // receive results:
      // is a result ready to be received?
      while ( resultsToReceive && MPI::COMM_WORLD.Iprobe( MPI::ANY_SOURCE, MESSAGE_TAG_RESULTS, status ) )
	{
	// receive result
	MPI::COMM_WORLD.Recv( msgBufferPtr, msgBufferSize, MPI::CHAR, status.Get_source(), MESSAGE_TAG_RESULTS, status );
		
	// what control point did we send this node for processing?
	const int fromNode = status.Get_source();
	const int receivedJob = taskByNode[fromNode];
	
	--resultsToReceive;
	taskByNode[fromNode] = -1;
	
	// determine last control point processed
	const size_t maxCpIndex = std::min<size_t>( receivedJob + controlPointsPerNode, numberOfControlPoints );
	this->ReorderGradientComponents( g.Elements, (Types::Coordinate*)msgBufferPtr, receivedJob, maxCpIndex );
	}
      }
    }

  // collect remaining results from worker nodes
  while ( resultsToReceive )
    {
    // receive result
    MPI::COMM_WORLD.Recv( msgBufferPtr, msgBufferSize, MPI::CHAR, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
    
    // what control point did we send this node for processing?
    const int fromNode = status.Get_source();
    const int receivedJob = taskByNode[fromNode];

    --resultsToReceive;
    taskByNode[fromNode] = -1;
    
    // determine last control point processed
    const size_t maxCpIndex = std::min<size_t>( receivedJob + controlPointsPerNode, numberOfControlPoints );
    this->ReorderGradientComponents( g.Elements, (Types::Coordinate*)msgBufferPtr, receivedJob, maxCpIndex );
    }  
 
  if ( this->m_RankMPI == 0 )
    { 
    // tell all worker nodes to finish their computation loops
    for ( int workerNode = 1; workerNode < this->m_SizeMPI; ++workerNode )
      {
      MPI::COMM_WORLD.Ssend( msgBufferPtr, 1, MPI::CHAR, workerNode, MESSAGE_TAG_FINISHED );
      }
    }
  
  // Distribute gradient estimate
  MPI::COMM_WORLD.Bcast( g.Elements, g.Dim * sizeof( g.Elements[0] ), MPI::CHAR, 0 /*root*/ );
  
  if ( this->m_PartialGradientMode )
    {
    const Types::Coordinate gthresh = g.MaxNorm() * this->m_PartialGradientThreshold;
    for ( size_t param = 0; param < g.Dim; ++param )
      {
      if ( fabs( g[param] ) < gthresh )
	{	
	g[param] = this->m_ParamStepArray[param] = 0.0;
	}
      }
    }

  if ( this->m_ForceZeroSum )
    {
    this->ForceZeroSumGradient( g );
    }

  return baseValue;
}

void
SplineWarpCongealingFunctional
::ReorderGradientComponents
( Types::Coordinate *const dst, const Types::Coordinate* src, const size_t fromCpIdx, const size_t toCpIdx )
{
  const size_t numberOfParameters = this->ParamVectorDim();
  
  // sort computed values into gradient vector
  size_t currentParameter = 0;
  for ( size_t cpIndex = fromCpIdx; cpIndex < toCpIdx; ++cpIndex )
    {
    for ( size_t cparam = 3*cpIndex; cparam < numberOfParameters; cparam += this->m_ParametersPerXform )
      {
      for ( size_t dim = 0; dim < 3; ++dim, ++currentParameter )
	{
	dst[cparam+dim] = src[currentParameter];
	}
      }
    }
}

CMTK_THREAD_RETURN_TYPE
SplineWarpCongealingFunctional
::EvaluateLocalGradientThreadFunc
( void* args )
{
  EvaluateLocalGradientThreadParameters* threadParameters = static_cast<EvaluateLocalGradientThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = This;
  const size_t threadID = threadParameters->ThisThreadIndex;

  const size_t numberOfXforms = ThisConst->m_XformVector.size();
  const size_t parametersPerXform = ThisConst->m_ParametersPerXform;
  const size_t paramVectorDim = ThisConst->ParamVectorDim();
  
  const byte paddingValue = ThisConst->m_PaddingValue;
  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const size_t numberOfImages = imagesTo - imagesFrom;
  const size_t numberOfImagesIncludingTemplate = ThisConst->m_UseTemplateData ? numberOfImages+1 : numberOfImages;

  Self::StaticThreadStorage* threadStorage = &(This->m_StaticThreadStorage[threadID]);
  if ( threadStorage->m_NeedToCopyXformParameters )
    {
    for ( size_t xi = 0; xi < numberOfXforms; ++xi )
      {
      threadStorage->m_Xforms[xi]->CopyParamVector( This->GetXformByIndex(xi) );
      threadStorage->m_NeedToCopyXformParameters = false;
      }
    }

  const UniformVolume* templateGrid = ThisConst->m_TemplateGrid;
  const size_t numberOfControlPoints = ThisConst->m_ParametersPerXform / 3;
  
  size_t cpIndex;
  
  while ( true )
    {
    This->m_ThreadWorkSemaphore.Wait();
    
    This->m_ControlPointIndexLock.Lock();
    while ( (cpIndex = This->m_ControlPointIndexNext) < This->m_ControlPointIndexLast )
      {
      ++This->m_ControlPointIndexNext;
      This->m_ControlPointIndexLock.Unlock();
      
      if ( !(cpIndex % 1000) )
	{
	std::cerr << cpIndex << " / " << numberOfControlPoints << "\r";
	}    
      
      std::vector<DataGrid::RegionType>::const_iterator voi = ThisConst->m_VolumeOfInfluenceArray.begin() + cpIndex;
      const size_t pixelsPerLineVOI = (voi->To()[0]-voi->From()[0]);
      
      std::fill( threadStorage->m_FPlus.begin(), threadStorage->m_FPlus.end(), 0 );
      std::fill( threadStorage->m_FMinus.begin(), threadStorage->m_FMinus.end(), 0 );
      std::fill( threadStorage->m_CountByParameterPlus.begin(), threadStorage->m_CountByParameterPlus.end(), 0 );
      std::fill( threadStorage->m_CountByParameterMinus.begin(), threadStorage->m_CountByParameterMinus.end(), 0 );
      
      for ( int z = voi->From()[2]; (z < voi->To()[2]); ++z ) 
	{
	for ( int y = voi->From()[1]; (y < voi->To()[1]); ++y )
	  {
	  // evaluate baseline entropies for this row
	  const size_t rowofs = templateGrid->GetOffsetFromIndex( voi->From()[0], y, z );
	  
	  size_t ofs = rowofs;
	  for ( size_t x = 0; x < pixelsPerLineVOI; ++x, ++ofs )
	    {
	    threadStorage->m_Histogram[x].Reset();
	    threadStorage->m_Count[x] = 0;
	  
	    const size_t kernelIdx = ThisConst->m_StandardDeviationByPixel[ofs];
	    const size_t kernelRadius = ThisConst->m_HistogramKernelRadius[kernelIdx];
	    const HistogramBinType* kernel = ThisConst->m_HistogramKernel[kernelIdx];
	  
	    if ( ThisConst->m_UseTemplateData )
	      {
	      const byte templateValue = ThisConst->m_TemplateData[ofs];
	      if ( templateValue != paddingValue )
		{
		threadStorage->m_Histogram[x].AddWeightedSymmetricKernel( templateValue, kernelRadius, kernel );
		++threadStorage->m_Count[x];
		}
	      }
	    
	    for ( size_t img = imagesFrom; img < imagesTo; ++img )
	      { 
	      const byte dataThisPixel = ThisConst->m_Data[img][ofs];
	      if ( dataThisPixel != paddingValue )
		{
		threadStorage->m_Histogram[x].AddWeightedSymmetricKernel( dataThisPixel, kernelRadius, kernel );
		++threadStorage->m_Count[x];
		}
	      }
	    }
	
	  size_t cparam = 3 * cpIndex;
	  size_t currentParameter = 0;
	  for ( size_t imageIdx = 0; imageIdx < numberOfXforms; ++imageIdx, cparam += parametersPerXform )
	    {
	    SplineWarpXform::SmartPtr xform = threadStorage->m_Xforms[imageIdx];
	    const UniformVolume* target = ThisConst->m_ImageVector[imageIdx];
	    const byte* dataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );
	  
	    for ( size_t dim = 0; dim < 3; ++dim, ++currentParameter )
	      {
	      const size_t cdparam = cparam + dim;
	      const size_t xfparam = 3 * cpIndex + dim;
	      const Types::Coordinate pStep = ThisConst->m_ParamStepArray[cdparam] * threadParameters->m_Step;
	    
	      if ( pStep > 0 )
		{
		const Types::Coordinate v0 = xform->GetParameter( xfparam );
		for ( int delta = 0; delta < 2; ++delta )
		  {
		  Types::Coordinate vTest = v0 + (2*delta-1) * pStep;
		  xform->SetParameter( xfparam, vTest );
		
		  xform->GetTransformedGridSequence( &(threadStorage->m_VectorList[0]), pixelsPerLineVOI, voi->From()[0], y, z );
		
		  size_t ofs = rowofs;
		  for ( size_t x = 0; x < pixelsPerLineVOI; ++x, ++ofs )
		    {
		    if ( threadStorage->m_Count[x] == numberOfImages ) // full count?
		      {
		      const byte baselineData = ThisConst->m_Data[imageIdx][ofs];
		      HistogramType& workingHistogram = threadStorage->m_Histogram[x];
		    
		      const size_t kernelIdx = ThisConst->m_StandardDeviationByPixel[ofs];
		      const size_t kernelRadius = ThisConst->m_HistogramKernelRadius[kernelIdx];
		      const HistogramBinType* kernel = ThisConst->m_HistogramKernel[kernelIdx];
		    
		      byte value;
		      if ( (baselineData != paddingValue) && target->ProbeData( value, dataPtr, threadStorage->m_VectorList[x] ) )
			{
			if ( value != paddingValue )
			  {
			  if ( delta )
			    ++threadStorage->m_CountByParameterPlus[currentParameter];
			  else
			    ++threadStorage->m_CountByParameterMinus[currentParameter];
			
			  if ( value != baselineData )
			    {
			    float fd = 0;
			    float invsamples = 1.0 / workingHistogram.SampleCount();
			  
			    const size_t idxMin = std::min( value, baselineData ) - kernelRadius;
			    const size_t idxMax = std::max( value, baselineData ) + kernelRadius + 1;
			  
			    for ( size_t idx = idxMin; idx < idxMax; ++idx )
			      {
			      const size_t idxV = abs(idx - value);
			      const float kernelIn = ( idxV < kernelRadius ) ? kernel[idxV] : 0;
			    
			      const size_t idxB = abs(idx - baselineData);
			      const float kernelOut = ( idxB < kernelRadius ) ? kernel[idxB] : 0;
			    
			      if ( (kernelIn > 0) || (kernelOut > 0) )
				{
				fd += MathUtil::plogp( invsamples * (workingHistogram[idx] - kernelOut + kernelIn) ) - MathUtil::plogp( invsamples * workingHistogram[idx] );
				}
			      }
			  
			    if ( delta ) 
			      threadStorage->m_FPlus[currentParameter] += fd;
			    else
			      threadStorage->m_FMinus[currentParameter] += fd;
			    }
			  }
			}
		      }
		    }
		  }
		xform->SetParameter( xfparam, v0 );
		}
	      }
	    }
	  }
	}
   
      // approximate gradient from upper and lower function evaluations  
      size_t currentParameterIdxLinear = 3 * numberOfXforms * (cpIndex - threadParameters->m_FirstIndexToCompute );
      size_t imageDimIdx = 0;
      for ( size_t cparam = 3*cpIndex; cparam < paramVectorDim; cparam += parametersPerXform )
	{
	for ( size_t dim = 0; dim < 3; ++dim, ++imageDimIdx, ++currentParameterIdxLinear )
	  {
	  threadParameters->m_Gradient[currentParameterIdxLinear] = 0.0;
	  assert( imageDimIdx < 3 * numberOfXforms );
	  if ( threadStorage->m_CountByParameterPlus[imageDimIdx] && threadStorage->m_CountByParameterMinus[imageDimIdx] )
	    {
	    double upper = threadStorage->m_FPlus[imageDimIdx] / threadStorage->m_CountByParameterPlus[imageDimIdx];
	    double lower = threadStorage->m_FMinus[imageDimIdx] / threadStorage->m_CountByParameterMinus[imageDimIdx];
	  
	    if ( (ThisConst->m_JacobianConstraintWeight > 0) || (ThisConst->m_BendingEnergyWeight > 0) )
	      {
	      const size_t cdparam = cparam+dim;
	      const size_t warpParam = cdparam % parametersPerXform;
	      const size_t imageIdx = cdparam / parametersPerXform;
	      const SplineWarpXform* xform = ThisConst->GetXformByIndex(imageIdx);
	      const Types::Coordinate pStep = ThisConst->m_ParamStepArray[cdparam] * threadParameters->m_Step;
	    
	      if ( ThisConst->m_JacobianConstraintWeight > 0 )
		{
		double upperJ, lowerJ;
		xform->GetJacobianConstraintDerivative( upperJ, lowerJ, warpParam, voi[warpParam], pStep );
	      
		upper -= ThisConst->m_JacobianConstraintWeight * upperJ;
		lower -= ThisConst->m_JacobianConstraintWeight * lowerJ;
		}
	    
	      if ( ThisConst->m_BendingEnergyWeight > 0 )
		{
		double upperE, lowerE;
		xform->GetGridEnergyDerivative( upperE, lowerE, warpParam, pStep );	      
		upper -= ThisConst->m_BendingEnergyWeight * upperE;
		lower -= ThisConst->m_BendingEnergyWeight * lowerE;
		}
	      }

	    if ( (upper > 0) || (lower > 0))
	      {
	      threadParameters->m_Gradient[currentParameterIdxLinear] = upper - lower;
	      }
	    }
	  }
	}

      This->m_ControlPointIndexLock.Lock();
      }
    This->m_ControlPointIndexLock.Unlock();

    This->m_ThreadReadySemaphore.Post();
    }
  
  return CMTK_THREAD_RETURN_VALUE;
}

} // namespace cmtk
