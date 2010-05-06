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

#include <mpi.h>

#ifdef DEBUG
#  define DEBUG_COMM
#endif
//#  define DEBUG_COMM

//#define DEBUG_GRADIENT

namespace
cmtk
{

/** \addtogroup Registration */
//@{

SplineWarpGroupwiseRegistrationRMIFunctional::ReturnType
SplineWarpGroupwiseRegistrationRMIFunctional::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
  const int tagCompute = 1;
  const int tagResults = 2;
  const int tagFinished = 3;

#ifdef DEBUG_GRADIENT
  g.SetAll( 12345 );
#endif

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  const size_t numberOfXforms = this->m_XformVector.size();
  
  const Self::ReturnType baseValue = this->EvaluateAt( v );

  if ( this->m_NeedsUpdateInformationByControlPoint )
    {
    this->UpdateInformationByControlPoint();
    }
  // allocate sufficiently many local thread data
  const size_t safeNumberOfThreads = 
    std::min( numberOfThreads, this->m_ControlPointScheduleOverlapFreeMaxLength );

  const size_t numberOfDOFsPerCell = 3 * numberOfXforms * safeNumberOfThreads;
  if ( this->m_ThreadSumOfProductsMatrix.size() < (2 * numberOfDOFsPerCell) )
    {
    this->m_ThreadSumOfProductsMatrix.resize( 2 * numberOfDOFsPerCell );
    }
  if ( this->m_ThreadSumsVector.size() < ( 2 * numberOfDOFsPerCell) )
    {
    this->m_ThreadSumsVector.resize( 2 * numberOfDOFsPerCell );
    }

  ThreadParameterArray<Self,EvaluateLocalGradientThreadParameters> threadParams( this, safeNumberOfThreads );

  const size_t numberOfControlPoints = this->m_ControlPointSchedule.size();
  // use more chunks than nodes to get more fine-grained parallelism
  const size_t controlPointsPerNode = 1 + static_cast<size_t>( numberOfControlPoints / (8 * this->m_SizeMPI) );
  // assign less work per chunk to root node to leave time to check for incoming results and
  // keep worker nodes maximally busy
  const size_t controlPointsRootNode = 1 + controlPointsPerNode / 8;

  MPI::Status status;

  size_t resultsToReceive = 0;
  std::vector<size_t> taskByNode( this->m_SizeMPI, -1 );

  const size_t msgBufferSize = 3 * sizeof( Types::Coordinate ) * controlPointsPerNode * numberOfXforms;
  std::vector<char> msgBuffer( msgBufferSize );
  char* msgBufferPtr = &msgBuffer[0];

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
#ifdef DEBUG_COMM
	    std::cerr << this->m_RankMPI << "\t" << "Sending task " << firstIndexToCompute << " to node " << dest << "\n";
#endif
	    MPI::COMM_WORLD.Ssend( &firstIndexToCompute, sizeof( firstIndexToCompute ), MPI::CHAR, dest, tagCompute );
	    ++resultsToReceive;
#ifdef DEBUG_COMM
	    std::cerr << this->m_RankMPI << "\t" << "Sent task " << firstIndexToCompute << " to node " << dest << "\n";
#endif

	    firstIndexToCompute += controlPointsPerNode;
	    }
	  }
	}
      }
    else
      {
      // receive job instructions      
#ifdef DEBUG_COMM
      std::cerr << this->m_RankMPI << "\t" << "Receiving on node " << this->m_RankMPI << "\n";
#endif
      MPI::COMM_WORLD.Recv( msgBufferPtr, msgBufferSize, MPI::CHAR, 0 /*root*/ , MPI::ANY_TAG, status );
      if ( status.Get_tag() == tagFinished )
	{
#ifdef DEBUG_COMM
	std::cerr << "Exiting on node " << this->m_RankMPI << "\n";
#endif
	allJobsSent = true;
	}
      else
	{
	int position = 0;
#ifdef DEBUG_COMM
	std::cerr << this->m_RankMPI << "\t" << "Unpacking task on node " << this->m_RankMPI << "\n";
#endif
	MPI::CHAR.Unpack( msgBufferPtr, msgBufferSize, &firstIndexToCompute, sizeof( firstIndexToCompute ), position, MPI::COMM_WORLD );
	}
      }
    
    if ( (firstIndexToCompute < numberOfControlPoints ) && ! allJobsSent )
      {
      for ( size_t thread = 0; thread < safeNumberOfThreads; ++thread )
	{
	threadParams[thread].m_ThreadStorageIndex = thread;
	threadParams[thread].m_Step = step;
	threadParams[thread].m_Gradient = (Types::Coordinate*)msgBufferPtr;
	threadParams[thread].m_MetricBaseValue = baseValue;
	threadParams[thread].m_FirstIndexToCompute = firstIndexToCompute;
	}

      const size_t cpsToCompute = 
	std::min( this->m_RankMPI ? controlPointsPerNode : controlPointsRootNode, 
		  (int)numberOfControlPoints - firstIndexToCompute );
#ifdef DEBUG_COMM
      std::cerr << this->m_RankMPI << "\t" << "Computing " << firstIndexToCompute << ":" << cpsToCompute << "\n";
#endif
      threadParams.RunInParallelFIFO
	( EvaluateLocalGradientThreadFunc, cpsToCompute, firstIndexToCompute );

      // on root node, simply sort results into gradient fields
      if ( this->m_RankMPI == 0 )
	{
	// determine last control point processed
	const size_t maxCpIndex = firstIndexToCompute + cpsToCompute;
#ifdef DEBUG_COMM
	StdErr.printf( "%d\tCopying local result %d:%d\n", this->m_RankMPI, firstIndexToCompute, maxCpIndex );
#endif
	// sort computed values into gradient vector
	this->ReorderGradientComponents( g.Elements, (Types::Coordinate*)msgBufferPtr, firstIndexToCompute, maxCpIndex );
	firstIndexToCompute += cpsToCompute;
	}
      else
	{
	// send results
	MPI::COMM_WORLD.Ssend( msgBufferPtr, 3 * numberOfXforms * cpsToCompute * sizeof( Types::Coordinate ), 
			      MPI::CHAR, 0, tagResults );
#ifdef DEBUG_COMM
	std::cerr << this->m_RankMPI << "\t" << "Sent results of size "
		  << 3 * numberOfXforms * cpsToCompute * sizeof( Types::Coordinate ) << "\n";
#endif
	}
      }
    
    // are we still waiting for results from worker nodes?
    // are the worker nodes still trying to send data?
    if ( this->m_RankMPI == 0 )
      {
      // receive results:
      // is a result ready to be received?
      while ( resultsToReceive && MPI::COMM_WORLD.Iprobe( MPI::ANY_SOURCE, tagResults, status ) )
	{
	// receive result
	MPI::COMM_WORLD.Recv( msgBufferPtr, msgBufferSize, MPI::CHAR, status.Get_source(), tagResults, status );
		
	// what control point did we send this node for processing?
	const int fromNode = status.Get_source();
	const int receivedJob = taskByNode[fromNode];
	
	if ( receivedJob != -1 )
	  {
	  --resultsToReceive;
	  taskByNode[fromNode] = -1;
	  
	  // determine last control point processed
	  const size_t maxCpIndex = std::min<size_t>( receivedJob + controlPointsPerNode, numberOfControlPoints );
	  this->ReorderGradientComponents( g.Elements, (Types::Coordinate*)msgBufferPtr, receivedJob, maxCpIndex );
#ifdef DEBUG_COMM
	  StdErr.printf( "%d\tReceived result %d:%d of size %d from node %d\n", this->m_RankMPI, receivedJob, maxCpIndex, status.Get_elements( MPI::CHAR ), fromNode );
#endif
	  }
	else
	  {
	  StdErr.printf( "WARNNING: Received extraneous result from node %d\n", fromNode );
	  }
	}
      }
    }

  // collect remaining results from worker nodes
  while ( resultsToReceive )
    {
    // receive result
    MPI::COMM_WORLD.Recv( msgBufferPtr, msgBufferSize, MPI::CHAR, MPI::ANY_SOURCE, MPI::ANY_TAG, status );
    
#ifdef DEBUG_COMM
    StdErr.printf( "%d\tFinalizing - received result from node %d\n", this->m_RankMPI, status.Get_source() );
#endif
    
    // what control point did we send this node for processing?
    const int fromNode = status.Get_source();
    const int receivedJob = taskByNode[fromNode];

    if ( receivedJob != -1 )
      {
      --resultsToReceive;
      taskByNode[fromNode] = -1;
    
      // determine last control point processed
      const size_t maxCpIndex = std::min<size_t>( receivedJob + controlPointsPerNode, numberOfControlPoints );
      this->ReorderGradientComponents( g.Elements, (Types::Coordinate*)msgBufferPtr, receivedJob, maxCpIndex );
      }
    else
      {
      StdErr.printf( "WARNING: Received extraneous result from node %d\n", fromNode );
      }
     }  
 
  if ( this->m_RankMPI == 0 )
    { 
    // tell all worker nodes to finish their computation loops
#ifdef DEBUG_COMM
    StdErr << this->m_RankMPI << "\t" << "Leaving computation loop\n";
#endif
    for ( int workerNode = 1; workerNode < this->m_SizeMPI; ++workerNode )
      {
      MPI::COMM_WORLD.Ssend( msgBufferPtr, 1, MPI::CHAR, workerNode, tagFinished );
      }
    }
  
#ifdef DEBUG_COMM
  StdErr << this->m_RankMPI << "\t" << "Broadcasting gradient estimate\n";
#endif
  
  // Distribute gradient estimate
  MPI::COMM_WORLD.Bcast( g.Elements, g.Dim * sizeof( g.Elements[0] ), MPI::CHAR, 0 /*root*/ );

#ifdef DEBUG_COMM
  StdErr << this->m_RankMPI << "\t" << "Done with gradient\n\n";
#endif
  
//  StdErr.printf( "gmax = %f\n" , g.MaxNorm() );
  
  if ( this->m_PartialGradientMode )
    {
    const Types::Coordinate gthresh = g.MaxNorm() * this->m_PartialGradientThreshold;
    for ( size_t param = 0; param < g.Dim; ++param )
      {
#ifdef DEBUG_GRADIENT
      if ( this->m_RankMPI == 0 )
	if ( g[param] == 12345 )
	  {
	  std::cerr << param << "\t";
	  }
#endif
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

#ifdef DEBUG_COMM
  MPI::COMM_WORLD.Barrier();
#endif
  return baseValue;
}

void
SplineWarpGroupwiseRegistrationRMIFunctional
::ReorderGradientComponents
( Types::Coordinate *const dst, const Types::Coordinate* src,
  const size_t fromCpIdx, const size_t toCpIdx )
{
  const size_t numberOfParameters = this->ParamVectorDim();

  // sort computed values into gradient vector
  size_t currentParameter = 0;
  for ( size_t cpIndex = fromCpIdx; cpIndex < toCpIdx; ++cpIndex )
    {
    const size_t cp = this->m_ControlPointSchedule[cpIndex];
    
    for ( size_t cparam = 3*cp; cparam < numberOfParameters; cparam += this->m_ParametersPerXform )
      {
      for ( size_t dim = 0; dim < 3; ++dim, ++currentParameter )
	{
	dst[cparam+dim] = src[currentParameter];
	}
      }
    }
}

CMTK_THREAD_RETURN_TYPE
SplineWarpGroupwiseRegistrationRMIFunctional::EvaluateLocalGradientThreadFunc
( void* args )
{
  EvaluateLocalGradientThreadParameters* threadParameters = static_cast<EvaluateLocalGradientThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  const size_t threadID = threadParameters->ThisThreadIndex;
  const size_t threadStorageIndex = threadParameters->m_ThreadStorageIndex;

#ifndef DEBUG_COMM  
  if ( !(threadID % 100) )
    {
    StdErr.printf( "%d / %d\r", (int)threadID, (int)ThisConst->m_ControlPointSchedule.size() );
    }
#endif
  
  const size_t cpIndex = ThisConst->m_ControlPointSchedule[threadID];
  const DataGrid::RegionType& voi = ThisConst->m_VolumeOfInfluenceArray[cpIndex];
  const size_t pixelsPerLineVOI = (voi.To()[0]-voi.From()[0]);
  std::vector<Vector3D> vectorList( pixelsPerLineVOI );
  std::vector<size_t> count( pixelsPerLineVOI );
  
  const size_t numberOfXforms = ThisConst->m_XformVector.size();
  std::vector<size_t> totalNumberOfSamples( 6 * numberOfXforms );
  std::fill( totalNumberOfSamples.begin(), totalNumberOfSamples.end(), ThisConst->m_TotalNumberOfSamples );

  const size_t parametersPerXform = ThisConst->m_ParametersPerXform;
  const size_t paramVectorDim = ThisConst->ParamVectorDim();

  const byte paddingValue = ThisConst->m_PaddingValue;
  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const size_t numberOfImages = imagesTo - imagesFrom;

  const UniformVolume* templateGrid = ThisConst->m_TemplateGrid;

  const size_t threadDataIdx = 6 * threadStorageIndex * numberOfXforms;
  for ( size_t image = 0; image < 6 * numberOfXforms; ++image )
    {
    const SumsAndProductsVectorType& srcSumOfProducts = ThisConst->m_SumOfProductsMatrix;
    SumsAndProductsVectorType& dstSumOfProducts = This->m_ThreadSumOfProductsMatrix[threadDataIdx + image];
    dstSumOfProducts.resize( srcSumOfProducts.size() );
    std::copy( srcSumOfProducts.begin(), srcSumOfProducts.end(), dstSumOfProducts.begin() );

    const SumsAndProductsVectorType& srcSumsVector = ThisConst->m_SumsVector;
    SumsAndProductsVectorType& dstSumsVector = This->m_ThreadSumsVector[threadDataIdx + image];
    dstSumsVector.resize( srcSumsVector.size() );
    std::copy( srcSumsVector.begin(), srcSumsVector.end(), dstSumsVector.begin() );
    }
  
  for ( int z = voi.From()[2]; z < voi.To()[2]; ++z ) 
    {
    for ( int y = voi.From()[1]; y < voi.To()[1]; ++y )
      {      
      // check which pixels in this row have a full sample count
      const size_t rowofs = templateGrid->GetOffsetFromIndex( voi.From()[0], y, z );

      std::fill( count.begin(), count.end(), 0 );
      for ( size_t img = 0; img < numberOfXforms; ++img )
	{ 
	const byte* dataPtr = ThisConst->m_Data[img]+rowofs;
	for ( size_t x = 0; x < pixelsPerLineVOI; ++x )
	  {
	  const byte dataThisPixel = dataPtr[x];
	  if ( dataThisPixel != paddingValue )
	    {
	    ++count[x];
	    }
	  }
	}
      
      size_t cparam = 3 * cpIndex;
      size_t currentParameter = 0;
      for ( size_t img = 0; img < numberOfXforms; ++img, cparam += parametersPerXform )
	{
	SplineWarpXform::SmartPtr xform = This->GetXformByIndex( img );
	const UniformVolume* target = ThisConst->m_ImageVector[img];
	const byte* targetDataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );
	
	for ( size_t dim = 0; dim < 3; ++dim )
	  {
	  const size_t cdparam = cparam + dim;
	  const size_t xfparam = 3 * cpIndex + dim;
	  const Types::Coordinate pStep = ThisConst->m_ParamStepArray[cdparam] * threadParameters->m_Step;

	  if ( pStep > 0 )
	    {
	    const Types::Coordinate v0 = xform->GetParameter( xfparam );
	    for ( int delta = 0; delta < 2; ++delta, ++currentParameter )
	      {
	      SumsAndProductsVectorType& dstSumOfProducts = This->m_ThreadSumOfProductsMatrix[threadDataIdx+currentParameter];
	      SumsAndProductsVectorType& dstSumsVector = This->m_ThreadSumsVector[threadDataIdx+currentParameter];
	      
	      Types::Coordinate vTest = v0 + (2*delta-1) * pStep;
	      xform->SetParameter( xfparam, vTest );
	      xform->GetTransformedGridSequence( &(vectorList[0]), pixelsPerLineVOI, voi.From()[0], y, z );
	      
	      byte* rowDataPtr = ThisConst->m_Data[img] + rowofs;
	      for ( size_t x = 0; x < pixelsPerLineVOI; ++x, ++rowDataPtr )
		{
		const int baselineData = *rowDataPtr;
		if ( (count[x] == numberOfImages) || ((count[x] == numberOfImages-1) && (baselineData == paddingValue) ) ) // full count?
		  {
		  byte newData;
		  if ( !target->ProbeData( newData, targetDataPtr, vectorList[x] ) )
		    newData = paddingValue;
		  
		  if ( newData != baselineData )
		    {
		    if ( baselineData != paddingValue )
		      {
		      dstSumsVector[img] -= baselineData;
		      size_t midx = 0;
		      for ( size_t img2 = 0; img2 < numberOfImages; ++img2 )
			{
			for ( size_t otherImg = 0; otherImg <= img2; ++otherImg, ++midx )
			  {
			  if ( img2 == img ) 
			    {
			    const int otherData = ThisConst->m_Data[otherImg][rowofs+x];
			    dstSumOfProducts[midx] -= baselineData * otherData;
			    }
			  else
			    {
			    if ( otherImg == img )
			      {
			      const int otherData = ThisConst->m_Data[img2][rowofs+x];
			      dstSumOfProducts[midx] -= baselineData * otherData;
			      }
			    }
			  }
			}
		      }
		    
		    if ( newData != paddingValue )
		      {
		      if ( count[x] == numberOfImages-1 )
			{
			++totalNumberOfSamples[currentParameter];
			}

		      dstSumsVector[img] += newData;
		      size_t midx = 0;
		      for ( size_t img2 = imagesFrom; img2 < imagesTo; ++img2 )
			{
			for ( size_t otherImg = imagesFrom; otherImg <= img2; ++otherImg, ++midx )
			  {
			  if ( img2 == img )
			    {
			    if ( otherImg == img )
			      {
			      dstSumOfProducts[midx] += newData * newData;
			      }
			    else
			      {
			      const int otherData = ThisConst->m_Data[otherImg][rowofs+x];
			      dstSumOfProducts[midx] += newData * otherData;
			      }
			    }
			  else
			    {
			    if ( otherImg == img )
			      {
			      const int otherData = ThisConst->m_Data[img2][rowofs+x];
			      dstSumOfProducts[midx] += newData * otherData;
			      }
			    }
			  }
			}
		      }
		    else
		      {
		      if ( count[x] == numberOfImages )
			{
			--totalNumberOfSamples[currentParameter];
			}
		      }
		    }
		  }
		}
	      }
	    xform->SetParameter( xfparam, v0 );
	    }
	  else
	    {
	    currentParameter += 2;
	    }
	  }
	}
      }
    }

  Matrix2D<Self::ReturnType> covarianceMatrix( numberOfImages, numberOfImages );
  
  // approximate gradient from upper and lower function evaluations  
  size_t img = 0, currentParameter = 0;
  const size_t parameterOffset = 3 * numberOfXforms * (threadID - threadParameters->m_FirstIndexToCompute);

  const float fBaseValue = threadParameters->m_MetricBaseValue;
  for ( size_t cparam = 3*cpIndex; cparam < paramVectorDim; cparam += parametersPerXform )
    {
    for ( size_t dim = 0; dim < 3; ++dim, ++img, currentParameter += 2 )
      {
      const Self::ReturnType fMinus = ThisConst->GetMetric( This->m_ThreadSumOfProductsMatrix[threadDataIdx+currentParameter], This->m_ThreadSumsVector[threadDataIdx+currentParameter], totalNumberOfSamples[currentParameter], covarianceMatrix );
      const Self::ReturnType fPlus = ThisConst->GetMetric( This->m_ThreadSumOfProductsMatrix[threadDataIdx+currentParameter+1], This->m_ThreadSumsVector[threadDataIdx+currentParameter+1], totalNumberOfSamples[currentParameter], covarianceMatrix );

      if ( (fPlus > fBaseValue) || (fMinus > fBaseValue) )
	{
	threadParameters->m_Gradient[parameterOffset + currentParameter / 2] = fPlus - fMinus;
	}
      else
	{
	threadParameters->m_Gradient[parameterOffset + currentParameter / 2] = 0.0;
	}
      }
    }

  return CMTK_THREAD_RETURN_VALUE;
}

} // namespace cmtk
