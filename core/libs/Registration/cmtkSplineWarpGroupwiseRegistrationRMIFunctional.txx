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

namespace
cmtk
{

/** \addtogroup Registration */
//@{

SplineWarpGroupwiseRegistrationRMIFunctional::ReturnType
SplineWarpGroupwiseRegistrationRMIFunctional::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
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

  if ( this->m_ThreadSumOfProductsMatrix.size() < (6 * numberOfXforms * safeNumberOfThreads) )
    {
    this->m_ThreadSumOfProductsMatrix.resize( 6 * numberOfXforms * safeNumberOfThreads );
    }
  if ( this->m_ThreadSumsVector.size() < (6 * numberOfXforms * safeNumberOfThreads) )
    {
    this->m_ThreadSumsVector.resize( 6 * numberOfXforms * safeNumberOfThreads );
    }

  ThreadParameterArray<Self,EvaluateLocalGradientThreadParameters> threadParams( this, safeNumberOfThreads );
  for ( size_t thread = 0; thread < safeNumberOfThreads; ++thread )
    {
    threadParams[thread].m_ThreadStorageIndex = thread;
    threadParams[thread].m_Step = step;
    threadParams[thread].m_Gradient = g.Elements;
    threadParams[thread].m_MetricBaseValue = baseValue;
    }

  threadParams.RunInParallelFIFO( EvaluateLocalGradientThreadFunc, this->m_ControlPointSchedule.size() );

  if ( this->m_PartialGradientMode )
    {
    const Types::Coordinate gthresh = g.MaxNorm() * this->m_PartialGradientThreshold;
    for ( size_t param = 0; param < g.Dim; ++param )
      if ( fabs( g[param] ) < gthresh )
	{
	g[param] = this->m_ParamStepArray[param] = 0.0;
	}
    }

  if ( this->m_ForceZeroSum )
    {
    this->ForceZeroSumGradient( g );
    }
  
  return baseValue;
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
  
  if ( !(threadID % 100) )
    {
    std::cerr << threadID << " / " << ThisConst->m_ControlPointSchedule.size() << "\r";
    }
  
  const size_t cpIndex = ThisConst->m_ControlPointSchedule[threadID];
  std::vector<DataGrid::RegionType>::const_iterator voi = ThisConst->m_VolumeOfInfluenceArray.begin() + cpIndex;

  const size_t pixelsPerLineVOI = (voi->To()[0]-voi->From()[0]);
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
  
  for ( int z = voi->From()[2]; (z < voi->To()[2]); ++z ) 
    {
    for ( int y = voi->From()[1]; (y < voi->To()[1]); ++y )
      {      
      // check which pixels in this row have a full sample count
      const size_t rowofs = templateGrid->GetOffsetFromIndex( voi->From()[0], y, z );

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
	SplineWarpXform::SmartPtr xform = This->GetXformByIndex(img);
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
	      xform->GetTransformedGridRow( pixelsPerLineVOI, &(vectorList[0]), voi->From()[0], y, z );
	      
	      byte* rowDataPtr = ThisConst->m_Data[img] + rowofs;
	      for ( size_t x = 0; x < pixelsPerLineVOI; ++x, ++rowDataPtr )
		{
		const int baselineData = *rowDataPtr;
		if ( (count[x] == numberOfImages) || 
		     ((count[x] == numberOfImages-1) && (baselineData == paddingValue) ) ) // full count?
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
		      for ( size_t img2 = imagesFrom; img2 < imagesTo; ++img2 )
			{
			for ( size_t otherImg = imagesFrom; otherImg <= img2; ++otherImg, ++midx )
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
  const Functional::ReturnType fBaseValue = threadParameters->m_MetricBaseValue;
  for ( size_t cparam = 3*cpIndex; cparam < paramVectorDim; cparam += parametersPerXform )
    {
    for ( size_t dim = 0; dim < 3; ++dim, ++img, currentParameter += 2 )
      {
      const Self::ReturnType fMinus = ThisConst->GetMetric( This->m_ThreadSumOfProductsMatrix[threadDataIdx+currentParameter], This->m_ThreadSumsVector[threadDataIdx+currentParameter], totalNumberOfSamples[currentParameter], covarianceMatrix );
      const Self::ReturnType fPlus = ThisConst->GetMetric( This->m_ThreadSumOfProductsMatrix[threadDataIdx+currentParameter+1], This->m_ThreadSumsVector[threadDataIdx+currentParameter+1], totalNumberOfSamples[currentParameter], covarianceMatrix );

      if ( (fPlus > fBaseValue) || (fMinus > fBaseValue) )
	{
	threadParameters->m_Gradient[cparam+dim] = fPlus - fMinus;
	}
      else
	{
	threadParameters->m_Gradient[cparam+dim] = 0.0;
	}
      }
    }
  
  return CMTK_THREAD_RETURN_VALUE;
}

} // namespace cmtk
