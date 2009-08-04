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
  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  this->m_ThreadHistograms.resize( numberOfThreads );
  
  const Self::ReturnType baseValue = this->EvaluateAt( v );

  this->m_ControlPointIndexNext = 0;
  this->m_ControlPointIndexLast = this->m_ParametersPerXform / 3;
  
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
      
  ThreadParameterArray<Self,EvaluateLocalGradientThreadParameters> threadParams( this, numberOfThreads );
  for ( size_t thread = 0; thread < numberOfThreads; ++thread )
    {
    threadParams[thread].m_Step = step;
    threadParams[thread].m_Gradient = g.Elements;
    }
  threadParams.RunInParallel( EvaluateLocalGradientThreadFunc );

  if ( this->m_PartialGradientMode )
    {
    const Types::Coordinate gthresh = g.MaxNorm() * this->m_PartialGradientThreshold;
#pragma omp parallel for
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
  
  Self::StaticThreadStorage* threadStorage = (&This->m_StaticThreadStorage[threadID]);
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

  This->m_ControlPointIndexLock.Lock();
  while ( (cpIndex = This->m_ControlPointIndexNext) < This->m_ControlPointIndexLast )
    {
    ++This->m_ControlPointIndexNext;
    This->m_ControlPointIndexLock.Unlock();

    if ( !threadID && !(cpIndex % 1000) )
      {
      std::cerr << cpIndex << " / " << numberOfControlPoints << "\r";
      }    
    
    const Rect3D* voi = &(ThisConst->m_VolumeOfInfluenceArray[cpIndex]);
    const size_t pixelsPerLineVOI = (voi->endX-voi->startX);

    std::fill( threadStorage->m_FPlus.begin(), threadStorage->m_FPlus.end(), static_cast<Self::ReturnType>( 0 ) );
    std::fill( threadStorage->m_FMinus.begin(), threadStorage->m_FMinus.end(), static_cast<Self::ReturnType>( 0 ) );
    std::fill( threadStorage->m_CountByParameterPlus.begin(), threadStorage->m_CountByParameterPlus.end(), 0 );
    std::fill( threadStorage->m_CountByParameterMinus.begin(), threadStorage->m_CountByParameterMinus.end(), 0 );
    
    for ( int z = voi->startZ; (z < voi->endZ); ++z ) 
      {
      for ( int y = voi->startY; (y < voi->endY); ++y )
	{
	// evaluate baseline entropies for this row
	const size_t rowofs = templateGrid->GetOffsetFromIndex( voi->startX, y, z );
	
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
		
		xform->GetTransformedGridSequence( &(threadStorage->m_VectorList[0]), pixelsPerLineVOI, voi->startX, y, z );
		
		size_t ofs = rowofs;
		for ( size_t x = 0; x < pixelsPerLineVOI; ++x, ++ofs )
		  {
		  if ( threadStorage->m_Count[x] == numberOfImagesIncludingTemplate ) // full count?
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
			  double fd = 0;
			  double invsamples = 1.0 / workingHistogram.SampleCount();
			  
			  const size_t idxMin = std::min( value, baselineData ) - kernelRadius;
			  const size_t idxMax = std::max( value, baselineData ) + kernelRadius + 1;
			  
			  for ( size_t idx = idxMin; idx < idxMax; ++idx )
			    {
			    const size_t idxV = abs((int)(idx - value));
			    const double kernelIn = ( idxV < kernelRadius ) ? kernel[idxV] : 0;
			    
			    const size_t idxB = abs((int)(idx - baselineData));
			    const double kernelOut = ( idxB < kernelRadius ) ? kernel[idxB] : 0;
			    
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
    
    size_t imageDimIdx = 0;
    for ( size_t cparam = 3*cpIndex; cparam < paramVectorDim; cparam += parametersPerXform )
      {
      for ( size_t dim = 0; dim < 3; ++dim, ++imageDimIdx )
	{
	const size_t cdparam = cparam+dim;
	threadParameters->m_Gradient[cdparam] = 0.0;
	if ( threadStorage->m_CountByParameterPlus[imageDimIdx] && threadStorage->m_CountByParameterMinus[imageDimIdx] )
	  {
	  double upper = threadStorage->m_FPlus[imageDimIdx] / threadStorage->m_CountByParameterPlus[imageDimIdx];
	  double lower = threadStorage->m_FMinus[imageDimIdx] / threadStorage->m_CountByParameterMinus[imageDimIdx];
	  
	  if ( (ThisConst->m_JacobianConstraintWeight > 0) || (ThisConst->m_BendingEnergyWeight > 0) )
	    {
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
	    threadParameters->m_Gradient[cparam+dim] = upper - lower;
	    }
	  }
	}
      }
    This->m_ControlPointIndexLock.Lock();
    }
  This->m_ControlPointIndexLock.Unlock();

  return CMTK_THREAD_RETURN_VALUE;
}

} // namespace cmtk
