/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include <System/cmtkDebugOutput.h>

#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkTypes.h>

#include <math.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TMetricFunctional>
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::SplineWarpMultiChannelRegistrationFunctional() 
  : m_AdaptiveFixEntropyThreshold( false ),
    m_AdaptiveFixThreshFactor( 0.5 ), 
    m_JacobianConstraintWeight( 0.0 ),
    m_NumberOfThreads( ThreadPool::GetGlobalThreadPool().GetNumberOfThreads() )
{
}

template<class TMetricFunctional>
template<class TAffineMetricFunctional>
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::SplineWarpMultiChannelRegistrationFunctional
( AffineMultiChannelRegistrationFunctional<TAffineMetricFunctional>& affineFunctional )
  : m_AdaptiveFixEntropyThreshold( false ),
    m_AdaptiveFixThreshFactor( 0.5 ),
    m_JacobianConstraintWeight( 0.0 ),
    m_NumberOfThreads( ThreadPool::GetGlobalThreadPool().GetNumberOfThreads() )
{
  this->SetInitialAffineTransformation( affineFunctional.GetTransformation() );
  this->AddReferenceChannels( affineFunctional.m_ReferenceChannels.begin(), affineFunctional.m_ReferenceChannels.end() );
  this->AddFloatingChannels( affineFunctional.m_FloatingChannels.begin(), affineFunctional.m_FloatingChannels.end() );
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::InitTransformation( const Vector3D& domain, const Types::Coordinate gridSpacing, const bool exact )
{
  if ( this->m_ReferenceChannels.size() == 0 )
    {
    StdErr << "ERROR: call to SplineWarpMultiChannelRegistrationFunctional::InitTransformation() before reference channel image was set.\n";
    exit( 1 );
    }

  this->m_Transformation.Init( domain, gridSpacing, &this->m_InitialAffineTransformation, exact );
  this->m_ThreadTransformations.resize( this->m_NumberOfThreads, SplineWarpXform::SmartPtr::Null() );
  for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
    {
    this->m_ThreadTransformations[thread] = SplineWarpXform::SmartPtr( new SplineWarpXform( domain, gridSpacing, &this->m_InitialAffineTransformation, exact ) );
    }
  this->UpdateTransformationData();
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::RefineTransformation()
{
  this->m_Transformation.Refine();
  for ( size_t thread = 0; thread < this->m_ThreadTransformations.size(); ++thread )
    {
    this->m_ThreadTransformations[thread]->Refine();
    }
  this->UpdateTransformationData();
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::UpdateTransformationData()
{
  this->m_Transformation.RegisterVolume( this->m_ReferenceChannels[0] );
  for ( size_t thread = 0; thread < this->m_ThreadTransformations.size(); ++thread )
    {
    this->m_ThreadTransformations[thread]->RegisterVolume( this->m_ReferenceChannels[0] );
    }
  
  this->m_StepScaleVector.resize( this->m_Transformation.VariableParamVectorDim() );
  for ( size_t idx = 0; idx < this->m_StepScaleVector.size(); ++idx ) 
    {
    this->m_StepScaleVector[idx] = this->GetParamStep( idx );
    }
  
  const size_t numberOfControlPoints = this->m_Transformation.VariableParamVectorDim() / 3;
  this->m_VolumeOfInfluenceVector.resize( numberOfControlPoints );
  
  const UniformVolume::CoordinateRegionType referenceDomain( this->m_ReferenceChannels[0]->m_Offset, this->m_ReferenceChannels[0]->Size );
  
  for ( size_t idx = 0; idx < numberOfControlPoints; ++idx ) 
    {
    this->m_VolumeOfInfluenceVector[idx] = this->m_ReferenceChannels[0]->GetGridRange( this->m_Transformation.GetVolumeOfInfluence( idx * 3, referenceDomain ) );
    }

  m_UpdateTransformationFixedControlPointsRequired = true;
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::UpdateTransformationFixedControlPoints()
{
  this->m_UpdateTransformationFixedControlPointsRequired = false;

  std::vector<Types::DataItemRange> valueRange( this->m_NumberOfChannels );

  size_t channel = 0;
  for ( size_t ref = 0; ref < this->m_ReferenceChannels.size(); ++ref, ++channel )
    {
    valueRange[channel] = this->m_ReferenceChannels[ref]->GetData()->GetRange();
    }
  for ( size_t flt = 0; flt < this->m_FloatingChannels.size(); ++flt, ++channel )
    {
    valueRange[channel] = this->m_FloatingChannels[flt]->GetData()->GetRange();
    }
  
  const size_t numberOfControlPoints = this->m_Transformation.VariableParamVectorDim() / 3;

  const UniformVolume::CoordinateRegionType referenceDomain( this->m_ReferenceChannels[0]->m_Offset, this->m_ReferenceChannels[0]->Size );

  std::vector<bool> active( numberOfControlPoints );
  std::fill( active.begin(), active.end(), true );

  if ( this->m_AdaptiveFixThreshFactor > 0 ) 
    {
    if ( this->m_AdaptiveFixEntropyThreshold )
      {
      Histogram<unsigned int> histogram( 128 );
      std::vector<float> channelEntropy( numberOfControlPoints * this->m_NumberOfChannels );
      std::vector<float> minEntropy( this->m_NumberOfChannels );
      std::vector<float> maxEntropy( this->m_NumberOfChannels );
      size_t channelEntropyIdx = 0;
      
      for ( size_t cp = 0; cp < numberOfControlPoints; ++cp ) 
	{
	// We cannot use the precomputed table of VOIs here because in "fast" mode, these VOIs are smaller than we want them here.
	const DataGrid::RegionType region = this->m_ReferenceChannels[0]->GetGridRange( this->m_Transformation.GetVolumeOfInfluence( 3 * cp, referenceDomain, false /*fast*/ ) );
	
	for ( size_t channel = 0; channel < this->m_NumberOfChannels; ++channel )
	  {
	  histogram.SetRange( valueRange[channel] );
	  
	  size_t r = region.From()[0] + this->m_ReferenceDims[0] * ( region.From()[1] + this->m_ReferenceDims[1] * region.From()[2] );
	  const int endOfLine = ( region.From()[0] + ( this->m_ReferenceDims[0]-region.To()[0]) );
	  const int endOfPlane = this->m_ReferenceDims[0] * ( region.From()[1] + (this->m_ReferenceDims[1]-region.To()[1]) );
	  
	  if ( channel < this->m_ReferenceChannels.size() )
	    {
	    const TypedArray* refChannel = this->m_ReferenceChannels[channel]->GetData();
	    for ( int pZ = region.From()[2]; pZ<region.To()[2]; ++pZ, r += endOfPlane ) 
	      for ( int pY = region.From()[1]; pY<region.To()[1]; ++pY, r += endOfLine ) 
		for ( int pX = region.From()[0]; pX<region.To()[0]; ++pX, ++r ) 
		  {
		  Types::DataItem refValue;
		  if ( refChannel->Get( refValue, r ) )
		    histogram.Increment( histogram.ValueToBin( refValue ) );
		  }
	    }
	  else
	    {
	    const float* fltChannel = &(this->m_ReformattedFloatingChannels[channel-this->m_ReferenceChannels.size()][0]);
	    for ( int pZ = region.From()[2]; pZ<region.To()[2]; ++pZ, r += endOfPlane ) 
	      for ( int pY = region.From()[1]; pY<region.To()[1]; ++pY, r += endOfLine ) 
		for ( int pX = region.From()[0]; pX<region.To()[0]; ++pX, ++r ) 
		  {
		  const float fltValue = fltChannel[r];
		  if ( finite( fltValue ) )
		    histogram.Increment( histogram.ValueToBin( fltValue ) );
		  }
	    }
	  channelEntropy[channelEntropyIdx++] = static_cast<float>( histogram.GetEntropy() );
	  }
	}
      
      size_t idx = 0;
      for ( size_t channel = 0; channel < this->m_NumberOfChannels; ++channel, ++idx ) 
	{
	minEntropy[channel] = channelEntropy[idx];
	maxEntropy[channel] = channelEntropy[idx];
	}
      for ( size_t cp = 1; cp < numberOfControlPoints; ++cp ) 
	{
	for ( size_t channel = 0; channel < this->m_NumberOfChannels; ++channel, ++idx ) 
	  {
	  minEntropy[channel] = std::min( minEntropy[channel], channelEntropy[idx] );
	  maxEntropy[channel] = std::max( maxEntropy[channel], channelEntropy[idx] );
	  }
	}
      
      for ( size_t channel = 0; channel < this->m_NumberOfChannels; ++channel, ++idx ) 
	{
	minEntropy[channel] += this->m_AdaptiveFixThreshFactor * (maxEntropy[channel] - minEntropy[channel]);
	}
      
      idx = 0;
      for ( size_t cp=0; cp<numberOfControlPoints; ++cp ) 
	{
	active[cp] = false;
	for ( size_t channel = 0; channel < this->m_NumberOfChannels; ++channel, ++idx ) 
	  {
	  if ( channelEntropy[idx] > minEntropy[channel] )
	    {
	    active[cp] = true;
	    }
	  }
	}
      }
    else // use intensity, not entropy threshold
      {
      // scale min values with threshold factor
      for ( size_t channel = 0; channel < this->m_NumberOfChannels; ++channel )
	{
	valueRange[channel].m_LowerBound = valueRange[channel].m_LowerBound + valueRange[channel].Width() * this->m_AdaptiveFixThreshFactor;
	}
      
      for ( size_t cp = 0; cp < numberOfControlPoints; ++cp ) 
	{
	// We cannot use the precomputed table of VOIs here because in "fast" mode, these VOIs are smaller than we want them here.
	const DataGrid::RegionType region = this->m_ReferenceChannels[0]->GetGridRange( this->m_Transformation.GetVolumeOfInfluence( 3 * cp, referenceDomain, false /*fast*/ ) );
	
	active[cp] = false;

	for ( size_t channel = 0; channel < this->m_NumberOfChannels; ++channel )
	  {
	  size_t r = region.From()[0] + this->m_ReferenceDims[0] * ( region.From()[1] + this->m_ReferenceDims[1] * region.From()[2] );
	  const int endOfLine = ( region.From()[0] + ( this->m_ReferenceDims[0]-region.To()[0]) );
	  const int endOfPlane = this->m_ReferenceDims[0] * ( region.From()[1] + (this->m_ReferenceDims[1]-region.To()[1]) );
	  
	  if ( channel < this->m_ReferenceChannels.size() )
	    {
	    const TypedArray* refChannel = this->m_ReferenceChannels[channel]->GetData();
	    for ( int pZ = region.From()[2]; pZ<region.To()[2]; ++pZ, r += endOfPlane ) 
	      for ( int pY = region.From()[1]; pY<region.To()[1]; ++pY, r += endOfLine ) 
		for ( int pX = region.From()[0]; pX<region.To()[0]; ++pX, ++r ) 
		  {
		  Types::DataItem refValue;
		  if ( refChannel->Get( refValue, r ) && (refValue > valueRange[channel].m_LowerBound ) )
		    {
		    /* found pixel over threshold; set flag and terminate loops */
		    active[cp] = true;

		    channel = this->m_NumberOfChannels;
		    pX = region.To()[0];
		    pY = region.To()[1];
		    pZ = region.To()[2];
		    }
		  }
	    }
	  else
	    {
	    const float* fltChannel = &(this->m_ReformattedFloatingChannels[channel-this->m_ReferenceChannels.size()][0]);
	    for ( int pZ = region.From()[2]; pZ<region.To()[2]; ++pZ, r += endOfPlane ) 
	      for ( int pY = region.From()[1]; pY<region.To()[1]; ++pY, r += endOfLine ) 
		for ( int pX = region.From()[0]; pX<region.To()[0]; ++pX, ++r ) 
		  {
		  const float fltValue = fltChannel[r];
		  if ( finite( fltValue ) && (fltValue > valueRange[channel].m_LowerBound) )
		    {
		    /* found pixel over threshold; set flag and terminate loops */
		    active[cp] = true;

		    channel = this->m_NumberOfChannels;
		    pX = region.To()[0];
		    pY = region.To()[1];
		    pZ = region.To()[2];
		    }
		  }
	    }
	  }
	}
      }
    }
  
  size_t inactive = 0;

  for ( size_t cp = 0; cp < numberOfControlPoints; ++cp )
    {
    size_t param = 3 * cp;
    if ( active[cp] )
      {
      for ( size_t dim = 0; dim<3; ++dim, ++param ) 
	{
	this->m_Transformation.SetParameterActive( param );
	this->m_StepScaleVector[param] = this->GetParamStep( param );
	}
      }
    else
      {
      for ( size_t dim = 0; dim<3; ++dim, ++param ) 
	{
	this->m_Transformation.SetParameterInactive( param );
	this->m_StepScaleVector[param] = 0;
	}
      inactive += 3;
      }
    }

  // now fix any fixed coordinate dimensions
  for ( std::list<int>::const_iterator it = this->m_FixedCoordinateDimensions.begin(); 
	it != this->m_FixedCoordinateDimensions.end(); ++it )
    {
    size_t param = *it;
    for ( size_t cp = 0; cp < numberOfControlPoints; ++cp, param += 3 )
      {
      this->m_Transformation.SetParameterInactive( param );
      this->m_StepScaleVector[param] = 0;

      if ( active[cp] )
	++inactive;
      }
    }

  DebugOutput( 1 ).GetStream().printf( "Deactivated %d out of %d control points.\n", (int)inactive / 3, (int)this->ParamVectorDim() / 3 );
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::ContinueMetricStoreReformatted( MetricData& metricData, const size_t rindex, const Vector3D& fvector )
{
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Types::DataItem values[ this->m_NumberOfChannels ];
#else
  std::vector<Types::DataItem> values( this->m_NumberOfChannels );
#endif
  
  size_t idx = 0;
  for ( size_t ref = 0; ref < this->m_ReferenceChannels.size(); ++ref )
    {
    if ( !this->m_ReferenceChannels[ref]->GetDataAt( values[idx++], rindex ) ) return;
    }
  
  for ( size_t flt = 0; flt < this->m_FloatingChannels.size(); ++flt )
    {
    if ( !this->m_FloatingInterpolators[flt]->GetDataAt( fvector, values[idx++] ) )
      {
      for ( size_t f = 0; f < this->m_FloatingChannels.size(); ++f ) this->m_ReformattedFloatingChannels[f][rindex] = MathUtil::GetFloatNaN();
      return;
      }
    }
  
  idx = this->m_ReferenceChannels.size();
  for ( size_t flt = 0; flt < this->m_FloatingChannels.size(); ++flt, ++idx )
    this->m_ReformattedFloatingChannels[flt][rindex] = static_cast<float>( values[idx] );
  
  metricData += &(values[0]);
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::BacktraceMetric( MetricData& metricData, const DataGrid::RegionType& region )
{
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Types::DataItem values[ this->m_NumberOfChannels ];
#else
  std::vector<Types::DataItem> values( this->m_NumberOfChannels );
#endif
  
  for ( int pZ = region.From()[2]; pZ < region.To()[2]; ++pZ ) 
    {
    for ( int pY = region.From()[1]; pY < region.To()[1]; ++pY ) 
      {
      size_t rindex = region.From()[0] + this->m_ReferenceDims[0] * ( pY + this->m_ReferenceDims[1] );
      for ( int pX = region.From()[0]; pX < region.To()[0]; ++pX, ++rindex ) 
	{
	bool allChannelsValid = true;

	size_t idx = 0;
	for ( size_t ref = 0; (ref < this->m_ReferenceChannels.size()) && allChannelsValid; ++ref )
	  {
	  if ( !this->m_ReferenceChannels[ref]->GetDataAt( values[idx++], rindex ) ) 
	    allChannelsValid = false;
	  }
	
	for ( size_t flt = 0; (flt < this->m_FloatingChannels.size()) && allChannelsValid; ++flt, ++idx )
	  {
	  values[idx] = this->m_ReformattedFloatingChannels[flt][rindex];
	  if ( !finite( values[idx] ) ) 
	    allChannelsValid = false;
	  }

	if ( allChannelsValid )
	  {
	  metricData -= &(values[0]);
	  }
	}
      }
    }
}

template<class TMetricFunctional>
typename SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>::ReturnType
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::Evaluate() 
{
  if ( this->m_ReformattedFloatingChannels.size() == 0 )
    {
    this->AllocateReformattedFloatingChannels();
    }

  this->m_MetricData.Init( this );

  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfThreads = threadPool.GetNumberOfThreads();
  const size_t numberOfTasks = 4 * numberOfThreads - 3;

  std::vector< ThreadParameters<Self> > threadParams( numberOfTasks );
  for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx )
    {
    threadParams[taskIdx].thisObject = this;
    }
  threadPool.Run( EvaluateThreadFunction, threadParams );
  
  typename Self::ReturnType costFunction = this->GetMetric( this->m_MetricData );
  if ( this->m_JacobianConstraintWeight > 0 )
    {
    costFunction -= this->m_JacobianConstraintWeight * this->m_Transformation.GetJacobianConstraint();
    }
  return costFunction;
}

template<class TMetricFunctional>
typename SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>::ReturnType
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::EvaluateIncremental
( const SplineWarpXform* transformation, MetricData& metricData, const DataGrid::RegionType& region )
{
  const size_t pixelsPerLineRegion = region.To()[0] - region.From()[0];
  std::vector<Vector3D> pFloating( pixelsPerLineRegion );

  const DataGrid::IndexType& dims = this->m_ReferenceDims;
  const int dimsX = dims[0], dimsY = dims[1];

  for ( int pZ = region.From()[2]; pZ < region.To()[2]; ++pZ ) 
    {
    for ( int pY = region.From()[1]; pY < region.To()[1]; ++pY ) 
      {
      transformation->GetTransformedGridRow( pixelsPerLineRegion, &pFloating[0], region.From()[0], pY, pZ );

      size_t r = region.From()[0] + dimsX * (pY + dimsY * pZ );
      for ( int pX = region.From()[0]; pX < region.To()[0]; ++pX, ++r ) 
	{
	// Continue metric computation.
	this->ContinueMetric( metricData, r, pFloating[pX-region.From()[0]] );
	}
      }
    }
  
  return this->GetMetric( metricData );
}

template<class TMetricFunctional>
typename SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>::ReturnType
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
  const typename Self::ReturnType current = this->EvaluateAt( v );

  // need to call EvaluateAt() first to make sure all reformatted floating channels are up to date.
  if ( this->m_UpdateTransformationFixedControlPointsRequired )
    this->UpdateTransformationFixedControlPoints();

  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfThreads = threadPool.GetNumberOfThreads();
  const size_t numberOfTasks = 4 * numberOfThreads - 3;  

  std::vector< EvaluateGradientThreadParameters > threadParams( numberOfTasks );
  for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx )
    {
    threadParams[taskIdx].thisObject = this;
    threadParams[taskIdx].m_Step = step;
    threadParams[taskIdx].m_ParameterVector = &v;
    threadParams[taskIdx].m_Gradient = g.Elements;
    threadParams[taskIdx].m_MetricBaseValue = current;    
    }
  threadPool.Run( EvaluateWithGradientThreadFunction, threadParams );

  return current;
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::AllocateReformattedFloatingChannels() 
{
  this->m_ReformattedFloatingChannels.resize( this->GetNumberOfFloatingChannels() );
  for ( size_t flt = 0; flt < this->m_ReformattedFloatingChannels.size(); ++flt )
    {
    this->m_ReformattedFloatingChannels[flt].resize( this->m_ReferenceChannels[0]->GetNumberOfPixels() );
    }
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
::ClearReformattedFloatingChannels() 
{
  this->m_ReformattedFloatingChannels.resize( 0 );
}

} // namespace cmtk
