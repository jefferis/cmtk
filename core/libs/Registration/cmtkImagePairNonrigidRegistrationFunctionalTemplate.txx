/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>

#include <System/cmtkDebugOutput.h>

#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

template<class VM>
void
cmtk::ImagePairNonrigidRegistrationFunctionalTemplate<VM>::MatchRefFltIntensities()
{
  const Types::DataItem paddingValue = DataTypeTraits<Types::DataItem>::ChoosePaddingValue();
  TypedArray::SmartPtr warpedArray( TypedArray::Create( TYPE_ITEM, this->m_WarpedVolume, this->m_FloatingGrid->GetNumberOfPixels(), true /*padding*/, &paddingValue ) );

  UniformVolume::SmartPtr floatingCopy( this->m_FloatingGrid->Clone() );
  floatingCopy->GetData()->ApplyFunctionObject( TypedArrayFunctionHistogramMatching( *warpedArray, *(this->m_ReferenceGrid->GetData()) ) );
  this->m_Metric->SetFloatingVolume( floatingCopy );
}

template<class VM>
void
cmtk::ImagePairNonrigidRegistrationFunctionalTemplate<VM>::UpdateWarpFixedParameters() 
{
  int numCtrlPoints = this->Dim / 3;
  
  std::vector<double> mapRef( numCtrlPoints );
  std::vector<double> mapMod( numCtrlPoints );

  int inactive = 0;

  const Types::DataItem unsetY = DataTypeTraits<Types::DataItem>::ChoosePaddingValue();

  if ( this->m_ReferenceDataClass == DATACLASS_LABEL ) 
    {
    if ( this->m_ActiveCoordinates )
      this->m_Warp->SetParametersActive( this->m_ActiveCoordinates );
    else
      this->m_Warp->SetParametersActive();
    
#pragma omp parallel for reduction(+:inactive)
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      /// We cannot use the precomputed table of VOIs here because in "fast" mode, these VOIs are smaller than we want them here.
      const DataGrid::RegionType voi = this->GetReferenceGridRange( this->m_Warp->GetVolumeOfInfluence( 3 * ctrl, this->m_ReferenceDomain, false /*disable fast mode*/ ) );
      
      int r = voi.From()[0] + this->m_DimsX * ( voi.From()[1] + this->m_DimsY * voi.From()[2] );
      
      bool active = false;
      for ( int pZ = voi.From()[2]; (pZ < voi.To()[2]) && !active; ++pZ ) 
	{
	for ( int pY = voi.From()[1]; (pY < voi.To()[1]) && !active; ++pY ) 
	  {
	  for ( int pX = voi.From()[0]; (pX < voi.To()[0]); ++pX, ++r ) 
	    {
	    if ( ( this->m_Metric->GetSampleX( r ) != 0 ) || ( ( this->m_WarpedVolume[r] != unsetY ) && ( this->m_WarpedVolume[r] != 0 ) ) ) 
	      {
	      active = true;
	      break;
	      }
	    }
	  r += ( voi.From()[0] + ( this->m_DimsX-voi.To()[0] ) );
	  }
	r += this->m_DimsX * ( voi.From()[1] + ( this->m_DimsY-voi.To()[1] ) );
	}
      
      if ( !active ) 
	{
	inactive += 3;
	
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  this->m_Warp->SetParameterInactive( dim );
	  }
	}
      }
    } 
  else
    {
#ifdef _OPENMP
    if ( this->m_ThreadConsistencyHistograms.size() != this->m_NumberOfThreads )
      {
      this->m_ThreadConsistencyHistograms.resize( this->m_NumberOfThreads );
      
      const unsigned int numSamplesX = this->m_Metric->GetNumberOfSamplesX();
      const Types::DataItemRange rangeX = this->m_Metric->GetDataRangeX();
      const unsigned int numBinsX = JointHistogramBase::CalcNumBins( numSamplesX, rangeX );
      
      const unsigned int numSamplesY = this->m_Metric->GetNumberOfSamplesY();
      const Types::DataItemRange rangeY = this->m_Metric->GetDataRangeY();
      const unsigned int numBinsY = JointHistogramBase::CalcNumBins( numSamplesY, rangeY );

      for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
	{
	this->m_ThreadConsistencyHistograms[thread] = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>() );
	
	this->m_ThreadConsistencyHistograms[thread]->Resize( numBinsX, numBinsY );
	this->m_ThreadConsistencyHistograms[thread]->SetRangeX( rangeX );
	this->m_ThreadConsistencyHistograms[thread]->SetRangeY( rangeY );
	}
      }
#else
    if ( !this->m_ConsistencyHistogram )
      {
      this->m_ConsistencyHistogram = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>() );
      const unsigned int numSamplesX = this->m_Metric->GetNumberOfSamplesX();
      const Types::DataItemRange rangeX = this->m_Metric->GetDataRangeX();
      const unsigned int numBinsX = JointHistogramBase::CalcNumBins( numSamplesX, rangeX );
      
      const unsigned int numSamplesY = this->m_Metric->GetNumberOfSamplesY();
      const Types::DataItemRange rangeY = this->m_Metric->GetDataRangeY();
      const unsigned int numBinsY = JointHistogramBase::CalcNumBins( numSamplesY, rangeY );
      
      this->m_ConsistencyHistogram->Resize( numBinsX, numBinsY );
      this->m_ConsistencyHistogram->SetRangeX( rangeX );
      this->m_ConsistencyHistogram->SetRangeY( rangeY );
      }
#endif

#pragma omp parallel for    
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
#ifdef _OPENMP
      JointHistogram<unsigned int>& histogram = *(this->m_ThreadConsistencyHistograms[ omp_get_thread_num() ]);
#else
      JointHistogram<unsigned int>& histogram = *(this->m_ConsistencyHistogram);
#endif
      histogram.Reset();
      
      // We cannot use the precomputed table of VOIs here because in "fast" mode, these VOIs are smaller than we want them here.
      const DataGrid::RegionType voi = this->GetReferenceGridRange( this->m_Warp->GetVolumeOfInfluence( 3 * ctrl, this->m_ReferenceDomain, false /*disable fast mode*/ ) );
      
      int r = voi.From()[0] + this->m_DimsX * ( voi.From()[1] + this->m_DimsY * voi.From()[2] );
      
      const int endOfLine = ( voi.From()[0] + ( this->m_DimsX-voi.To()[0]) );
      const int endOfPlane = this->m_DimsX * ( voi.From()[1] + (this->m_DimsY-voi.To()[1]) );
      
      for ( int pZ = voi.From()[2]; pZ<voi.To()[2]; ++pZ ) 
	{
	for ( int pY = voi.From()[1]; pY<voi.To()[1]; ++pY ) 
	  {
	  for ( int pX = voi.From()[0]; pX<voi.To()[0]; ++pX, ++r ) 
	    {
	    // Continue metric computation.
	    if ( this->m_WarpedVolume[r] != unsetY ) 
	      {
	      histogram.Increment( histogram.ValueToBinX( this->m_Metric->GetSampleX( r ) ), histogram.ValueToBinY( this->m_WarpedVolume[r] ) );
	      }
	    }
	  r += endOfLine;
	  }
	r += endOfPlane;
	}
      histogram.GetMarginalEntropies( mapRef[ctrl], mapMod[ctrl] );
      }
    
    double refMin = HUGE_VAL, refMax = -HUGE_VAL;
    double modMin = HUGE_VAL, modMax = -HUGE_VAL;
    for ( int ctrl=0; ctrl<numCtrlPoints; ++ctrl ) 
      {
      if ( mapRef[ctrl] < refMin ) refMin = mapRef[ctrl];
      if ( mapRef[ctrl] > refMax ) refMax = mapRef[ctrl];
      if ( mapMod[ctrl] < modMin ) modMin = mapMod[ctrl];
      if ( mapMod[ctrl] > modMax ) modMax = mapMod[ctrl];
      }
    
    const double refThresh = refMin + this->m_AdaptiveFixThreshFactor * (refMax - refMin);
    const double modThresh = modMin + this->m_AdaptiveFixThreshFactor * (modMax - modMin);
    
    if ( this->m_ActiveCoordinates )
      this->m_Warp->SetParametersActive( this->m_ActiveCoordinates );
    else
      this->m_Warp->SetParametersActive();
      
    for ( int ctrl=0; ctrl<numCtrlPoints; ++ctrl ) 
      {
      if (  ( mapRef[ctrl] < refThresh ) && ( mapMod[ctrl] < modThresh ) ) 
	{
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  this->m_Warp->SetParameterInactive( dim );
	  }
	inactive += 3;
	}
      }
    }

  for ( size_t idx = 0; idx < this->Dim; ++idx ) 
    {
    if ( this->m_Warp->GetParameterActive( idx ) )
      {
      this->m_StepScaleVector[idx] = this->GetParamStep( idx );
      }
    else
      {
      this->m_StepScaleVector[idx] = 0;
      }
    }
  
  DebugOutput( 1 ).GetStream().printf( "Deactivated %d out of %d parameters.\n", inactive, (int)this->Dim );
  
  this->WarpNeedsFixUpdate = false;
}

