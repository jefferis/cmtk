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

#include <omp.h>

#include <System/cmtkDebugOutput.h>

template<class VM> 
void
cmtk::ParallelElasticFunctional<VM>::UpdateWarpFixedParameters() 
{
  int numCtrlPoints = this->Dim / 3;
  
  std::vector<double> mapRef( numCtrlPoints );
  std::vector<double> mapMod( numCtrlPoints );

  int inactive = 0;

  const typename VM::Exchange unsetY = this->Metric->DataY.padding();

  if ( this->ReferenceDataClass == DATACLASS_LABEL ) 
    {
    if ( this->m_ActiveCoordinates )
      this->Warp->SetParametersActive( this->m_ActiveCoordinates );
    else
      this->Warp->SetParametersActive();

#pragma omp parallel for reduction(+:inactive)
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      /// We cannot use the precomputed table of VOIs here because in "fast" mode, these VOIs are smaller than we want them here.
      const DataGrid::RegionType voi = this->GetReferenceGridRange( this->Warp->GetVolumeOfInfluence( 3 * ctrl, this->m_ReferenceDomain, false /* disable fast mode */ ) );
      
      int r = voi.From()[0] + this->DimsX * ( voi.From()[1] + this->DimsY * voi.From()[2] );
      
      bool active = false;
      for ( int pZ = voi.From()[2]; (pZ < voi.To()[2]) && !active; ++pZ ) 
	{
	for ( int pY = voi.From()[1]; (pY < voi.To()[1]) && !active; ++pY ) 
	  {
	  for ( int pX = voi.From()[0]; (pX < voi.To()[0]); ++pX, ++r ) 
	    {
	    if ( ( this->Metric->GetSampleX( r ) != 0 ) || ( ( this->WarpedVolume[r] != unsetY ) && ( this->WarpedVolume[r] != 0 ) ) ) 
	      {
	      active = true;
	      break;
	      }
	    }
	  r += ( voi.From()[0] + ( this->DimsX-voi.To()[0] ) );
	  }
	r += this->DimsX * ( voi.From()[1] + ( this->DimsY-voi.To()[1] ) );
	}
      
      if ( !active ) 
	{
	inactive += 3;
	
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  this->Warp->SetParameterInactive( dim );
	  }
	}
      }
    } 
  else
    {
    if ( this->m_ThreadConsistencyHistograms.size() != this->m_NumberOfThreads )
      {
      this->m_ThreadConsistencyHistograms.resize( this->m_NumberOfThreads );
      
      const unsigned int numSamplesX = this->Metric->DataX.NumberOfSamples;
      const Types::DataItemRange rangeX = this->Metric->DataX.GetValueRange();
      const unsigned int numBinsX = JointHistogramBase::CalcNumBins( numSamplesX, rangeX );

      const unsigned int numSamplesY = this->Metric->DataY.NumberOfSamples;
      const Types::DataItemRange rangeY = this->Metric->DataY.GetValueRange();
      const unsigned int numBinsY = JointHistogramBase::CalcNumBins( numSamplesY, rangeY );
      
      for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
	{
	this->m_ThreadConsistencyHistograms[thread] = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>() );
	
	this->m_ThreadConsistencyHistograms[thread]->Resize( numBinsX, numBinsY );
	this->m_ThreadConsistencyHistograms[thread]->SetRangeX( rangeX );
	this->m_ThreadConsistencyHistograms[thread]->SetRangeY( rangeY );
	}
      }

#pragma omp parallel for
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      JointHistogram<unsigned int>& threadHistogram = *(this->m_ThreadConsistencyHistograms[ omp_get_thread_num() ]);
      threadHistogram.Reset();
      
      // We cannot use the precomputed table of VOIs here because in "fast" mode, these VOIs are smaller than we want them here.
      const DataGrid::RegionType voi = this->GetReferenceGridRange( this->Warp->GetVolumeOfInfluence( 3 * ctrl, this->m_ReferenceDomain, false /* disable fast mode */ ) );
      
      int r = voi.From()[0] + this->DimsX * ( voi.From()[1] + this->DimsY * voi.From()[2] );
      
      const int endOfLine = ( voi.From()[0] + ( this->DimsX-voi.To()[0]) );
      const int endOfPlane = this->DimsX * ( voi.From()[1] + (this->DimsY-voi.To()[1]) );
      
      for ( int pZ = voi.From()[2]; pZ<voi.To()[2]; ++pZ ) 
	{
	for ( int pY = voi.From()[1]; pY<voi.To()[1]; ++pY ) 
	  {
	  for ( int pX = voi.From()[0]; pX<voi.To()[0]; ++pX, ++r ) 
	    {
	    // Continue metric computation.
	    if ( this->WarpedVolume[r] != unsetY ) 
	      {
	      threadHistogram.Increment( threadHistogram.ValueToBinX( this->Metric->GetSampleX( r ) ), threadHistogram.ValueToBinY( this->WarpedVolume[r] ) );
	      }
	    }
	  r += endOfLine;
	  }
	r += endOfPlane;
	}
      threadHistogram.GetMarginalEntropies( mapRef[ctrl], mapMod[ctrl] );
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
      this->Warp->SetParametersActive( this->m_ActiveCoordinates );
    else
      this->Warp->SetParametersActive();
          
    for ( int ctrl=0; ctrl<numCtrlPoints; ++ctrl ) 
      {
      if (  ( mapRef[ctrl] < refThresh ) && ( mapMod[ctrl] < modThresh ) ) 
	{
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  this->Warp->SetParameterInactive( dim );
	  }
	inactive += 3;
	}
      }
    }
  
  for ( size_t idx = 0; idx < this->Dim; ++idx ) 
    {

    if ( this->Warp->GetParameterActive( idx ) )
      {
      this->StepScaleVector[idx] = this->GetParamStep( idx );
      }
    else
      {
      this->StepScaleVector[idx] = 0;
      }
    }
  
  DebugOutput( 1 ).GetStream().printf( "Deactivated %d out of %d parameters.\n", inactive, (int)this->Dim );
  
  this->WarpNeedsFixUpdate = false;
}
