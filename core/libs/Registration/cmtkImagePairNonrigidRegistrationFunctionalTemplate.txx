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

#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>

#include <vector>

template<class VM>
void
cmtk::ImagePairNonrigidRegistrationFunctionalTemplate<VM>::MatchRefFltIntensities()
{
  const Types::DataItem paddingValue = DataTypeTraits<Types::DataItem>::ChoosePaddingValue();
  TypedArray::SmartPtr warpedArray( TypedArray::Create( TYPE_ITEM, this->m_WarpedVolume, this->m_FloatingGrid->GetNumberOfPixels(), false /*freeArray*/, true /*padding*/, &paddingValue ) );

  UniformVolume::SmartPtr floatingCopy( this->m_FloatingGrid->Clone( true /*copyData*/ ) );
  floatingCopy->GetData()->ApplyFunctionObject( TypedArrayFunctionHistogramMatching( *warpedArray, *(this->m_ReferenceGrid->GetData()) ) );
  this->m_Metric->SetFloatingVolume( floatingCopy );
}

template<class VM>
void
cmtk::ImagePairNonrigidRegistrationFunctionalTemplate<VM>::UpdateWarpFixedParameters() 
{
  if ( !this->m_ConsistencyHistogram ) 
    {
    this->m_ConsistencyHistogram = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>() );
    const unsigned int numSamplesX = this->m_Metric->GetNumberOfSamplesX();
    const Types::DataItemRange rangeX = this->m_Metric->GetDataRangeX();
    const unsigned int numBinsX = this->m_ConsistencyHistogram->CalcNumBins( numSamplesX, rangeX );
    
    const unsigned int numSamplesY = this->m_Metric->GetNumberOfSamplesY();
    const Types::DataItemRange rangeY = this->m_Metric->GetDataRangeY();
    const unsigned int numBinsY = this->m_ConsistencyHistogram->CalcNumBins( numSamplesY, rangeY );
    
    this->m_ConsistencyHistogram->Resize( numBinsX, numBinsY );
    this->m_ConsistencyHistogram->SetRangeX( rangeX );
    this->m_ConsistencyHistogram->SetRangeY( rangeY );
    }
  
  int numCtrlPoints = this->Dim / 3;
  
  std::vector<double> mapRef( numCtrlPoints );
  std::vector<double> mapMod( numCtrlPoints );

  DataGrid::RegionType voi;
  Vector3D fromVOI, toVOI;
  int pX, pY, pZ;

  int inactive = 0;

  const Types::DataItem unsetY = DataTypeTraits<Types::DataItem>::ChoosePaddingValue();

  if ( this->m_ReferenceDataClass == DATACLASS_LABEL ) 
    {
    this->m_Warp->SetParameterActive();
    
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      /// We cannot use the precomputed table of VOIs here because in "fast"
      /// mode, these VOIs are smaller than we want them here.
      this->m_Warp->GetVolumeOfInfluence( 3 * ctrl, this->ReferenceFrom, this->ReferenceTo, fromVOI, toVOI, 0 );
      voi = this->GetReferenceGridRange( fromVOI, toVOI );
      
      int r = voi.From()[0] + this->m_DimsX * ( voi.From()[1] + this->m_DimsY * voi.From()[2] );
      
      bool active = false;
      for ( pZ = voi.From()[2]; (pZ < voi.To()[2]) && !active; ++pZ ) 
	{
	for ( pY = voi.From()[1]; (pY < voi.To()[1]) && !active; ++pY ) 
	  {
	  for ( pX = voi.From()[0]; (pX < voi.To()[0]); ++pX, ++r ) 
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
	  this->m_StepScaleVector[dim] = this->GetParamStep( dim );
	  }
	}
      }
    } 
  else
    {
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      this->m_ConsistencyHistogram->Reset();
      
      // We cannot use the precomputed table of VOIs here because in "fast"
      // mode, these VOIs are smaller than we want them here.
      this->m_Warp->GetVolumeOfInfluence( 3 * ctrl, this->ReferenceFrom, this->ReferenceTo, fromVOI, toVOI, 0 );
      voi = this->GetReferenceGridRange( fromVOI, toVOI );
      
      int r = voi.From()[0] + this->m_DimsX * ( voi.From()[1] + this->m_DimsY * voi.From()[2] );
      
      const int endOfLine = ( voi.From()[0] + ( this->m_DimsX-voi.To()[0]) );
      const int endOfPlane = this->m_DimsX * ( voi.From()[1] + (this->m_DimsY-voi.To()[1]) );
      
      for ( pZ = voi.From()[2]; pZ<voi.To()[2]; ++pZ ) 
	{
	for ( pY = voi.From()[1]; pY<voi.To()[1]; ++pY ) 
	  {
	  for ( pX = voi.From()[0]; pX<voi.To()[0]; ++pX, ++r ) 
	    {
	    // Continue metric computation.
	    if ( this->m_WarpedVolume[r] != unsetY ) 
	      {
	      this->m_ConsistencyHistogram->Increment( this->m_ConsistencyHistogram->ValueToBinX( this->m_Metric->GetSampleX( r ) ), this->m_ConsistencyHistogram->ValueToBinY( this->m_WarpedVolume[r] ) );
	      }
	    }
	  r += endOfLine;
	  }
	r += endOfPlane;
	}
      this->m_ConsistencyHistogram->GetMarginalEntropies( mapRef[ctrl], mapMod[ctrl] );
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
      
    this->m_Warp->SetParameterActive();
      
    for ( int ctrl=0; ctrl<numCtrlPoints; ++ctrl ) 
      {
      if (  ( mapRef[ctrl] < refThresh ) && ( mapMod[ctrl] < modThresh ) ) 
	{
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  this->m_Warp->SetParameterInactive( dim );
	  this->m_StepScaleVector[dim] = this->GetParamStep( dim );
	  }
	inactive += 3;
	}
      }
    }
  
  fprintf( stderr, "Deactivated %d out of %d parameters.\n", inactive, (int)this->Dim );
  
  this->WarpNeedsFixUpdate = false;
}

