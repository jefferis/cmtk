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

#include "cmtkVoxelMatchingElasticFunctional.h"

#include <System/cmtkDebugOutput.h>

#include <Registration/cmtkVoxelMatchingMutInf.h>
#include <Registration/cmtkVoxelMatchingNormMutInf.h>
#include <Registration/cmtkVoxelMatchingCorrRatio.h>
#include <Registration/cmtkVoxelMatchingMeanSquaredDifference.h>
#include <Registration/cmtkVoxelMatchingCrossCorrelation.h>

#include <Base/cmtkInterpolator.h>
#include <Base/cmtkSplineWarpXform.h>

#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

VoxelMatchingElasticFunctional::VoxelMatchingElasticFunctional
( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
  : VoxelMatchingFunctional( reference, floating ),
    m_ActiveCoordinates( NULL )
{
  Dim = 0;

  this->m_ReferenceDomain = UniformVolume::CoordinateRegionType( UniformVolume::CoordinateVectorType( UniformVolume::CoordinateVectorType::Init( 0 ) ), reference->Size );

  this->m_AdaptiveFixParameters = false;
  this->m_AdaptiveFixThreshFactor = 0.5;

  VectorCache = Memory::ArrayC::Allocate<Vector3D>( ReferenceDims[0] );
  VolumeOfInfluence = NULL;
}

VoxelMatchingElasticFunctional::~VoxelMatchingElasticFunctional()
{
  Memory::ArrayC::Delete( VectorCache );
}

template<class W>
void
VoxelMatchingElasticFunctional_WarpTemplate<W>::WeightedDerivative
( double& lower, double& upper, W& warp, const int param, const Types::Coordinate step ) const
{
  if ( this->m_JacobianConstraintWeight > 0 )
    {
    double lowerConstraint = 0, upperConstraint = 0;
    warp.GetJacobianConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step );
    lower -= this->m_JacobianConstraintWeight * lowerConstraint;
    upper -= this->m_JacobianConstraintWeight * upperConstraint;
    } 

  if ( this->m_RigidityConstraintWeight > 0 )
    {
    double lowerConstraint = 0, upperConstraint = 0;

    if ( this->m_RigidityConstraintMap )
      {
      warp.GetRigidityConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step, this->m_RigidityConstraintMap );
      }
    else
      {
      warp.GetRigidityConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step );
      }
    lower -= this->m_RigidityConstraintWeight * lowerConstraint;
    upper -= this->m_RigidityConstraintWeight * upperConstraint;
    } 
  
  if ( this->m_GridEnergyWeight > 0 ) 
    {
    double lowerEnergy = 0, upperEnergy = 0;
    warp.GetGridEnergyDerivative( lowerEnergy, upperEnergy, param, step );
    lower -= this->m_GridEnergyWeight * lowerEnergy;
    upper -= this->m_GridEnergyWeight * upperEnergy;
    }

  // Catch infinite values that result from a folding grid. Effectively
  // prevent this by setting the gradient term to 0.
  if ( !finite(upper) || !finite(lower) ) 
    {
    lower = upper = 0;
    }
  else
    {
    if ( this->m_MatchedLandmarkList.GetPtr() ) 
      {
      double lowerMSD, upperMSD;
      warp.GetDerivativeLandmarksMSD( lowerMSD, upperMSD, this->m_MatchedLandmarkList, param, step );
      lower -= this->m_LandmarkErrorWeight * lowerMSD;
      upper -= this->m_LandmarkErrorWeight * upperMSD;
      }
    if ( InverseTransformation ) 
      {
      double lowerIC, upperIC;
      warp.GetDerivativeInverseConsistencyError( lowerIC, upperIC, this->InverseTransformation, this->ReferenceGrid, &(this->VolumeOfInfluence[param]), param, step );
      lower -= InverseConsistencyWeight * lowerIC;
      upper -= InverseConsistencyWeight * upperIC;
      }
    }
}

template<class W> 
void
VoxelMatchingElasticFunctional_WarpTemplate<W>::SetWarpXform
( typename W::SmartPtr& warp )
{
  Warp = W::SmartPtr::DynamicCastFrom( warp );
  if ( Warp )
    {
    Warp->RegisterVolume( ReferenceGrid );
    
    if ( Dim != Warp->VariableParamVectorDim() ) 
      {
      if ( VolumeOfInfluence ) 
	Memory::ArrayC::Delete( VolumeOfInfluence );
      Dim = Warp->VariableParamVectorDim();
      this->StepScaleVector.resize( Dim );
      this->VolumeOfInfluence = Memory::ArrayC::Allocate<DataGrid::RegionType>( Dim );
      }
    
    DataGrid::RegionType *VOIptr = this->VolumeOfInfluence;
    for ( size_t dim=0; dim<Dim; ++dim, ++VOIptr ) 
      {
      StepScaleVector[dim] = this->GetParamStep( dim );      
      *VOIptr = this->GetReferenceGridRange( Warp->GetVolumeOfInfluence( dim, this->m_ReferenceDomain ) );
      }
    
    WarpNeedsFixUpdate = true;
  }
}

template<class VM>
void
VoxelMatchingElasticFunctional_Template<VM>::UpdateWarpFixedParameters() 
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
#ifdef _OPENMP
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
#else
      const unsigned int numSamplesX = this->Metric->DataX.NumberOfSamples;
      const Types::DataItemRange rangeX = this->Metric->DataX.GetValueRange();
      const unsigned int numBinsX = JointHistogramBase::CalcNumBins( numSamplesX, rangeX );

      const unsigned int numSamplesY = this->Metric->DataY.NumberOfSamples;
      const Types::DataItemRange rangeY = this->Metric->DataY.GetValueRange();
      const unsigned int numBinsY = JointHistogramBase::CalcNumBins( numSamplesY, rangeY );
      
      this->m_ConsistencyHistogram = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>() );
      
      this->m_ConsistencyHistogram->Resize( numBinsX, numBinsY );
      this->m_ConsistencyHistogram->SetRangeX( rangeX );
      this->m_ConsistencyHistogram->SetRangeY( rangeY );
#endif
      
#pragma omp parallel for
      for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
	{
#ifdef _OPENMP
	JointHistogram<unsigned int>& threadHistogram = *(this->m_ThreadConsistencyHistograms[ omp_get_thread_num() ]);
#else
	JointHistogram<unsigned int>& threadHistogram = *(this->m_ConsistencyHistogram);
#endif
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

VoxelMatchingElasticFunctional* 
CreateElasticFunctional
( const int metric,
  UniformVolume::SmartPtr& refVolume, 
  UniformVolume::SmartPtr& fltVolume )
{
  switch ( fltVolume->GetData()->GetDataClass() ) 
    {
    case DATACLASS_UNKNOWN :
    case DATACLASS_GREY :
      switch ( metric ) 
	{
	case 0:
	  return new VoxelMatchingElasticFunctional_Template< VoxelMatchingNormMutInf_Trilinear>( refVolume, fltVolume );
	case 1:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_Trilinear>( refVolume, fltVolume );
	case 2:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_Trilinear>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMeanSquaredDifference>( refVolume, fltVolume );
	case 5:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCrossCorrelation>( refVolume, fltVolume );
	default:
	  return NULL;
	}
    case DATACLASS_LABEL:    
      switch ( metric ) 
	{
	case 0:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingNormMutInf_NearestNeighbor>( refVolume, fltVolume );
	case 1:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_NearestNeighbor>( refVolume, fltVolume );
	case 2:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_NearestNeighbor>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMeanSquaredDifference>( refVolume, fltVolume );
	case 5:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCrossCorrelation>( refVolume, fltVolume );
      default:
	return NULL;
	}
    }
  
  return NULL;
}

template class VoxelMatchingElasticFunctional_WarpTemplate<SplineWarpXform>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingNormMutInf_Trilinear>;
template class VoxelMatchingElasticFunctional_Template<VoxelMatchingNormMutInf_NearestNeighbor>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_Trilinear>;
template class VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_NearestNeighbor>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_Trilinear>;
template class VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_NearestNeighbor>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingCrossCorrelation>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingMeanSquaredDifference>;

} // namespace cmtk
