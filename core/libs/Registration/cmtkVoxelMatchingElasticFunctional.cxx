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

#include <cmtkVoxelMatchingElasticFunctional.h>

#ifdef CMTK_BUILD_SMP
#  include <cmtkParallelElasticFunctional.h>
#endif
#include <cmtkSplineWarpXform.h>

#include <cmtkVoxelMatchingMutInf.h>
#include <cmtkVoxelMatchingNormMutInf.h>
#include <cmtkVoxelMatchingCorrRatio.h>
#include <cmtkVoxelMatchingMeanSquaredDifference.h>
#include <cmtkVoxelMatchingCrossCorrelation.h>

#include <cmtkInterpolator.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

VoxelMatchingElasticFunctional::VoxelMatchingElasticFunctional
( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
  : VoxelMatchingFunctional( reference, floating )
{
  Dim = 0;

  ReferenceFrom.Set( 0,0,0 );
  ReferenceTo.Set( reference->Size );

  AdaptiveFixParameters = false;
  AdaptiveFixThreshFactor = 0.5;
  StepScaleVector = NULL;

  VectorCache = Memory::AllocateArray<Vector3D>( ReferenceDims[0] );
  VolumeOfInfluence = NULL;
}

VoxelMatchingElasticFunctional::~VoxelMatchingElasticFunctional()
{
  if ( VectorCache ) delete[] VectorCache;
  if ( StepScaleVector ) delete[] StepScaleVector;
}

template<class W>
void
VoxelMatchingElasticFunctional_WarpTemplate<W>::WeightedDerivative
( double& lower, double& upper, typename W::SmartPtr& warp, 
  const int param, const Types::Coordinate step ) const
{
  if ( JacobianConstraintWeight > 0 )
    {
    double lowerConstraint = 0, upperConstraint = 0;
    warp->GetJacobianConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step );
    lower -= JacobianConstraintWeight * lowerConstraint;
    upper -= JacobianConstraintWeight * upperConstraint;
    } 

  if ( RigidityConstraintWeight > 0 )
    {
    double lowerConstraint = 0, upperConstraint = 0;

    if ( this->RigidityConstraintMap )
      {
      warp->GetRigidityConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step,
					     this->RigidityConstraintMap );
      }
    else
      {
      warp->GetRigidityConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step );
      }
    lower -= RigidityConstraintWeight * lowerConstraint;
    upper -= RigidityConstraintWeight * upperConstraint;
    } 
  
  if ( GridEnergyWeight > 0 ) 
    {
    double lowerEnergy = 0, upperEnergy = 0;
    warp->GetGridEnergyDerivative( lowerEnergy, upperEnergy, param, step );
    lower -= GridEnergyWeight * lowerEnergy;
    upper -= GridEnergyWeight * upperEnergy;
    }

  // Catch infinite values that result from a folding grid. Effectively
  // prevent this by setting the gradient term to 0.
  if ( !finite(upper) || !finite(lower) ) 
    {
    lower = upper = 0;
    }
  else
    {
    if ( MatchedLandmarkList.GetPtr() ) 
      {
      double lowerMSD, upperMSD;
      warp->GetDerivativeLandmarksMSD( lowerMSD, upperMSD, MatchedLandmarkList, param, step );
      lower -= LandmarkErrorWeight * lowerMSD;
      upper -= LandmarkErrorWeight * upperMSD;
      }
    if ( InverseTransformation ) 
      {
      double lowerIC, upperIC;
      warp->GetDerivativeInverseConsistencyError
	( lowerIC, upperIC, this->InverseTransformation, this->ReferenceGrid, &(this->VolumeOfInfluence[param]), param, step );
      lower -= InverseConsistencyWeight * lowerIC;
      upper -= InverseConsistencyWeight * upperIC;
      }
    }
}

template<class W> 
void
VoxelMatchingElasticFunctional_WarpTemplate<W>::SetWarpXform
( WarpXform::SmartPtr& warp )
{
  Warp = W::SmartPtr::DynamicCastFrom( warp );
  if ( Warp )
    {
    Warp->RegisterVolume( ReferenceGrid );
    Warp->SetIncompressibilityMap( IncompressibilityMap );
    
    if ( Dim != Warp->VariableParamVectorDim() ) 
      {
      if ( StepScaleVector ) delete[] StepScaleVector;
      if ( VolumeOfInfluence ) delete[] VolumeOfInfluence;
      Dim = Warp->VariableParamVectorDim();
      StepScaleVector = Memory::AllocateArray<Types::Coordinate>( Dim );
      VolumeOfInfluence = Memory::AllocateArray<Rect3D>( Dim );
      }
    
    Rect3D *VOIptr = VolumeOfInfluence;
    Vector3D fromVOI, toVOI;
    for ( size_t dim=0; dim<Dim; ++dim, ++VOIptr ) 
      {
      StepScaleVector[dim] = this->GetParamStep( dim );
      Warp->GetVolumeOfInfluence( dim, ReferenceFrom, ReferenceTo, fromVOI, toVOI );
      this->GetReferenceGridRange( fromVOI, toVOI, *VOIptr );
      }
    
    WarpNeedsFixUpdate = true;
  }
}

template<class VM, class W>
void
VoxelMatchingElasticFunctional_Template<VM,W>::UpdateWarpFixedParameters() 
{
  if ( this->ConsistencyHistogram.IsNull() ) 
    {
    this->ConsistencyHistogram = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>() );
    unsigned int numSamplesX = this->Metric->DataX.NumberOfSamples;
    Types::DataItem fromX, toX;
    this->Metric->DataX.GetValueRange( fromX, toX );
    unsigned int numBinsX = this->ConsistencyHistogram->CalcNumBins( numSamplesX, fromX, toX );
    
    unsigned int numSamplesY = this->Metric->DataY.NumberOfSamples;
    Types::DataItem fromY, toY;
    this->Metric->DataY.GetValueRange( fromY, toY );
    unsigned int numBinsY = this->ConsistencyHistogram->CalcNumBins( numSamplesY, fromY, toY );
    
    this->ConsistencyHistogram->SetNumBins( numBinsX, numBinsY );
    this->ConsistencyHistogram->SetRangeX( fromX, toX );
    this->ConsistencyHistogram->SetRangeY( fromY, toY );
    }
  
  int numCtrlPoints = this->Dim / 3;
  
  double *mapRef = Memory::AllocateArray<double>( numCtrlPoints );
  double *mapMod = Memory::AllocateArray<double>( numCtrlPoints );

  Rect3D voi;
  Vector3D fromVOI, toVOI;
  int pX, pY, pZ;

  int inactive = 0;

  const typename VM::Exchange unsetY = this->Metric->DataY.padding();

  if ( this->ReferenceDataClass == DATACLASS_LABEL ) 
    {
    this->Warp->SetParameterActive();
    
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      /// We cannot use the precomputed table of VOIs here because in "fast"
      /// mode, these VOIs are smaller than we want them here.
      this->Warp->GetVolumeOfInfluence( 3 * ctrl, this->ReferenceFrom, this->ReferenceTo, fromVOI, toVOI, 0 );
      this->GetReferenceGridRange( fromVOI, toVOI, voi );
      
      int r = voi.startX + this->DimsX * ( voi.startY + this->DimsY * voi.startZ );
      
      bool active = false;
      for ( pZ = voi.startZ; (pZ < voi.endZ) && !active; ++pZ ) 
	{
	for ( pY = voi.startY; (pY < voi.endY) && !active; ++pY ) 
	  {
	  for ( pX = voi.startX; (pX < voi.endX); ++pX, ++r ) 
	    {
	    if ( ( this->Metric->GetSampleX( r ) != 0 ) || ( ( this->WarpedVolume[r] != unsetY ) && ( this->WarpedVolume[r] != 0 ) ) ) 
	      {
	      active = true;
	      break;
	      }
	    }
	  r += ( voi.startX + ( this->DimsX-voi.endX ) );
	  }
	r += this->DimsX * ( voi.startY + ( this->DimsY-voi.endY ) );
	}
      
      if ( !active ) 
	{
	inactive += 3;
	
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  this->Warp->SetParameterInactive( dim );
	  this->StepScaleVector[dim] = this->GetParamStep( dim );
	  }
	}
      }
    } 
  else
    {
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      this->ConsistencyHistogram->Reset();
      
      // We cannot use the precomputed table of VOIs here because in "fast"
      // mode, these VOIs are smaller than we want them here.
      this->Warp->GetVolumeOfInfluence( 3 * ctrl, this->ReferenceFrom, this->ReferenceTo, fromVOI, toVOI, 0 );
      this->GetReferenceGridRange( fromVOI, toVOI, voi );
      
      int r = voi.startX + this->DimsX * ( voi.startY + this->DimsY * voi.startZ );
      
      const int endOfLine = ( voi.startX + ( this->DimsX-voi.endX) );
      const int endOfPlane = this->DimsX * ( voi.startY + (this->DimsY-voi.endY) );
      
      for ( pZ = voi.startZ; pZ<voi.endZ; ++pZ ) 
	{
	for ( pY = voi.startY; pY<voi.endY; ++pY ) 
	  {
	  for ( pX = voi.startX; pX<voi.endX; ++pX, ++r ) 
	    {
	    // Continue metric computation.
	    if ( WarpedVolume[r] != unsetY ) 
	      {
	      this->ConsistencyHistogram->Increment
		( this->ConsistencyHistogram->ValueToBinX( this->Metric->GetSampleX( r ) ), 
		  this->ConsistencyHistogram->ValueToBinY( this->WarpedVolume[r] ) );
	      }
	    }
	  r += endOfLine;
	  }
	r += endOfPlane;
	}
      this->ConsistencyHistogram->GetMarginalEntropies( mapRef[ctrl], mapMod[ctrl] );
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
    
    const double refThresh = refMin + this->AdaptiveFixThreshFactor * (refMax - refMin);
    const double modThresh = modMin + this->AdaptiveFixThreshFactor * (modMax - modMin);
      
    this->Warp->SetParameterActive();
      
    for ( int ctrl=0; ctrl<numCtrlPoints; ++ctrl ) 
      {
      if (  ( mapRef[ctrl] < refThresh ) && ( mapMod[ctrl] < modThresh ) ) 
	{
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  this->Warp->SetParameterInactive( dim );
	  this->StepScaleVector[dim] = this->GetParamStep( dim );
	  }
	inactive += 3;
	}
      }
    }
  
  fprintf( stderr, "Deactivated %d out of %d parameters.\n", inactive, (int)this->Dim );
  
  delete[] mapRef;
  delete[] mapMod;

  this->WarpNeedsFixUpdate = false;
}

VoxelMatchingElasticFunctional* 
CreateElasticFunctional
( const int metric,
  UniformVolume::SmartPtr& refVolume, 
  UniformVolume::SmartPtr& fltVolume )
{
#ifdef CMTK_BUILD_SMP
  switch ( fltVolume->GetData()->GetDataClass() ) 
    {
    case DATACLASS_UNKNOWN :
    case DATACLASS_BINARY :
      // not a lot we can do here.
      return NULL;
    case DATACLASS_GREY :
      switch ( metric ) 
	{
	case 0:
	  return new ParallelElasticFunctional< VoxelMatchingNormMutInf_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 1:
	  return new ParallelElasticFunctional<VoxelMatchingMutInf_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 2:
	  return new ParallelElasticFunctional<VoxelMatchingCorrRatio_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new ParallelElasticFunctional<VoxelMatchingMeanSquaredDifference,SplineWarpXform>( refVolume, fltVolume );
	case 5:
	  return new ParallelElasticFunctional<VoxelMatchingCrossCorrelation,SplineWarpXform>( refVolume, fltVolume );
	default:
	  return NULL;
	}
    case DATACLASS_LABEL:
      switch ( metric ) 
	{
	case 0:
	  return new ParallelElasticFunctional< VoxelMatchingNormMutInf_NearestNeighbor, SplineWarpXform >( refVolume, fltVolume );
	case 1:
	  return new ParallelElasticFunctional<VoxelMatchingMutInf_NearestNeighbor,SplineWarpXform>( refVolume, fltVolume );
	case 2:
	  return new ParallelElasticFunctional<VoxelMatchingCorrRatio_NearestNeighbor,SplineWarpXform>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new ParallelElasticFunctional<VoxelMatchingMeanSquaredDifference,SplineWarpXform>( refVolume, fltVolume );
	case 5:
	  return new ParallelElasticFunctional<VoxelMatchingCrossCorrelation,SplineWarpXform>( refVolume, fltVolume );
	default:
	  return NULL;
	}
    }

#else // single processing

  switch ( fltVolume->GetData()->GetDataClass() ) 
    {
    case DATACLASS_UNKNOWN :
    case DATACLASS_BINARY :
      // not a lot we can do here.
      return NULL;
    case DATACLASS_GREY :
      switch ( metric ) 
	{
	case 0:
	  return new VoxelMatchingElasticFunctional_Template< VoxelMatchingNormMutInf_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 1:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 2:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMeanSquaredDifference,SplineWarpXform>( refVolume, fltVolume );
	case 5:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCrossCorrelation,SplineWarpXform>( refVolume, fltVolume );
	default:
	  return NULL;
	}
    case DATACLASS_LABEL:    
      switch ( metric ) 
	{
	case 0:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingNormMutInf_NearestNeighbor, SplineWarpXform>( refVolume, fltVolume );
	case 1:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_NearestNeighbor,SplineWarpXform>( refVolume, fltVolume );
	case 2:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_NearestNeighbor,SplineWarpXform>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingMeanSquaredDifference,SplineWarpXform>( refVolume, fltVolume );
	case 5:
	  return new VoxelMatchingElasticFunctional_Template<VoxelMatchingCrossCorrelation,SplineWarpXform>( refVolume, fltVolume );
      default:
	return NULL;
	}
    }
#endif
  
  return NULL;
}

template class VoxelMatchingElasticFunctional_WarpTemplate<SplineWarpXform>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingNormMutInf_Trilinear, SplineWarpXform>;
template class VoxelMatchingElasticFunctional_Template<VoxelMatchingNormMutInf_NearestNeighbor, SplineWarpXform>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_Trilinear, SplineWarpXform>;
template class VoxelMatchingElasticFunctional_Template<VoxelMatchingMutInf_NearestNeighbor, SplineWarpXform>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_Trilinear, SplineWarpXform>;
template class VoxelMatchingElasticFunctional_Template<VoxelMatchingCorrRatio_NearestNeighbor, SplineWarpXform>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingCrossCorrelation, SplineWarpXform>;

template class VoxelMatchingElasticFunctional_Template<VoxelMatchingMeanSquaredDifference, SplineWarpXform>;

} // namespace cmtk
