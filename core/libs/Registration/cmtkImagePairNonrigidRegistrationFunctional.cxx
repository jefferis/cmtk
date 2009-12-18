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

#include <cmtkImagePairNonrigidRegistrationFunctional.h>
#include <cmtkImagePairNonrigidRegistrationFunctionalTemplate.h>

#include <cmtkSplineWarpXform.h>

#include <cmtkImagePairSimilarityMeasureCR.h>
#include <cmtkImagePairSimilarityMeasureMSD.h>
#include <cmtkImagePairSimilarityMeasureNCC.h>
#include <cmtkImagePairSimilarityMeasureNMI.h>
#include <cmtkImagePairSimilarityMeasureMI.h>

#include <cmtkInterpolator.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairNonrigidRegistrationFunctional::ImagePairNonrigidRegistrationFunctional
( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
  : ImagePairRegistrationFunctional( reference, floating )
{
  Dim = 0;

  ReferenceFrom.Set( 0,0,0 );
  ReferenceTo.Set( reference->Size );

  this->m_AdaptiveFixParameters = false;
  this->m_AdaptiveFixThreshFactor = 0.5;
  StepScaleVector = NULL;

  VectorCache = Memory::AllocateArray<Vector3D>( ReferenceDims[0] );
  VolumeOfInfluence = NULL;
}

ImagePairNonrigidRegistrationFunctional::~ImagePairNonrigidRegistrationFunctional()
{
  if ( VectorCache ) delete[] VectorCache;
  if ( StepScaleVector ) delete[] StepScaleVector;
}

template<class VM,class W>
void
ImagePairNonrigidRegistrationFunctionalTemplate<VM,W>::WeightedDerivative
( double& lower, double& upper, typename W::SmartPtr& warp, 
  const int param, const Types::Coordinate step ) const
{
  if ( this->m_JacobianConstraintWeight > 0 )
    {
    double lowerConstraint = 0, upperConstraint = 0;
    warp->GetJacobianConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step );
    lower -= this->m_JacobianConstraintWeight * lowerConstraint;
    upper -= this->m_JacobianConstraintWeight * upperConstraint;
    } 

  if ( this->m_RigidityConstraintWeight > 0 )
    {
    double lowerConstraint = 0, upperConstraint = 0;

    if ( this->m_RigidityConstraintMap )
      {
      warp->GetRigidityConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step, this->m_RigidityConstraintMap );
      }
    else
      {
      warp->GetRigidityConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step );
      }
    lower -= this->m_RigidityConstraintWeight * lowerConstraint;
    upper -= this->m_RigidityConstraintWeight * upperConstraint;
    } 
  
  if ( this->m_GridEnergyWeight > 0 ) 
    {
    double lowerEnergy = 0, upperEnergy = 0;
    warp->GetGridEnergyDerivative( lowerEnergy, upperEnergy, param, step );
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
      warp->GetDerivativeLandmarksMSD( lowerMSD, upperMSD, this->m_MatchedLandmarkList, param, step );
      lower -= this->m_LandmarkErrorWeight * lowerMSD;
      upper -= this->m_LandmarkErrorWeight * upperMSD;
      }
    if ( InverseTransformation ) 
      {
      double lowerIC, upperIC;
      warp->GetDerivativeInverseConsistencyError( lowerIC, upperIC, this->InverseTransformation, this->ReferenceGrid, &(this->VolumeOfInfluence[param]), param, step );
      lower -= InverseConsistencyWeight * lowerIC;
      upper -= InverseConsistencyWeight * upperIC;
      }
    }
}

template<class VM,class W> 
void
ImagePairNonrigidRegistrationFunctionalTemplate<VM,W>::SetWarpXform
( WarpXform::SmartPtr& warp )
{
  this->Warp = W::SmartPtr::DynamicCastFrom( warp );
  if ( this->Warp )
    {
    Warp->RegisterVolume( ReferenceGrid );
    Warp->SetIncompressibilityMap( this->m_IncompressibilityMap );
    
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

    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread ) 
      {
      if ( thread ) 
	{
	this->m_ThreadWarp[thread] = typename W::SmartPtr( dynamic_cast<W*>( this->Warp->Clone() ) );
	this->m_ThreadWarp[thread]->RegisterVolume( this->ReferenceGrid );
	} 
      else 
	{
	this->m_ThreadWarp[thread] = W::SmartPtr::DynamicCastFrom( this->Warp );
	}
      } 
    }
}

template<class VM, class W>
void
ImagePairNonrigidRegistrationFunctionalTemplate<VM,W>::UpdateWarpFixedParameters() 
{
  if ( this->m_ConsistencyHistogram.IsNull() ) 
    {
    this->m_ConsistencyHistogram = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>() );
    unsigned int numSamplesX = this->m_Metric->GetNumberOfSamplesX();
    Types::DataItem fromX, toX;
    this->m_Metric->GetDataRangeX( fromX, toX );
    unsigned int numBinsX = this->m_ConsistencyHistogram->CalcNumBins( numSamplesX, fromX, toX );
    
    unsigned int numSamplesY = this->m_Metric->GetNumberOfSamplesY();
    Types::DataItem fromY, toY;
    this->m_Metric->GetDataRangeY( fromY, toY );
    unsigned int numBinsY = this->m_ConsistencyHistogram->CalcNumBins( numSamplesY, fromY, toY );
    
    this->m_ConsistencyHistogram->SetNumBins( numBinsX, numBinsY );
    this->m_ConsistencyHistogram->SetRangeX( fromX, toX );
    this->m_ConsistencyHistogram->SetRangeY( fromY, toY );
    }
  
  int numCtrlPoints = this->Dim / 3;
  
  double *mapRef = Memory::AllocateArray<double>( numCtrlPoints );
  double *mapMod = Memory::AllocateArray<double>( numCtrlPoints );

  Rect3D voi;
  Vector3D fromVOI, toVOI;
  int pX, pY, pZ;

  int inactive = 0;

  const Types::DataItem unsetY = DataTypeTraits<Types::DataItem>::ChoosePaddingValue();

  if ( this->ReferenceDataClass == DATACLASS_LABEL ) 
    {
    this->Warp->SetParameterActive();
    
    for ( int ctrl = 0; ctrl < numCtrlPoints; ++ctrl ) 
      {
      /// We cannot use the precomputed table of VOIs here because in "fast"
      /// mode, these VOIs are smaller than we want them here.
      this->Warp->GetVolumeOfInfluence( 3 * ctrl, this->ReferenceFrom, this->ReferenceTo, fromVOI, toVOI, 0 );
      this->GetReferenceGridRange( fromVOI, toVOI, voi );
      
      int r = voi.startX + this->m_DimsX * ( voi.startY + this->m_DimsY * voi.startZ );
      
      bool active = false;
      for ( pZ = voi.startZ; (pZ < voi.endZ) && !active; ++pZ ) 
	{
	for ( pY = voi.startY; (pY < voi.endY) && !active; ++pY ) 
	  {
	  for ( pX = voi.startX; (pX < voi.endX); ++pX, ++r ) 
	    {
	    if ( ( this->m_Metric->GetSampleX( r ) != 0 ) || ( ( this->m_WarpedVolume[r] != unsetY ) && ( this->m_WarpedVolume[r] != 0 ) ) ) 
	      {
	      active = true;
	      break;
	      }
	    }
	  r += ( voi.startX + ( this->m_DimsX-voi.endX ) );
	  }
	r += this->m_DimsX * ( voi.startY + ( this->m_DimsY-voi.endY ) );
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
      this->m_ConsistencyHistogram->Reset();
      
      // We cannot use the precomputed table of VOIs here because in "fast"
      // mode, these VOIs are smaller than we want them here.
      this->Warp->GetVolumeOfInfluence( 3 * ctrl, this->ReferenceFrom, this->ReferenceTo, fromVOI, toVOI, 0 );
      this->GetReferenceGridRange( fromVOI, toVOI, voi );
      
      int r = voi.startX + this->m_DimsX * ( voi.startY + this->m_DimsY * voi.startZ );
      
      const int endOfLine = ( voi.startX + ( this->m_DimsX-voi.endX) );
      const int endOfPlane = this->m_DimsX * ( voi.startY + (this->m_DimsY-voi.endY) );
      
      for ( pZ = voi.startZ; pZ<voi.endZ; ++pZ ) 
	{
	for ( pY = voi.startY; pY<voi.endY; ++pY ) 
	  {
	  for ( pX = voi.startX; pX<voi.endX; ++pX, ++r ) 
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

ImagePairNonrigidRegistrationFunctional* 
ImagePairNonrigidRegistrationFunctional::Create
( const int metric,
  UniformVolume::SmartPtr& refVolume, 
  UniformVolume::SmartPtr& fltVolume,
  const Interpolators::InterpolationEnum interpolation )
{
  switch ( metric ) 
    {
    case 0:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNMI,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 1:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMI,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 2:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureCR,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 3:
      return NULL; // masked NMI retired
    case 4:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMSD,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 5:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNCC,SplineWarpXform>( refVolume, fltVolume, interpolation );
    default:
      return NULL;
    }

  return NULL;
}

} // namespace cmtk
