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

#include <cmtkSymmetricElasticFunctional.h>

#include <cmtkVoxelMatchingMutInf.h>
#include <cmtkVoxelMatchingNormMutInf.h>
#include <cmtkVoxelMatchingCorrRatio.h>
#include <cmtkVoxelMatchingMeanSquaredDifference.h>
#include <cmtkVoxelMatchingCrossCorrelation.h>

#include <cmtkSplineWarpXform.h>

#include <cmtkInterpolator.h>

#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class VM, class W>
void
SymmetricElasticFunctional_Template<VM,W>::SetWarpXform
( WarpXform::SmartPtr& warpFwd, WarpXform::SmartPtr& warpBwd ) 
{
  this->FwdFunctional.SetWarpXform( warpFwd );
  this->FwdFunctional.SetInverseTransformation( warpBwd );
  
  this->BwdFunctional.SetWarpXform( warpBwd );
  this->BwdFunctional.SetInverseTransformation( warpFwd );
}

template<class VM, class W>
typename SymmetricElasticFunctional_Template<VM,W>::ReturnType
SymmetricElasticFunctional_Template<VM,W>
::EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
  CoordinateVector vFwd( this->FwdFunctional.ParamVectorDim(), v.Elements, false /*freeElements*/ );
  CoordinateVector gFwd( this->FwdFunctional.ParamVectorDim(), g.Elements, false /*freeElements*/ );

  CoordinateVector vBwd( this->BwdFunctional.ParamVectorDim(), v.Elements+this->FwdFunctional.ParamVectorDim(), false /*freeElements*/ );
  CoordinateVector gBwd( this->BwdFunctional.ParamVectorDim(), g.Elements+this->FwdFunctional.ParamVectorDim(), false /*freeElements*/ );

  const typename Self::ReturnType result = this->FwdFunctional.EvaluateWithGradient( vFwd, gFwd, step ) + this->BwdFunctional.EvaluateWithGradient( vBwd, gBwd, step );
  return result;
}

SymmetricElasticFunctional* 
CreateSymmetricElasticFunctional( const int metric, UniformVolume::SmartPtr& refVolume,  UniformVolume::SmartPtr& fltVolume )
{
  switch ( fltVolume->GetData()->GetDataClass() ) 
    {
    case DATACLASS_UNKNOWN :
    case DATACLASS_GREY :
      switch ( metric ) 
	{
	case 0:
	  return new SymmetricElasticFunctional_Template< VoxelMatchingNormMutInf_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 1:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingMutInf_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 2:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingCorrRatio_Trilinear,SplineWarpXform>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked NMI retired
	case 4:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingMeanSquaredDifference,SplineWarpXform>( refVolume, fltVolume );
	case 5:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingCrossCorrelation,SplineWarpXform>( refVolume, fltVolume );
	default:
	  return NULL;
	}
    case DATACLASS_LABEL:
      switch ( metric ) 
	{
	case 0:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingNormMutInf_NearestNeighbor, SplineWarpXform>( refVolume, fltVolume );
	case 1:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingMutInf_NearestNeighbor,SplineWarpXform>( refVolume, fltVolume );
	case 2:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingCorrRatio_NearestNeighbor,SplineWarpXform>( refVolume, fltVolume );
	case 3:
	  return NULL; // masked NMI retired
	case 4:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingMeanSquaredDifference,SplineWarpXform>( refVolume, fltVolume );
	case 5:
	  return new SymmetricElasticFunctional_Template<VoxelMatchingCrossCorrelation,SplineWarpXform>( refVolume, fltVolume );
	default:
	  return NULL;
	}
    }
  
  return NULL;
}

} // namespace cmtk
