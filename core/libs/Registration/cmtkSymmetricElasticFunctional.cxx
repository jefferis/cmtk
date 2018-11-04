/*
//
//  Copyright 2016 Google, Inc.
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

#include <Registration/cmtkSymmetricElasticFunctional.h>

#include <Registration/cmtkVoxelMatchingCorrRatio.h>
#include <Registration/cmtkVoxelMatchingCrossCorrelation.h>
#include <Registration/cmtkVoxelMatchingMeanSquaredDifference.h>
#include <Registration/cmtkVoxelMatchingMutInf.h>
#include <Registration/cmtkVoxelMatchingNormMutInf.h>

#include <Base/cmtkInterpolator.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkUniformVolume.h>

namespace cmtk {

/** \addtogroup Registration */
//@{

template <class VM>
void SymmetricElasticFunctional_Template<VM>::SetWarpXform(
    SplineWarpXform::SmartPtr &warpFwd, SplineWarpXform::SmartPtr &warpBwd) {
  this->FwdFunctional.SetWarpXform(warpFwd);
  this->FwdFunctional.SetInverseTransformation(warpBwd);

  this->BwdFunctional.SetWarpXform(warpBwd);
  this->BwdFunctional.SetInverseTransformation(warpFwd);
}

template <class VM>
typename SymmetricElasticFunctional_Template<VM>::ReturnType
SymmetricElasticFunctional_Template<VM>::EvaluateWithGradient(
    CoordinateVector &v, CoordinateVector &g, const Types::Coordinate step) {
  CoordinateVector vFwd(this->FwdFunctional.ParamVectorDim(), v.Elements,
                        false /*freeElements*/);
  CoordinateVector gFwd(this->FwdFunctional.ParamVectorDim(), g.Elements,
                        false /*freeElements*/);

  CoordinateVector vBwd(this->BwdFunctional.ParamVectorDim(),
                        v.Elements + this->FwdFunctional.ParamVectorDim(),
                        false /*freeElements*/);
  CoordinateVector gBwd(this->BwdFunctional.ParamVectorDim(),
                        g.Elements + this->FwdFunctional.ParamVectorDim(),
                        false /*freeElements*/);

  return this->FwdFunctional.EvaluateWithGradient(vFwd, gFwd, step) +
         this->BwdFunctional.EvaluateWithGradient(vBwd, gBwd, step);
}

SymmetricElasticFunctional *CreateSymmetricElasticFunctional(
    const int metric, UniformVolume::SmartPtr &refVolume,
    UniformVolume::SmartPtr &fltVolume) {
  switch (fltVolume->GetData()->GetDataClass()) {
    case DATACLASS_UNKNOWN:
    case DATACLASS_GREY:
      switch (metric) {
        case 0:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingNormMutInf_Trilinear>(refVolume, fltVolume);
        case 1:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingMutInf_Trilinear>(refVolume, fltVolume);
        case 2:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingCorrRatio_Trilinear>(refVolume, fltVolume);
        case 3:
          return NULL;  // masked NMI retired
        case 4:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingMeanSquaredDifference>(refVolume, fltVolume);
        case 5:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingCrossCorrelation>(refVolume, fltVolume);
        default:
          return NULL;
      }
    case DATACLASS_LABEL:
      switch (metric) {
        case 0:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingNormMutInf_NearestNeighbor>(refVolume, fltVolume);
        case 1:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingMutInf_NearestNeighbor>(refVolume, fltVolume);
        case 2:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingCorrRatio_NearestNeighbor>(refVolume, fltVolume);
        case 3:
          return NULL;  // masked NMI retired
        case 4:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingMeanSquaredDifference>(refVolume, fltVolume);
        case 5:
          return new SymmetricElasticFunctional_Template<
              VoxelMatchingCrossCorrelation>(refVolume, fltVolume);
        default:
          return NULL;
      }
  }

  return NULL;
}

}  // namespace cmtk
