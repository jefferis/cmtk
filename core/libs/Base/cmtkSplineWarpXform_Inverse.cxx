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

#include "cmtkSplineWarpXform.h"

#include <Base/cmtkMathUtil.h>

namespace cmtk {

/** \addtogroup Base */
//@{

bool SplineWarpXform::ApplyInverse(const Self::SpaceVectorType &v,
                                   Self::SpaceVectorType &u,
                                   const Types::Coordinate accuracy) const {
  return this->ApplyInverseWithInitial(v, u, this->FindClosestControlPoint(v),
                                       accuracy);
}

SplineWarpXform::SpaceVectorType SplineWarpXform::FindClosestControlPoint(
    const Self::SpaceVectorType &v) const {
  // find closest control point -- we'll go from there.
  Types::Coordinate closestDistance = FLT_MAX;
  Types::Coordinate idx[3];
  for (int dim = 0; dim < 3; ++dim) idx[dim] = 0.5 * this->m_Dims[dim];

  for (Types::Coordinate step = 0.25 * MathUtil::Min(3, idx); step > 0.01;
       step *= 0.5) {
    bool improved = true;
    while (improved) {
      improved = false;
      int closestDim = 0, closestDir = 0;

      for (int dim = 0; dim < 3; ++dim) {
        for (int dir = -1; dir < 2; dir += 2) {
          const Types::Coordinate oldIdx = idx[dim];
          idx[dim] += dir * step;
          if ((idx[dim] > 0) && (idx[dim] <= this->m_Dims[dim] - 2)) {
            Self::SpaceVectorType cp = this->Apply(
                this->GetOriginalControlPointPosition(idx[0], idx[1], idx[2]));
            cp -= v;
            const Types::Coordinate distance = cp.RootSumOfSquares();
            if (distance < closestDistance) {
              closestDistance = distance;
              closestDim = dim;
              closestDir = dir;
              improved = true;
            }
          }
          idx[dim] = oldIdx;
        }
      }

      if (improved) {
        idx[closestDim] += closestDir * step;
      }
    }
  }

  assert((idx[0] <= this->m_Dims[0] - 1) && (idx[1] <= this->m_Dims[1] - 1) &&
         (idx[2] <= this->m_Dims[2] - 1));
  assert(idx[0] >= 0 && idx[1] >= 0 && idx[2] >= 0);

  return this->GetOriginalControlPointPosition(idx[0], idx[1], idx[2]);
}

}  // namespace cmtk
