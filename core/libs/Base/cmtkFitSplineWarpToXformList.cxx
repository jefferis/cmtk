/*
//
//  Copyright 2012, 2014 SRI International
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

#include "cmtkFitSplineWarpToXformList.h"

#include <Base/cmtkRegionIndexIterator.h>

#include <System/cmtkDebugOutput.h>

void cmtk::FitSplineWarpToXformList::ComputeResiduals(
    const SplineWarpXform &splineWarp) {
  const DataGrid::IndexType dims = this->m_XformField.m_Dims;

  this->m_Residuals.resize(dims.Product());

#pragma omp parallel for
  for (int z = 0; z < dims[2]; ++z) {
    size_t ofs = z * dims[0] * dims[1];

    for (int y = 0; y < dims[1]; ++y) {
      for (int x = 0; x < dims[0]; ++x, ++ofs) {
        if (this->m_XformValidAt[ofs])
          this->m_Residuals[ofs] =
              this->m_XformField[ofs] - splineWarp.GetTransformedGrid(x, y, z);
      }
    }
  }
}

cmtk::SplineWarpXform::SmartPtr cmtk::FitSplineWarpToXformList::Fit(
    const SplineWarpXform::ControlPointIndexType &finalDims, const int nLevels,
    const bool fitAffineFirst) {
  AffineXform::SmartPtr initialAffine;
  if (fitAffineFirst) {
    initialAffine = this->Superclass::Fit();
  } else {
    initialAffine = AffineXform::SmartPtr(new AffineXform);
  }

  // we may have to adjust nLevels downwards
  int numberOfLevels = nLevels;

  SplineWarpXform::ControlPointIndexType initialDims = finalDims;
  for (int level = 1; level < numberOfLevels; ++level) {
    if ((initialDims[0] & 1) && (initialDims[1] & 1) &&
        (initialDims[2] & 1) &&  // check that all dims are still odd numbers
        (initialDims.MinValue() > 4))  // check that at least 4 control points
                                       // will be left in each dimension
    {
      // apply the inverse of the refinement formula used in
      // SplineWarpXform::Refine()
      initialDims.AddScalar(+3);
      initialDims /= 2;
    } else {
      numberOfLevels = level;

      DebugOutput(2) << "INFO: adjusted number of levels to " << numberOfLevels
                     << " from " << nLevels
                     << " to ensure sufficient number of control points\n";
    }
  }

  SplineWarpXform *splineWarp =
      new SplineWarpXform(this->m_XformField.m_Size, initialDims,
                          CoordinateVector::SmartPtr::Null(), initialAffine);
  this->FitSpline(*splineWarp, numberOfLevels);

  return cmtk::SplineWarpXform::SmartPtr(splineWarp);
}

cmtk::SplineWarpXform::SmartPtr cmtk::FitSplineWarpToXformList::Fit(
    const Types::Coordinate finalSpacing, const int nLevels,
    const bool fitAffineFirst) {
  AffineXform::SmartPtr initialAffine;
  if (fitAffineFirst) {
    initialAffine = this->Superclass::Fit();
  } else {
    initialAffine = AffineXform::SmartPtr(new AffineXform);
  }

  // compute the start spacing of multi-level approximation by doubling final
  // spacing until user-defined initial spacing is exceeded.
  Types::Coordinate spacing = finalSpacing * (1 << (nLevels - 1));

  // initialize B-spline transformation
  SplineWarpXform *splineWarp =
      new SplineWarpXform(this->m_XformField.m_Size, spacing, initialAffine);
  this->FitSpline(*splineWarp, nLevels);

  return cmtk::SplineWarpXform::SmartPtr(splineWarp);
}

void cmtk::FitSplineWarpToXformList::FitSpline(SplineWarpXform &splineWarp,
                                               const int nLevels) {
  // loop until final control point spacing
  for (int level = 0; level < nLevels; ++level) {
    DebugOutput(5) << "Multi-resolution spline fitting level " << level + 1
                   << " out of " << nLevels << "\n";
    // refine control point grid unless this is first iteration
    if (level) {
      splineWarp.Refine();
    }

    DebugOutput(6) << "  Control point grid is " << splineWarp.m_Dims[0] << "x"
                   << splineWarp.m_Dims[1] << "x" << splineWarp.m_Dims[2]
                   << "\n";

    // compute residuals
    splineWarp.RegisterVolume(this->m_XformField);
    this->ComputeResiduals(splineWarp);

    // loop over all control points to compute deltas as the spline coefficients
    // that fit current residuals
    std::vector<FixedVector<3, Types::Coordinate>> delta(
        splineWarp.m_NumberOfControlPoints,
        FixedVector<3, Types::Coordinate>(0.0));
    std::vector<Types::Coordinate> weight(splineWarp.m_NumberOfControlPoints,
                                          0.0);

    for (RegionIndexIterator<DataGrid::RegionType> voxelIt(
             this->m_XformField.GetWholeImageRegion());
         voxelIt != voxelIt.end(); ++voxelIt) {
      const DataGrid::IndexType voxelIdx = voxelIt.Index();

      const size_t voxelOfs = this->m_XformField.GetOffsetFromIndex(voxelIdx);
      if (this->m_XformValidAt[voxelOfs]) {
        Types::Coordinate sumOfSquares = 0;
        Types::Coordinate wklm[4][4][4], w2klm[4][4][4];
        for (int m = 0; m < 4; ++m) {
          for (int l = 0; l < 4; ++l) {
            const Types::Coordinate wlm =
                splineWarp.m_GridSpline[1][4 * voxelIdx[1] + l] *
                splineWarp.m_GridSpline[2][4 * voxelIdx[2] + m];
            for (int k = 0; k < 4; ++k) {
              sumOfSquares +=
                  (w2klm[m][l][k] = MathUtil::Square(
                       wklm[m][l][k] =
                           splineWarp.m_GridSpline[0][4 * voxelIdx[0] + k] *
                           wlm));
            }
          }
        }

        for (int m = 0; m < 4; ++m) {
          const size_t mOfs = splineWarp.m_Dims[1] *
                              (splineWarp.m_GridIndexes[2][voxelIdx[2]] + m);
          for (int l = 0; l < 4; ++l) {
            const size_t mlOfs =
                splineWarp.m_Dims[0] *
                (splineWarp.m_GridIndexes[1][voxelIdx[1]] + l + mOfs);
            for (int k = 0; k < 4; ++k) {
              const size_t cpOfs =
                  splineWarp.m_GridIndexes[0][voxelIdx[0]] + k + mlOfs;
              delta[cpOfs] += w2klm[m][l][k] * wklm[m][l][k] / sumOfSquares *
                              this->m_Residuals[voxelOfs];
              weight[cpOfs] += w2klm[m][l][k];
            }
          }
        }
      }
    }

    // apply delta
#pragma omp parallel for
    for (int cp = 0; cp < static_cast<int>(splineWarp.m_NumberOfControlPoints);
         ++cp) {
      if (weight[cp] != 0) {
        delta[cp] /= weight[cp];
        splineWarp.SetShiftedControlPointPositionByOffset(
            splineWarp.GetShiftedControlPointPositionByOffset(cp) + delta[cp],
            cp);
      } else {
        // nothing to do - keep control point where it is.
      }
    }
  }
}
