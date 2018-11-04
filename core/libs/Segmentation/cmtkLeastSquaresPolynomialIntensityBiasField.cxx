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

#include "cmtkLeastSquaresPolynomialIntensityBiasField.h"

#include <Base/cmtkLeastSquares.h>
#include <Base/cmtkMatrix.h>
#include <Base/cmtkPolynomial.h>
#include <Base/cmtkRegionIndexIterator.h>

namespace cmtk {

LeastSquaresPolynomialIntensityBiasField::
    LeastSquaresPolynomialIntensityBiasField(const UniformVolume &image,
                                             const std::vector<bool> &mask,
                                             const int degree) {
  const UniformVolume::CoordinateVectorType center =
      image.GetCenterCropRegion();

  // first, compute average intensity over masked region
  Types::DataItem avg = 0;
  size_t nPixelsMask = 0;

  const DataGrid::RegionType region = image.GetWholeImageRegion();
  for (RegionIndexIterator<DataGrid::RegionType> it(region); it != it.end();
       ++it) {
    const size_t ofs = image.GetOffsetFromIndex(it.Index());
    if (mask[ofs]) {
      avg += fabs(image.GetDataAt(ofs));
      ++nPixelsMask;
    }
  }

  if (!nPixelsMask) throw Self::EmptyMaskException();

  avg /= nPixelsMask;

  // set up least-squares problem
  const size_t nVars = PolynomialHelper::GetNumberOfMonomials(degree);
  if (nVars <
      2) {  // we ignore constant term, so we need at least 2 variables to be
            // able to do anything meningful. for fewer, just use original data
    this->m_CorrectedData = image.GetData();
    return;
  }

  std::vector<Types::DataItem> dataVector(nPixelsMask);
  Matrix2D<Types::DataItem> uMatrix(
      nPixelsMask, nVars - 1);  // nVars-1 because we ignore zero-order term

  size_t cntPx = 0;
  for (RegionIndexIterator<DataGrid::RegionType> it(region); it != it.end();
       ++it) {
    const size_t ofs = image.GetOffsetFromIndex(it.Index());

    if (mask[ofs]) {
      const UniformVolume::CoordinateVectorType xyz = ComponentDivide(
          image.GetGridLocation(it.Index()) - center, image.m_Size);
      dataVector[cntPx] = image.GetDataAt(ofs) / avg - 1.0;
      for (size_t n = 1; n < nVars; ++n) {
        uMatrix[cntPx][n - 1] =
            Polynomial<4, Types::DataItem>::EvaluateMonomialAt(n, xyz[0],
                                                               xyz[1], xyz[2]);
      }
      ++cntPx;
    }
  }

  // solve least-squares problem
  const std::vector<Types::DataItem> params =
      LeastSquares<Types::DataItem>(uMatrix).Solve(dataVector);

  // apply solution
  this->m_CorrectedData =
      TypedArray::Create(image.GetData()->GetType(), image.GetNumberOfPixels());
  this->m_BiasData = TypedArray::Create(TYPE_ITEM, image.GetNumberOfPixels());

  for (RegionIndexIterator<DataGrid::RegionType> it(region); it != it.end();
       ++it) {
    const size_t ofs = image.GetOffsetFromIndex(it.Index());

    const UniformVolume::CoordinateVectorType xyz = ComponentDivide(
        image.GetGridLocation(it.Index()) - center, image.m_Size);

    Types::DataItem bias = 1.0;
    for (size_t n = 1; n < nVars; ++n) {
      bias +=
          params[n - 1] * Polynomial<4, Types::DataItem>::EvaluateMonomialAt(
                              n, xyz[0], xyz[1], xyz[2]);
    }

    this->m_BiasData->Set(bias, ofs);

    Types::DataItem value;
    if (image.GetData()->Get(value, ofs)) {
      this->m_CorrectedData->Set(value / bias, ofs);
    } else {
      this->m_CorrectedData->SetPaddingAt(ofs);
    }
  }
}

}  // namespace cmtk
