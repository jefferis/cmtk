/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#include "cmtkFilterVolume.h"

#include <Base/cmtkFilterMask.h>
#include <Base/cmtkHistogram.h>
#include <Base/cmtkJointHistogram.h>
#include <Base/cmtkMathUtil.h>

#include <System/cmtkProgress.h>
#include <System/cmtkThreads.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cmtk {

/** \addtogroup Base */
//@{

TypedArray::SmartPtr FilterVolume::GaussianFilter(
    const UniformVolume *volume, const Units::GaussianSigma &kernelWidth,
    const Types::Coordinate radius, const TypedArray *maskData) {
  const TypedArray *inputData = volume->GetData();
  if (!inputData) return TypedArray::SmartPtr(NULL);

  TypedArray::SmartPtr filtered =
      TypedArray::Create(inputData->GetType(), inputData->GetDataSize());

  const DataGrid::IndexType &dims = volume->m_Dims;
  FilterMask<3> filter(dims, volume->Deltas(), radius,
                       FilterMask<3>::Gaussian(kernelWidth));

  const Types::GridIndexType dimsX = dims[AXIS_X];
  const Types::GridIndexType dimsY = dims[AXIS_Y];
  const Types::GridIndexType dimsZ = dims[AXIS_Z];

  Progress::Begin(0, dimsZ, 1, "Gaussian Filter");

#pragma omp parallel for
  for (Types::GridIndexType z = 0; z < dimsZ; ++z) {
    size_t offset = z * dimsX * dimsY;
    Progress::SetProgress(z);
    for (Types::GridIndexType y = 0; y < dimsY; ++y)
      for (Types::GridIndexType x = 0; x < dimsX; ++x, ++offset) {
        Types::DataItem average = 0.0, weight = 0.0;

        Types::DataItem maskValue = 0.0;
        if (maskData) {
          maskData->Get(maskValue, offset);
        } else {
          maskValue = 1.0;
        }

        if (maskValue) {
          FilterMask<3>::const_iterator it = filter.begin();
          while (it != filter.end()) {
            const Types::GridIndexType xx = x + it->Location[0];
            const Types::GridIndexType yy = y + it->Location[1];
            const Types::GridIndexType zz = z + it->Location[2];

            if ((xx >= 0) && (yy >= 0) && (zz >= 0) &&
                (xx < (Types::GridIndexType)dimsX) &&
                (yy < (Types::GridIndexType)dimsY) &&
                (zz < (Types::GridIndexType)dimsZ)) {
              const size_t srcOffset = volume->GetOffsetFromIndex(xx, yy, zz);
              Types::DataItem value;
              if (inputData->Get(value, srcOffset)) {
                average += it->Coefficient * value;
                weight += it->Coefficient;
              }
            }
            ++it;
          }
        }
        if (weight > 0.0) {
          filtered->Set(average / weight, offset);
        } else {
          filtered->SetPaddingAt(offset);
        }
      }
  }

  Progress::Done();

  return filtered;
}

TypedArray::SmartPtr FilterVolume ::RohlfingFilter(
    const UniformVolume *volume, const TypedArray *subjectData,
    const TypedArray *maskData, const Units::GaussianSigma &iFilterSigma,
    const Units::GaussianSigma &filterWidth,
    const Types::Coordinate filterRadius) {
  const TypedArray *inputData = volume->GetData();
  if (!inputData) return TypedArray::SmartPtr(NULL);

  const Types::DataItemRange rangeSubj = subjectData->GetRange();

  const size_t numBins = 1024;
#ifdef _OPENMP
  const size_t maxThreads = omp_get_max_threads();
  std::vector<Histogram<Types::DataItem>::SmartPtr> histograms(maxThreads);
  for (size_t thread = 0; thread < maxThreads; ++thread) {
    histograms[thread] = Histogram<Types::DataItem>::SmartPtr(
        new Histogram<Types::DataItem>(numBins));
    histograms[thread]->SetRange(rangeSubj);
  }
#else   // #ifdef _OPENMP
  Histogram<Types::DataItem> histogram(numBins);
  histogram.SetRange(rangeSubj);
#endif  // #ifdef _OPENMP

  const size_t iKernelRadius =
      1 + static_cast<size_t>(2 * iFilterSigma.Value() * numBins);
  std::vector<Types::DataItem> iKernel(iKernelRadius);
  if (iKernelRadius > 1) {
    const Types::DataItem normFactor = static_cast<Types::DataItem>(
        1.0 / (sqrt(2 * M_PI) * iFilterSigma.Value() *
               numBins));  // not really necessary since we normalize during
                           // convolution
    for (size_t i = 0; i < iKernelRadius; ++i) {
      iKernel[i] = static_cast<Types::DataItem>(
          normFactor *
          exp(-MathUtil::Square(1.0 * i / (iFilterSigma.Value() * numBins)) /
              2));
    }
  } else {
    iKernel[0] = 1.0;
  }

  TypedArray::SmartPtr filtered =
      TypedArray::Create(inputData->GetType(), inputData->GetDataSize());

  const DataGrid::IndexType &dims = volume->GetDims();
  FilterMask<3> filter(dims, volume->Deltas(), filterRadius,
                       FilterMask<3>::Gaussian(filterWidth));

  const Types::GridIndexType dimsX = dims[AXIS_X];
  const Types::GridIndexType dimsY = dims[AXIS_Y];
  const Types::GridIndexType dimsZ = dims[AXIS_Z];

  Progress::Begin(0, dimsZ, 1, "Rohlfing Intensity-Consistent Filter");

#pragma omp parallel for
  for (Types::GridIndexType z = 0; z < static_cast<Types::GridIndexType>(dimsZ);
       ++z) {
    size_t offset = z * dimsX * dimsY;

#ifdef _OPENMP
    const size_t threadIdx = omp_get_thread_num();
    Histogram<Types::DataItem> &histogram = *(histograms[threadIdx]);
    if (!threadIdx)
#endif  // #ifdef _OPENMP
      Progress::SetProgress(z);

    for (Types::GridIndexType y = 0; y < dimsY; ++y)
      for (Types::GridIndexType x = 0; x < dimsX; ++x, ++offset) {
        Types::DataItem average = 0.0, weight = 0.0;

        Types::DataItem maskValue = 1.0;
        if (maskData) maskData->Get(maskValue, offset);

        Types::DataItem valueSubjCenter;
        if (maskValue && subjectData->Get(valueSubjCenter, offset)) {
          histogram.Reset();
          histogram.AddWeightedSymmetricKernel(
              histogram.ValueToBin(valueSubjCenter), iKernelRadius,
              &(iKernel[0]));

          for (FilterMask<3>::const_iterator it = filter.begin();
               it != filter.end(); ++it) {
            const Types::GridIndexType xx = x + it->Location[0];
            const Types::GridIndexType yy = y + it->Location[1];
            const Types::GridIndexType zz = z + it->Location[2];

            if ((xx >= 0) && (yy >= 0) && (zz >= 0) &&
                (xx < (Types::GridIndexType)dimsX) &&
                (yy < (Types::GridIndexType)dimsY) &&
                (zz < (Types::GridIndexType)dimsZ)) {
              Types::DataItem value;
              const size_t srcOffset = it->RelativeIndex + offset;
              if (inputData->Get(value, srcOffset)) {
                Types::DataItem valueSubj;
                if (subjectData->Get(valueSubj, srcOffset)) {
                  const size_t bin = histogram.ValueToBin(valueSubj);
                  const Types::DataItem prob = it->Coefficient * histogram[bin];

                  average += value * prob;
                  weight += prob;
                }
              }
            }
          }
        }

        if (weight > 0.0) {
          filtered->Set(average / weight, offset);
        } else {
          filtered->SetPaddingAt(offset);
        }
      }
  }

  Progress::Done();

  return filtered;
}

TypedArray::SmartPtr FilterVolume::StudholmeFilter(
    const UniformVolume *volume, const TypedArray *subjectData,
    const TypedArray *averageData, const TypedArray *maskData,
    std::list<TypedArray::SmartPtr> imgList, const Types::DataItem binWidth,
    const Units::GaussianSigma &filterWidth,
    const Types::Coordinate filterRadius) {
  const TypedArray *inputData = volume->GetData();
  if (!inputData) return TypedArray::SmartPtr(NULL);

  const Types::DataItemRange range = averageData->GetRange();
  const size_t numBins =
      std::min(128, 1 + static_cast<int>((range.Width()) / binWidth));

  TypedArray::SmartPtr filtered =
      TypedArray::Create(inputData->GetType(), inputData->GetDataSize());

  const DataGrid::IndexType &dims = volume->GetDims();
  const Types::GridIndexType dimsX = dims[AXIS_X];
  const Types::GridIndexType dimsY = dims[AXIS_Y];
  const Types::GridIndexType dimsZ = dims[AXIS_Z];
  const Types::GridIndexType numberOfRows = dimsY * dimsZ;

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  std::vector<JointHistogram<Types::DataItem>> histogramByThread(
      numberOfThreads);
  std::vector<FilterMask<3>::SmartPtr> filterByThread(numberOfThreads);

  for (size_t idx = 0; idx < numberOfThreads; ++idx) {
    histogramByThread[idx].Resize(numBins, numBins);
    histogramByThread[idx].SetRangeX(range);
    histogramByThread[idx].SetRangeY(range);

    FilterMask<3>::SmartPtr filter(
        new FilterMask<3>(dims, volume->Deltas(), filterRadius,
                          FilterMask<3>::Gaussian(filterWidth)));
    filterByThread[idx] = filter;
  }

  Progress::Begin(0, numberOfRows, 1, "Studholme Intensity-Consistent Filter");
#pragma omp parallel for
  for (Types::GridIndexType row = 0;
       row < static_cast<Types::GridIndexType>(numberOfRows); ++row) {
    const Types::GridIndexType y = row % dimsY;
    const Types::GridIndexType z = row / dimsY;

    Progress::SetProgress(z);
    size_t offset = row * dimsX;

#ifdef _OPENMP
    const int thread = omp_get_thread_num();
#else
    const int thread = 0;
#endif

    JointHistogram<Types::DataItem> &histogram = histogramByThread[thread];
    FilterMask<3> &filter = *(filterByThread[thread]);

    for (Types::GridIndexType x = 0; x < dimsX; ++x, ++offset) {
      Types::DataItem average = 0.0, weight = 0.0;
      histogram.Reset();

      Types::DataItem maskValue = 1.0;
      if (maskData) maskData->Get(maskValue, offset);

      Types::DataItem valueAvg;
      if (maskValue && averageData->Get(valueAvg, offset)) {
        // first iteration over filter: compute consistency histogram
        FilterMask<3>::iterator it = filter.begin();
        for (; it != filter.end(); ++it) {
          const Types::GridIndexType xx = x + it->Location[0];
          const Types::GridIndexType yy = y + it->Location[1];
          const Types::GridIndexType zz = z + it->Location[2];

          if ((xx >= 0) && (yy >= 0) && (zz >= 0) &&
              (xx < (Types::GridIndexType)dimsX) &&
              (yy < (Types::GridIndexType)dimsY) &&
              (zz < (Types::GridIndexType)dimsZ)) {
            it->Valid = true;

            const size_t srcOffset = it->RelativeIndex + offset;
            it->PixelIndex = srcOffset;

            Types::DataItem valueAvgSrc, valueSubj;
            if (averageData->Get(valueAvgSrc, srcOffset)) {
              const size_t binAvg = histogram.ValueToBinX(valueAvgSrc);
              for (std::list<TypedArray::SmartPtr>::iterator itImg =
                       imgList.begin();
                   itImg != imgList.end(); ++itImg) {
                if ((*itImg)->Get(valueSubj, srcOffset)) {
                  histogram.Increment(binAvg, histogram.ValueToBinY(valueSubj));
                }
              }
            }
          } else {
            it->Valid = false;
          }
        }

        const size_t binX = histogram.ValueToBinX(valueAvg);
        const Types::DataItem avgHistogramValueInv =
            static_cast<Types::DataItem>(1.0 / histogram.ProjectToX(binX));

        for (it = filter.begin(); it != filter.end(); ++it) {
          if (it->Valid) {
            Types::DataItem value;
            if (inputData->Get(value, it->PixelIndex)) {
              Types::DataItem valueSubj;
              if (subjectData->Get(valueSubj, it->PixelIndex)) {
                const size_t binY = histogram.ValueToBinY(valueSubj);

                const Types::DataItem prob = static_cast<Types::DataItem>(
                    it->Coefficient * avgHistogramValueInv *
                    histogram.GetBin(binX, binY));

                average += value * prob;
                weight += prob;
              }
            }
          }
        }
      }

      if (weight > 0.0) {
        filtered->Set(average / weight, offset);
      } else {
        filtered->SetPaddingAt(offset);
      }
    }
  }

  Progress::Done();

  return filtered;
}

TypedArray::SmartPtr FilterVolume::StudholmeFilter(
    const UniformVolume *volume, std::list<TypedArray::SmartPtr> subjectData,
    const TypedArray *averageData, const TypedArray *maskData,
    std::list<TypedArray::SmartPtr> imgList, const Types::DataItem binWidth,
    const Units::GaussianSigma &filterWidth,
    const Types::Coordinate filterRadius) {
  const TypedArray *inputData = volume->GetData();
  if (!inputData) return TypedArray::SmartPtr(NULL);

  const Types::DataItemRange range = averageData->GetRange();

  const size_t numBins =
      std::min(128, 1 + static_cast<int>(range.Width() / binWidth));
  JointHistogram<Types::DataItem> histogram(numBins, numBins);
  histogram.SetRangeX(range);
  histogram.SetRangeY(range);

  TypedArray::SmartPtr filtered =
      TypedArray::Create(inputData->GetType(), inputData->GetDataSize());

  const DataGrid::IndexType &dims = volume->GetDims();
  FilterMask<3> filter(dims, volume->Deltas(), filterRadius,
                       FilterMask<3>::Gaussian(filterWidth));

  const Types::GridIndexType dimsX = dims[AXIS_X];
  const Types::GridIndexType dimsY = dims[AXIS_Y];
  const Types::GridIndexType dimsZ = dims[AXIS_Z];

  Progress::Begin(0, dimsZ, 1, "Studholme Intensity-Consistent Filter");

  size_t offset = 0;
  for (Types::GridIndexType z = 0; z < dimsZ; ++z) {
    Progress::SetProgress(z);

    for (Types::GridIndexType y = 0; y < dimsY; ++y)
      for (Types::GridIndexType x = 0; x < dimsX; ++x, ++offset) {
        Types::DataItem average = 0.0, weight = 0.0;
        histogram.Reset();

        Types::DataItem maskValue = 1.0;
        if (maskData) maskData->Get(maskValue, offset);

        Types::DataItem valueAvg;
        if (maskValue && averageData->Get(valueAvg, offset)) {
          // first iteration over filter: compute consistency histogram
          for (FilterMask<3>::iterator it = filter.begin(); it != filter.end();
               ++it) {
            const Types::GridIndexType xx = x + it->Location[0];
            const Types::GridIndexType yy = y + it->Location[1];
            const Types::GridIndexType zz = z + it->Location[2];

            if ((xx < dimsX) && (yy < dimsY) && (zz < dimsZ)) {
              it->Valid = true;
              // since xx, yy, zz are unsigned, we need not check
              // for >= 0; this is taken care of by overflow (we
              // hope ;)
              const size_t srcOffset = volume->GetOffsetFromIndex(xx, yy, zz);
              Types::DataItem valueAvgSrc, valueSubj;
              if (averageData->Get(valueAvgSrc, srcOffset)) {
                const size_t binAvg = histogram.ValueToBinX(valueAvgSrc);
                for (std::list<TypedArray::SmartPtr>::iterator itImg =
                         imgList.begin();
                     itImg != imgList.end(); ++itImg) {
                  if ((*itImg)->Get(valueSubj, srcOffset)) {
                    histogram.Increment(binAvg,
                                        histogram.ValueToBinY(valueSubj));
                  }
                }
              }
            }
          }

          Histogram<Types::DataItem> *avgHistogram = histogram.GetMarginalX();
          const size_t binX = histogram.ValueToBinX(valueAvg);

          for (FilterMask<3>::iterator it = filter.begin(); it != filter.end();
               ++it) {
            const Types::GridIndexType xx = x + it->Location[0];
            const Types::GridIndexType yy = y + it->Location[1];
            const Types::GridIndexType zz = z + it->Location[2];

            if (it->Valid) {
              it->Valid = false;
              // since xx, yy, zz are unsigned, we need not check for
              // >= 0; this is taken care of by overflow (we hope ;)
              const size_t srcOffset = volume->GetOffsetFromIndex(xx, yy, zz);
              Types::DataItem value;
              if (inputData->Get(value, srcOffset)) {
                float prob = static_cast<float>(it->Coefficient);

                std::list<TypedArray::SmartPtr>::iterator subjectIt =
                    subjectData.begin();
                while (subjectIt != subjectData.end()) {
                  Types::DataItem valueSubj;
                  if ((*subjectIt)->Get(valueSubj, srcOffset)) {
                    const size_t binY = histogram.ValueToBinY(valueSubj);
                    prob *= static_cast<float>(histogram.GetBin(binX, binY) /
                                               (*avgHistogram)[binX]);
                  }
                  ++subjectIt;
                }

                average += value * prob;
                weight += prob;
              }
            }
          }

          delete avgHistogram;
        }

        if (weight > 0.0) {
          filtered->Set(average / weight, offset);
        } else {
          filtered->SetPaddingAt(offset);
        }
      }
  }

  Progress::Done();

  return filtered;
}

}  // namespace cmtk
