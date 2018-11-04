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

#ifndef __cmtkHistogramBase_h_included_
#define __cmtkHistogramBase_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkTypes.h>

#include <System/cmtkSmartPtr.h>

#include <algorithm>

#ifndef CMTK_HISTOGRAM_AUTOBINS
/// Constant for number of bins to be determined automatically.
#define CMTK_HISTOGRAM_AUTOBINS 0
#endif

/// Constant for affirmative fractional bin flag.
#define CMTK_HISTOGRAM_FRACTIONAL true

/// Constant for affirmative fractional bin flag.
#define CMTK_HISTOGRAM_DISCRETE false

/// Constant for resetting histogram bins on instantiation.
#define CMTK_HISTOGRAM_RESET true

/// Constant for not resetting histogram bins on instantiation.
#define CMTK_HISTOGRAM_NORESET false

/// Constant for copying histogram bin values on object duplication.
#define CMTK_HISTOGRAM_COPY true

/// Constant for not copying histogram bin values on object duplication.
#define CMTK_HISTOGRAM_NOCOPY false

namespace cmtk {

/** \addtogroup Base */
//@{

/** Common (non-template) base class for all 1-D histograms.
 */
class HistogramBase {
 protected:
  /// Width of data bins.
  Types::DataItem m_BinWidth;

  /// Lower value bound of data bins.
  Types::DataItem m_BinsLowerBound;

  /// Upper value bound of data bins.
  Types::DataItem m_BinsUpperBound;

 public:
  /// This class.
  typedef HistogramBase Self;

  /// Default constructor.
  HistogramBase() {
    this->m_BinWidth = 1.0;
    this->m_BinsLowerBound = this->m_BinsUpperBound = 0.0;
  }

  /// Virtual destructor.
  virtual ~HistogramBase() {}

  /// Return number of histogram bins.
  virtual size_t GetNumberOfBins() const = 0;

  /** Set data range corresponding to this histogram.
   */
  void SetRange(const Types::DataItemRange &range) {
    this->m_BinsLowerBound = range.m_LowerBound;
    this->m_BinsUpperBound = range.m_UpperBound;
    this->m_BinWidth = range.Width() / (this->GetNumberOfBins() - 1);
  }

  /** Set data range corresponding to this histogram with upper and lower bound
   * centered in first and last bin.
   */
  void SetRangeCentered(const Types::DataItemRange &range) {
    this->m_BinWidth = range.Width() / (this->GetNumberOfBins() - 1);
    this->m_BinsLowerBound = static_cast<Types::DataItem>(
        range.m_LowerBound - 0.5 * this->m_BinWidth);
    this->m_BinsUpperBound = static_cast<Types::DataItem>(
        range.m_UpperBound + 0.5 * this->m_BinWidth);
  }

  /** Get value range of the distribution.
   */
  const Types::DataItemRange GetRange() const {
    return Types::DataItemRange(this->m_BinsLowerBound, this->m_BinsUpperBound);
  }

  /** Get value range of a given bin.
   */
  virtual const Types::DataItemRange GetRangeBin(const size_t bin) const {
    const Types::DataItem from =
        this->m_BinsLowerBound + this->m_BinWidth * bin;
    return Types::DataItemRange(from, from + this->m_BinWidth);
  }

  /// Get bin width.
  Types::DataItem GetBinWidth() const { return this->m_BinWidth; }

  /** Return bin corresponding to a certain value of the distribution.
   *\param value A value from the distribution.
   *\return The index of the bin corresponding to the given value.
   */
  virtual size_t ValueToBin(const Types::DataItem value) const {
    const size_t binIndex = static_cast<size_t>(
        (value - this->m_BinsLowerBound) / this->m_BinWidth);
    return std::max<size_t>(0, std::min(this->GetNumberOfBins() - 1, binIndex));
  }

  /** Return fractional bin corresponding to a value of the distribution.
   *\param value A value from the distribution.
   *\return The index of the fractional bin index corresponding to the given
   * value. This value is an integer if and only if the given value is
   * identical to the lower bound of a bin.
   */
  virtual Types::DataItem ValueToBinFractional(
      const Types::DataItem value) const {
    const Types::DataItem binIndex =
        (value - this->m_BinsLowerBound) / this->m_BinWidth;
    return std::max<Types::DataItem>(
        0, std::min<Types::DataItem>(
               static_cast<Types::DataItem>(this->GetNumberOfBins() - 1),
               binIndex));
  }

  /** Return center of values represented by a certain bin.
   *\param bin Index of a bin from the distribution.
   *\return Average of upper and lower margin values of the given bin.
   */
  virtual Types::DataItem BinToValue(const size_t bin) const {
    return static_cast<Types::DataItem>(this->m_BinsLowerBound +
                                        (bin + 0.5) * this->m_BinWidth);
  }
};

//@}

}  // namespace cmtk

#endif  // #ifndef __cmtkHistogramBase_h_included_
