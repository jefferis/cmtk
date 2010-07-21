/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkTypedArrayFunctionHistogramMatching_h_included_
#define __cmtkTypedArrayFunctionHistogramMatching_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkTypedArrayFunction.h"
#include "Base/cmtkTypedArray.h"
#include "Base/cmtkHistogram.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Lookup class for histogram intensity matching.
 * This class provides a lookup table that is computed from the histograms of two
 * cmtk::TypedArray objects. The lookup can then be applied to the "variable" array
 * so that its distribution afterwards matches, as closely as possible, the distribution
 * of the "fixed" array.
 *
 * To apply histogram matching to "variableArray" based on the distribution of "fixedArray",
 * use the following:
 *\code
 *  variableArray->ApplyFunction( cmtk::TypedArrayFunctionHistogramMatching( variableArray, fixedArray ) );
 *\endcode
 * The variable array for setting up the matching function need not be the same as the array
 * that the function is applied to:
 *\code
 *  variableArray->ApplyFunction( cmtk::TypedArrayFunctionHistogramMatching( testArray, fixedArray ) );
 *\endcode
 */
class
TypedArrayFunctionHistogramMatching
/// Inherit from base class.
  : public TypedArrayFunction
{
public:
  /// This class.
  typedef TypedArrayFunctionHistogramMatching Self;

  /// Default number of histogram bins.
  static const size_t DefaultNumberOfHistogramBins = 1024;

  /// Histogram type.
  typedef Histogram<unsigned int> HistogramType;

  /// Constructor: build lookup from two data arrays.
  TypedArrayFunctionHistogramMatching( const TypedArray& variableArray, const TypedArray& fixedArray, const size_t numberOfHistogramBins = Self::DefaultNumberOfHistogramBins );

  /// Constructor: build lookup from two histograms.
  TypedArrayFunctionHistogramMatching( const Self::HistogramType& variableHistogram, const Self::HistogramType& fixedHistogram );

  /// Map a single value from the variable array to its new value.
  virtual Types::DataItem operator()( const Types::DataItem valueIn ) const;

private:
  /// Fixed array histogram.
  HistogramType::SmartPtr m_FixedArrayHistogram;

  /// Variable array histogram.
  HistogramType::SmartPtr m_VariableArrayHistogram;

  /// Lookup table that translates between the two histograms.
  std::vector<unsigned int> m_Lookup;

  /// Create lookup from histograms.
  void CreateLookup();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypedArrayFunctionHistogramMatching_h_included_
