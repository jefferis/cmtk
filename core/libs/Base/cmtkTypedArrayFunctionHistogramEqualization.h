/*
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkTypedArrayFunctionHistogramEqualization_h_included_
#define __cmtkTypedArrayFunctionHistogramEqualization_h_included_

#include <cmtkconfig.h>

#include <cmtkTypedArrayFunction.h>

#include <cmtkTypedArray.h>
#include <cmtkHistogram.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Lookup class for histogram equalization.
 */
class
TypedArrayFunctionHistogramEqualization
/// Inherit from base class.
  : public TypedArrayFunction
{
public:
  /// Constructor: build lookup.
  TypedArrayFunctionHistogramEqualization( const TypedArray& variableArray, const size_t numberOfHistogramBins );

  /// Map a single value from the variable array to its new value.
  virtual Types::DataItem operator()( const Types::DataItem valueIn ) const;

private:
  /// Data histogram.
  Histogram<unsigned int>::SmartPtr m_Histogram;

  /// Scale factor from cumulative distribution to histogram.
  Types::DataItem m_ScaleFactor;

  /// Minimum data value.
  Types::DataItem m_MinValue;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypedArrayFunctionHistogramEqualization_h_included_
