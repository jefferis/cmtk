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

#ifndef __cmtkTypedArrayHistogramEqualizationLookup_h_included_
#define __cmtkTypedArrayHistogramEqualizationLookup_h_included_

#include <cmtkconfig.h>

#include <cmtkTypedArrayLookup.h>

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
TypedArrayHistogramEqualizationLookup
{
public:
  /// Constructor: build lookup.
  TypedArrayHistogramEqualizationLookup( const TypedArray& variableArray, const size_t numberOfHistogramBins );

  /// Map a single value from the variable array to its new value.
  virtual Types::DataItem MapSingleValue( const Types::DataItem valueIn ) const;

private:
  /// Data histogram.
  Histogram<unsigned int>::SmartPtr m_Histogram;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypedArrayHistogramEqualizationLookup_h_included_
