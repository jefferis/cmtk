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

#include "cmtkTypedArrayFunctionHistogramEqualization.h"

#include <vector>
#include <algorithm>

cmtk::TypedArrayFunctionHistogramEqualization
::TypedArrayFunctionHistogramEqualization
( const TypedArray& variableArray, const size_t numberOfHistogramBins )
{  
  this->m_Histogram = Histogram<unsigned int>::SmartPtr( variableArray.GetHistogram( numberOfHistogramBins ) );
  (*this->m_Histogram)[0] = 0; // this effectively stretches the distribution
  this->m_Histogram->ConvertToCumulative();

  const Types::DataItemRange range = variableArray.GetRange();
  this->m_MinValue = range.m_LowerBound;
  this->m_ScaleFactor = 1.0 * range.Width() / (*this->m_Histogram)[numberOfHistogramBins-1];
}

cmtk::Types::DataItem 
cmtk::TypedArrayFunctionHistogramEqualization
::operator()( const cmtk::Types::DataItem valueIn ) const
{
  return this->m_MinValue + this->m_ScaleFactor * (*this->m_Histogram)[ this->m_Histogram->ValueToBin( valueIn ) ];
}
