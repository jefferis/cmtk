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

#include <vector>

#include <cmtkTypedArrayFunctionHistogramMatching.h>

cmtk::TypedArrayFunctionHistogramMatching
::TypedArrayFunctionHistogramMatching
( const TypedArray& variableArray, const TypedArray& fixedArray, const size_t numberOfHistogramBins )
  : m_Lookup( numberOfHistogramBins )
{  
  this->m_FixedArrayHistogram = Histogram<unsigned int>::SmartPtr( fixedArray.GetHistogram( numberOfHistogramBins ) );
  this->m_FixedArrayHistogram->ConvertToCumulative();
  
  this->m_VariableArrayHistogram = Histogram<unsigned int>::SmartPtr( variableArray.GetHistogram( numberOfHistogramBins ) );
  this->m_VariableArrayHistogram->ConvertToCumulative();
  
  this->CreateLookup();
}

void
cmtk::TypedArrayFunctionHistogramMatching
::CreateLookup()
{
  const size_t variableNumBins = this->m_VariableArrayHistogram->GetNumBins();
  std::vector<double> normalizedVariableHistogram( variableNumBins );
  for ( size_t l = 0; l < variableNumBins; ++l )
    {
    normalizedVariableHistogram[l] =  1.0 * this->m_VariableArrayHistogram->GetBin(l) / this->m_VariableArrayHistogram->GetBin( variableNumBins-1 );
    }
  
  const size_t fixedNumBins = this->m_FixedArrayHistogram->GetNumBins();
  std::vector<double> normalizedFixedHistogram( fixedNumBins );
  for ( size_t l = 0; l < fixedNumBins; ++l )
    {
    normalizedFixedHistogram[l] = 1.0 * this->m_FixedArrayHistogram->GetBin(l) / this->m_FixedArrayHistogram->GetBin( fixedNumBins-1 );
    }
  
  size_t j = 0;
  for ( size_t i = 0; i < variableNumBins; ++i )
    {
    while ((j < fixedNumBins) && (normalizedFixedHistogram[j] < normalizedVariableHistogram[i]))
      {
      ++j;
      }
    this->m_Lookup[i] = j;
    }
}
  
cmtk::Types::DataItem 
cmtk::TypedArrayFunctionHistogramMatching
::operator()( const cmtk::Types::DataItem valueIn ) const
{
  return this->m_FixedArrayHistogram->BinToValue( this->m_Lookup[ this->m_VariableArrayHistogram->ValueToBin( valueIn ) ] );
}
