/*
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkTypedArrayFunctionHistogramMatching.h"

#include <vector>

cmtk::TypedArrayFunctionHistogramMatching
::TypedArrayFunctionHistogramMatching
( const TypedArray& variableArray, const TypedArray& fixedArray, const size_t numberOfHistogramBins )
  : m_Lookup( numberOfHistogramBins )
{  
  this->m_FixedArrayHistogram = fixedArray.GetHistogram( numberOfHistogramBins, true /*centeredBins*/ );
  this->m_FixedArrayHistogram->ConvertToCumulative();
  
  this->m_VariableArrayHistogram = variableArray.GetHistogram( numberOfHistogramBins, true /*centeredBins*/ );
  this->m_VariableArrayHistogram->ConvertToCumulative();
  
  this->CreateLookup();
}

cmtk::TypedArrayFunctionHistogramMatching
::TypedArrayFunctionHistogramMatching( const Self::HistogramType& variableHistogram, const Self::HistogramType& fixedHistogram )
  : m_Lookup( variableHistogram.GetNumberOfBins() )
{
  this->m_FixedArrayHistogram = Self::HistogramType::SmartPtr( fixedHistogram.Clone() );
  this->m_FixedArrayHistogram->ConvertToCumulative();
  
  this->m_VariableArrayHistogram = Self::HistogramType::SmartPtr( variableHistogram.Clone() );
  this->m_VariableArrayHistogram->ConvertToCumulative();
  
  this->CreateLookup();
}

void
cmtk::TypedArrayFunctionHistogramMatching
::CreateLookup()
{
  const size_t variableNumBins = this->m_VariableArrayHistogram->GetNumberOfBins();
  std::vector<double> normalizedVariableHistogram( variableNumBins );
  for ( size_t l = 0; l < variableNumBins; ++l )
    {
    normalizedVariableHistogram[l] =  1.0 * (*(this->m_VariableArrayHistogram))[l] / (*(this->m_VariableArrayHistogram))[variableNumBins-1];
    }
  
  const size_t fixedNumBins = this->m_FixedArrayHistogram->GetNumberOfBins();
  std::vector<double> normalizedFixedHistogram( fixedNumBins );
  for ( size_t l = 0; l < fixedNumBins; ++l )
    {
    normalizedFixedHistogram[l] = 1.0 * (*(this->m_FixedArrayHistogram))[l] / (*(this->m_FixedArrayHistogram))[fixedNumBins-1];
    }
  
  this->m_Lookup[0] = 0;
  
  size_t j = 0;
  for ( size_t i = 1; i < variableNumBins; ++i )
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
