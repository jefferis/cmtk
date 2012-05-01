/*
//
//  Copyright 2012 SRI International
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

template<class THistogram>
cmtk::HistogramThresholdByVolume<THistogram>::HistogramThresholdByVolume( const typename Self::HistogramType& histogram, const typename Self::HistogramType::BinType volumeAbove )
{
  typename Self::HistogramType cumulativeHistogram = histogram;
  cumulativeHistogram.ConvertToCumulative();

  const size_t nBins = cumulativeHistogram.GetNumberOfBins();
  const typename Self::HistogramType::BinType volumeBelow = cumulativeHistogram[nBins-1] - volumeAbove;
  
  // reverse cumulative histogram - volume we want is above threshold, not below.
  for ( size_t i = 0; i < nBins; ++i )
    {
    if ( cumulativeHistogram[i] >= volumeBelow )
      {
      this->m_Threshold = histogram.BinToValue( i );
      return;
      }
    }
  
  // as a fall-back, return upper end of value range
  this->m_Threshold = histogram.GetRange().m_UpperBound;
}
