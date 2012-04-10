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

#include <Base/cmtkMathUtil.h>

#include <vector>

template<class THistogram>
cmtk::HistogramOtsuThreshold<THistogram>::HistogramOtsuThreshold( const Self::HistogramType& histogram )
{
  const size_t nBins = histogram.GetNumberOfBins();

  std::vector<Types::DataItem> cumulativeProb( nBins );
  std::vector<Types::DataItem> cumulativeMean( nBins );

  // create cumulative mean and probability tables
  const Types::DataItem invTotal = 1.0 / histogram.SampleCount();
  cumulativeProb[0] = invTotal * histogram[0];
  cumulativeMean[0] = cumulativeProb[0] * histogram.BinToValue( 0 );
  for ( size_t i = 1; i < nBins; ++i )
    {
    const Types::DataItem p = invTotal * histogram[i];
    cumulativeProb[i] = cumulativeProb[i-1] + p;
    cumulativeMean[i] = cumulativeMean[i-1] + (p * i);
    }

  // find threshold the maximizes inter-class variance
  Types::DataItem maxSigma = 0;
  size_t maxIndex = 0;
  for ( size_t i = 0; i < nBins-1; ++i )
    {
    const Types::DataItem w1 = cumulativeProb[i];
    const Types::DataItem w2 = 1-cumulativeProb[i];

    const Types::DataItem mu1 = cumulativeMean[i] / w1;
    const Types::DataItem mu2 = (cumulativeMean[nBins-1] - cumulativeMean[i]) / w2;
    const Types::DataItem muT = cumulativeMean[nBins-1];

    const Types::DataItem sigma = w1*MathUtil::Square(mu1-muT) + w2*MathUtil::Square(mu2-muT); 
    if ( sigma > maxSigma )
      {
      maxSigma = sigma;
      maxIndex = i;
      }
    }

  this->m_Threshold = histogram.BinToValue( maxIndex );
}
