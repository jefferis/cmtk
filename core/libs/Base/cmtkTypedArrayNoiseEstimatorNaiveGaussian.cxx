/*
//
//  Copyright 2008-2011 SRI International
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

#include "cmtkTypedArrayNoiseEstimatorNaiveGaussian.h"

namespace
cmtk
{

TypedArrayNoiseEstimatorNaiveGaussian::TypedArrayNoiseEstimatorNaiveGaussian
( const TypedArray& data, const size_t histogramBins )
{
  Histogram<unsigned int>::SmartPtr histogram( data.GetHistogram( histogramBins ) );
  
  // find first maximum
  size_t idx = 0;
  while ( (idx < (histogramBins-1)) && ( (*histogram)[idx] <= (*histogram)[idx+1] ) )
    {
    ++idx;
    }
  
  const Types::DataItem noiseMean = histogram->BinToValue( idx );
  
  // now find following minimum
  while ( (idx < (histogramBins-1)) && ( (*histogram)[idx] > (*histogram)[idx+1] ) )
    {
    ++idx;
    }
  
  // then, compute standard deviation of all values below that threshold from
  // first maximum.
  this->m_Threshold = histogram->BinToValue( idx );

  Types::DataItem sdev = 0;
  size_t count = 0;
  for ( size_t i = 0; i < data.GetDataSize(); ++i )
    {
    Types::DataItem value;
    if ( data.Get( value, i ) )
      {
      if ( value <= this->m_Threshold )
	{
        sdev += static_cast<Types::DataItem>( MathUtil::Square( value - noiseMean ) );
	++count;
	}
      }
    }

  if ( count )
    this->m_NoiseLevelSigma = static_cast<Types::DataItem>( sqrt( sdev/count ) );
  else
    this->m_NoiseLevelSigma = 0.0;
}

} // namespace cmtk
