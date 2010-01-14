/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkHistogram.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T> 
Histogram<T>*
Histogram<T>
::Clone ( const bool copyData ) const 
{
  Histogram<T> *clone = 
    new Histogram<T>( this->m_NumBins, CMTK_HISTOGRAM_NORESET );
  
  if ( copyData )
    memcpy( clone->Bins, this->Bins, this->m_NumBins * sizeof( T ) );
  else
    clone->Reset();
  
  return clone;
}

template<class T> 
size_t
Histogram<T>
::GetMaximumBinIndex () const 
{
  T maximum = this->Bins[0];
  size_t maximumIndex = 0;
  
  for ( size_t i = 0; i<this->m_NumBins; ++i ) 
    {
    if ( this->Bins[ i ] > maximum ) 
      {
      maximum = this->Bins[ i ];
      maximumIndex = i;
      }
    }
  return maximumIndex;
}

template<class T>
double
Histogram<T>
::GetEntropy() const 
{
  double H = 0;
  
  const T sampleCount = this->SampleCount();
  if ( ! sampleCount ) 
    return CMTK_DOUBLE_NAN;
  
  for ( size_t i=0; i<this->m_NumBins; ++i ) 
    {
    if ( this->Bins[i] ) 
      {
      const double pX = ((double)this->Bins[i]) / sampleCount;
      H -= pX*log(pX);
      }
    }
  return H;
}

template<class T> 
void
Histogram<T>
::AddHistogram
( const Self& other )
{
  assert( this->m_NumBins == other.m_NumBins );
  
  for ( size_t i = 0; i<this->m_NumBins; ++i ) 
    {
    this->Bins[i] += other.Bins[i];
    }
}

template<class T> 
void 
Histogram<T>
::RemoveHistogram
( const Self& other ) 
{
  assert( this->m_NumBins == other.m_NumBins );
  
  for ( size_t i = 0; i<this->m_NumBins; ++i ) 
    {
    assert( this->Bins[i] >= other.Bins[i] );
    this->Bins[i] -= other.Bins[i];
    }
}

template<class T>
void
Histogram<T>
::Normalize
( const T normalizeTo ) 
{
  T sampleCount = this->SampleCount();
  for ( size_t i = 0; i < this->m_NumBins; ++i )
    ( this->Bins[ i ] *= normalizeTo ) /= sampleCount;
}

template<class T>
void
Histogram<T>
::NormalizeMaximum
( const T normalizeTo ) 
{
  T maximum = this->GetMaximumBinValue();
  for ( size_t i = 0; i < this->m_NumBins; ++i )
    ( this->Bins[ i ] *= normalizeTo ) /= maximum;
}

template<class T>
Types::DataItem 
Histogram<T>
::GetPercentile( const Types::DataItem percentile ) const
{
  const Types::DataItem threshold = percentile * this->SampleCount();
  Types::DataItem cumulative = 0;
  for ( size_t i = 0; i < this->m_NumBins; ++i )
    {
    cumulative += this->GetBin( i );
    if ( cumulative >= threshold )
      return this->BinToValue( i );
    }

  return this->m_BinsLowerBound + this->m_BinWidth * (this->m_NumBins - 1);
}


template class Histogram<int>;
template class Histogram<unsigned int>;
template class Histogram<long int>;
template class Histogram<float>;
template class Histogram<double>;

} // namespace cmtk
