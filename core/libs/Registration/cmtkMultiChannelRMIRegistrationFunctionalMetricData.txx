/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkMultiChannelRMIRegistrationFunctional.h"

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TRealType,class TDataType,class TInterpolator>
void
MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData
::Init( Parent *const parent )
{
  const size_t nref = parent->m_ReferenceChannels.size();
  const size_t nflt = parent->m_FloatingChannels.size();

  this->m_Sums.resize( nref + nflt );
  std::fill( this->m_Sums.begin(), this->m_Sums.end(), static_cast<TRealType>( 0.0 ) );

  this->m_Products.resize( ((nref+nflt) * (nref+nflt+1)) / 2 );
  std::fill( this->m_Products.begin(), this->m_Products.end(), static_cast<TRealType>( 0.0 ) );

  this->m_CovarianceMatrix.Resize( nref+nflt, nref+nflt ); // needs no reset
  this->m_CovarianceMatrixRef.Resize( nref, nref ); // needs no reset
  this->m_CovarianceMatrixFlt.Resize( nflt, nflt ); // needs no reset
  
  this->m_TotalNumberOfSamples = 0;
}

template<class TRealType,class TDataType,class TInterpolator>
typename MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData& 
MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData::operator=
( const typename MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData& source )
{
  this->m_Sums.resize( source.m_Sums.size() );
  std::copy( source.m_Sums.begin(), source.m_Sums.end(), this->m_Sums.begin() );

  this->m_Products.resize( source.m_Products.size() );
  std::copy( source.m_Products.begin(), source.m_Products.end(), this->m_Products.begin() );
    
  // covariance matrices need not be copied as they only provide temporary storage for GetMetric()

  this->m_TotalNumberOfSamples = source.m_TotalNumberOfSamples;
  
  return *this;
}

template<class TRealType,class TDataType,class TInterpolator>
typename MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData& 
MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData::operator+=
( const typename MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData& other )
{
  assert( this->m_Sums.size() == other.m_Sums.size() );
  assert( this->m_Products.size() == other.m_Products.size() );
  
  for ( size_t idx = 0; idx < this->m_Sums.size(); ++idx )
    this->m_Sums[idx] += other.m_Sums[idx];

  for ( size_t idx = 0; idx < this->m_Products.size(); ++idx )
    this->m_Products[idx] += other.m_Products[idx];

  // covariance matrices need not be treated as they only provide temporary storage for GetMetric()

  this->m_TotalNumberOfSamples += other.m_TotalNumberOfSamples;
  
  return *this;
}

template<class TRealType,class TDataType,class TInterpolator>
typename MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData& 
MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData::operator-=
( const typename MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData& other )
{
  assert( this->m_Sums.size() == other.m_Sums.size() );
  assert( this->m_Products.size() == other.m_Products.size() );
  
  for ( size_t idx = 0; idx < this->m_Sums.size(); ++idx )
    this->m_Sums[idx] -= other.m_Sums[idx];

  for ( size_t idx = 0; idx < this->m_Products.size(); ++idx )
    this->m_Products[idx] -= other.m_Products[idx];

  // covariance matrices need not be treated as they only provide temporary storage for GetMetric()

  this->m_TotalNumberOfSamples -= other.m_TotalNumberOfSamples;
  
  return *this;
}

template<class TRealType,class TDataType,class TInterpolator>
void
MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData::operator+=
( const Types::DataItem* values )
{
  const size_t numberOfChannels = this->m_Sums.size();
  for ( size_t j = 0; j < numberOfChannels; ++j )
    {
    this->m_Sums[j] += values[j];
    }
    
  size_t idx = 0;
  for ( size_t j = 0; j < numberOfChannels; ++j )
    {
    for ( size_t i = 0; i <= j; ++i, ++idx )
      {
      this->m_Products[idx] += values[i] * values[j];
      }
    }

  ++this->m_TotalNumberOfSamples;
}

template<class TRealType,class TDataType,class TInterpolator>
void
MultiChannelRMIRegistrationFunctional<TRealType,TDataType,TInterpolator>::MetricData::operator-=
( const Types::DataItem* values )
{
  const size_t numberOfChannels = this->m_Sums.size();
  for ( size_t j = 0; j < numberOfChannels; ++j )
    {
    this->m_Sums[j] -= values[j];
    }
    
  size_t idx = 0;
  for ( size_t j = 0; j < numberOfChannels; ++j )
    {
    for ( size_t i = 0; i <= j; ++i, ++idx )
      {
      this->m_Products[idx] -= values[i] * values[j];
      }
    }

  --this->m_TotalNumberOfSamples;
}

} // namespace cmtk
