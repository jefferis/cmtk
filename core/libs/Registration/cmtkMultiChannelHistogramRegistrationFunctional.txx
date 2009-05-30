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

#include <cmtkException.h>
#include <cmtkMathUtil.h>
#include <cmtkTypes.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
void
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>
::ClearAllChannels()
{
  this->m_HashKeyScaleRef.resize( 0 );
  this->m_HashKeyOffsRef.resize( 0 );

  this->m_HashKeyScaleFlt.resize( 0 );
  this->m_HashKeyOffsFlt.resize( 0 );

  this->Superclass::ClearAllChannels();
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
void
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>
::AddReferenceChannel( UniformVolume::SmartPtr& channel )
{
  Types::DataItem min, max;
  const Types::DataItem maxBinIndex = (1<<NBitsPerChannel) - 1;

  channel->GetData()->GetRange( min, max );

  const Types::DataItem scale = maxBinIndex / (max-min);
  const Types::DataItem offset = -(min/scale);

  this->m_HashKeyScaleRef.push_back( scale );
  this->m_HashKeyOffsRef.push_back( offset );
  this->m_HashKeyShiftRef = NBitsPerChannel*this->m_ReferenceChannels.size();

  this->Superclass::AddReferenceChannel( channel );

  const size_t hashKeyBits = 8 * sizeof( THashKeyType );
  if ( this->m_NumberOfChannels * NBitsPerChannel > hashKeyBits )
    {
    StdErr << "ERROR in MultiChannelHistogramRegistrationFunctional:\n"
	      << "  Cannot represent total of " << this->m_NumberOfChannels << " channels with "
	      << NBitsPerChannel << " bits per channel using hash key type with "
	      << hashKeyBits << "bits.\n";
    exit( 1 );
    }
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
void
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>
::AddFloatingChannel( UniformVolume::SmartPtr& channel )
{
  Types::DataItem min, max;
  const Types::DataItem maxBinIndex = (1<<NBitsPerChannel) - 1;

  channel->GetData()->GetRange( min, max );

  const Types::DataItem scale = maxBinIndex / (max-min);
  const Types::DataItem offset = -(min/scale);

  this->m_HashKeyScaleFlt.push_back( scale );
  this->m_HashKeyOffsFlt.push_back( offset );

  this->Superclass::AddFloatingChannel( channel );

  const size_t hashKeyBits = 8 * sizeof( THashKeyType );
  if ( this->m_NumberOfChannels * NBitsPerChannel > hashKeyBits )
    {
    StdErr << "ERROR in MultiChannelHistogramRegistrationFunctional:\n"
	      << "  Cannot represent total of " << this->m_NumberOfChannels << " channels with "
	      << this->m_HistogramBitsPerChannel << " bits per channel using hash key type with "
	      << hashKeyBits << "bits.\n";
    exit( 1 );
    }
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
void
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>
::ContinueMetric( MetricData& metricData, const size_t rindex, const Vector3D& fvector )
{
#ifdef CMTK_VAR_AUTO_ARRAYSIZE
  Types::DataItem values[ this->m_NumberOfChannels ];
#else
  std::vector<Types::DataItem> values( this->m_NumberOfChannels );
#endif
  
  size_t idx = 0;
  for ( size_t ref = 0; ref < this->m_ReferenceChannels.size(); ++ref )
    {
    if ( !this->m_ReferenceChannels[ref]->GetDataAt( values[idx++], rindex ) ) return;
    }
  
  for ( size_t flt = 0; flt < this->m_FloatingChannels.size(); ++flt )
    {
    if ( !this->m_FloatingInterpolators[flt]->GetDataAt( fvector, values[idx++] ) ) return;
    }

  metricData += &(values[0]);
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
Functional::ReturnType
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>
::GetMetric( const MetricData& metricData ) const
{
  if ( metricData.m_TotalNumberOfSamples )
    {
    const double norm = 1.0 / metricData.m_TotalNumberOfSamples;

    double hXY = 0;
    typename MetricData::HashTableType::const_iterator it = metricData.m_JointHash.begin();
    for ( ; it != metricData.m_JointHash.end(); ++it )
      {
      if ( it->second )
	{
	const double p = norm * it->second;
	hXY -= p * log( p );
	}
      }

    double hX = 0;
    it = metricData.m_ReferenceHash.begin();
    for ( ; it != metricData.m_ReferenceHash.end(); ++it )
      {
      if ( it->second )
	{
	const double p = norm * it->second;
	hX -= p * log( p );
	}
      }

    double hY = 0;
    it = metricData.m_FloatingHash.begin();
    for ( ; it != metricData.m_FloatingHash.end(); ++it )
      {
      if ( it->second )
	{
	const double p = norm * it->second;
	hY -= p * log( p );
	}
      }
    
    if ( this->m_NormalizedMI )
      return static_cast<Functional::ReturnType>( (hX + hY) / hXY );
    else
      return static_cast<Functional::ReturnType>( hX + hY - hXY );
    }
  
  return static_cast<Functional::ReturnType>( -FLT_MAX );
}

} // namespace cmtk
