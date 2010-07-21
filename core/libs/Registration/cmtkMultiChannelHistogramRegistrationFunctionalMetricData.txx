/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkMultiChannelHistogramRegistrationFunctional.h"

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
void
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData
::Init( Parent *const parent )
{
  this->m_Parent = parent;

  this->m_JointHash.clear();
  this->m_ReferenceHash.clear();
  this->m_FloatingHash.clear();

  this->m_TotalNumberOfSamples = 0;
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
typename MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData& 
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData::operator=
( const typename MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData& source )
{
  this->m_JointHash = source.m_JointHash;
  this->m_ReferenceHash = source.m_ReferenceHash;
  this->m_FloatingHash = source.m_FloatingHash;

  this->m_TotalNumberOfSamples = source.m_TotalNumberOfSamples;
  
  return *this;
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
typename MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData& 
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData::operator+=
( const typename MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData& other )
{
  for ( typename HashTableType::const_iterator it = other.m_JointHash.begin(); it != other.m_JointHash.end(); ++it )
    {
    this->m_JointHash[it->first] += it->second;
    }
  for ( typename HashTableType::const_iterator it = other.m_ReferenceHash.begin(); it != other.m_ReferenceHash.end(); ++it )
    {
    this->m_ReferenceHash[it->first] += it->second;
    }
  for ( typename HashTableType::const_iterator it = other.m_FloatingHash.begin(); it != other.m_FloatingHash.end(); ++it )
    {
    this->m_FloatingHash[it->first] += it->second;
    }
    
  this->m_TotalNumberOfSamples += other.m_TotalNumberOfSamples;
  
  return *this;
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
typename MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData& 
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData::operator-=
( const typename MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData& other )
{
  for ( typename HashTableType::const_iterator it = other.m_JointHash.begin(); it != other.m_JointHash.end(); ++it )
    {
    this->m_JointHash[it->first] -= it->second;
    }
  for ( typename HashTableType::const_iterator it = other.m_ReferenceHash.begin(); it != other.m_ReferenceHash.end(); ++it )
    {
    this->m_ReferenceHash[it->first] -= it->second;
    }
  for ( typename HashTableType::const_iterator it = other.m_FloatingHash.begin(); it != other.m_FloatingHash.end(); ++it )
    {
    this->m_FloatingHash[it->first] -= it->second;
    }
    
  this->m_TotalNumberOfSamples -= other.m_TotalNumberOfSamples;
  
  return *this;
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
void
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData::operator+=
( const Types::DataItem* values )
{
  THashKeyType hashKeyRef = 0, hashKeyFlt = 0;
  size_t idx = 0;
  for ( size_t ref = 0; ref < m_Parent->m_ReferenceChannels.size(); ++ref, ++idx )
    {
    hashKeyRef |= static_cast<THashKeyType>(m_Parent->m_HashKeyScaleRef[ref] * values[idx] + m_Parent->m_HashKeyOffsRef[ref] ) << (NBitsPerChannel*ref);
    }
  for ( size_t flt = 0; flt < m_Parent->m_FloatingChannels.size(); ++flt, ++idx )
    {
    hashKeyFlt |= static_cast<THashKeyType>(m_Parent->m_HashKeyScaleFlt[flt] * values[idx] + m_Parent->m_HashKeyOffsFlt[flt] ) <<( NBitsPerChannel*flt);
    }
  
  THashKeyType hashKeyJnt = (hashKeyFlt << m_Parent->m_HashKeyShiftRef) + hashKeyRef;

  ++this->m_ReferenceHash[hashKeyRef];
  ++this->m_FloatingHash[hashKeyFlt];
  ++this->m_JointHash[hashKeyJnt];

  ++this->m_TotalNumberOfSamples;
}

template<class TDataType,class TInterpolator,class THashKeyType,char NBitsPerChannel>
void
MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel>::MetricData::operator-=
( const Types::DataItem* values )
{
  THashKeyType hashKeyRef = 0, hashKeyFlt = 0;
  size_t idx = 0;
  for ( size_t ref = 0; ref < m_Parent->m_ReferenceChannels.size(); ++ref, ++idx )
    {
    hashKeyRef |= static_cast<THashKeyType>(m_Parent->m_HashKeyScaleRef[ref] * values[idx] + m_Parent->m_HashKeyOffsRef[ref] ) << (NBitsPerChannel*ref);
    }
  for ( size_t flt = 0; flt < m_Parent->m_FloatingChannels.size(); ++flt, ++idx )
    {
    hashKeyFlt |= static_cast<THashKeyType>(m_Parent->m_HashKeyScaleFlt[flt] * values[idx] + m_Parent->m_HashKeyOffsFlt[flt] ) << (NBitsPerChannel*flt);
    }
  
  THashKeyType hashKeyJnt = (hashKeyFlt << m_Parent->m_HashKeyShiftRef) + hashKeyRef;

  --this->m_ReferenceHash[hashKeyRef];
  --this->m_FloatingHash[hashKeyFlt];
  --this->m_JointHash[hashKeyJnt];

  --this->m_TotalNumberOfSamples;
}

} // namespace cmtk
