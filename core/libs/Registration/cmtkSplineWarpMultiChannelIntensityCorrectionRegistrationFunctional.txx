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

#include <limits>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TMetricFunctional>
void
SplineWarpMultiChannelIntensityCorrectionRegistrationFunctional<TMetricFunctional>
::ContinueMetric( MetricData& metricData, const size_t rindex, const Vector3D& fvector )
{
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Types::DataItem values[ this->m_NumberOfChannels ];
#else
  std::vector<Types::DataItem> values( this->m_NumberOfChannels );
#endif
  
  size_t idx = 0;
  for ( size_t ref = 0; ref < this->m_ReferenceChannels.size(); ++ref )
    {
    if ( !this->m_ReferenceChannels[ref]->GetDataAt( values[idx++], rindex ) ) return;
    }
  
  const int planeSize = this->m_ReferenceDims[0] * this->m_ReferenceDims[1];
  const int z = rindex / planeSize;
  const int y = (rindex % planeSize) / this->m_ReferenceDims[0];
  const int x = rindex % this->m_ReferenceDims[0];

  const Types::DataItem jacobian = static_cast<Types::DataItem>( this->m_Transformation.GetJacobianDeterminant( x, y, z ) );

  for ( size_t flt = 0; flt < this->m_FloatingChannels.size(); ++flt, ++idx )
    {
    if ( !this->m_FloatingInterpolators[flt]->GetDataAt( fvector, values[idx] ) )
      return;
    values[idx] *= jacobian;
    }
  
  metricData += &(values[0]);
}

template<class TMetricFunctional>
void
SplineWarpMultiChannelIntensityCorrectionRegistrationFunctional<TMetricFunctional>
::ContinueMetricStoreReformatted( MetricData& metricData, const size_t rindex, const Vector3D& fvector )
{
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Types::DataItem values[ this->m_NumberOfChannels ];
#else
  std::vector<Types::DataItem> values( this->m_NumberOfChannels );
#endif
  
  size_t idx = 0;
  for ( size_t ref = 0; ref < this->m_ReferenceChannels.size(); ++ref )
    {
    if ( !this->m_ReferenceChannels[ref]->GetDataAt( values[idx++], rindex ) ) return;
    }
  
  const int planeSize = this->m_ReferenceDims[0] * this->m_ReferenceDims[1];
  const int z = rindex / planeSize;
  const int y = (rindex % planeSize) / this->m_ReferenceDims[0];
  const int x = rindex % this->m_ReferenceDims[0];

  const Types::DataItem jacobian = static_cast<Types::DataItem>( this->m_Transformation.GetJacobianDeterminant( x, y, z ) );

  for ( size_t flt = 0; flt < this->m_FloatingChannels.size(); ++flt, ++idx )
    {
    if ( !this->m_FloatingInterpolators[flt]->GetDataAt( fvector, values[idx] ) )
      {
      for ( size_t f = 0; f < this->m_FloatingChannels.size(); ++f )
	this->m_ReformattedFloatingChannels[f][rindex] = std::numeric_limits<float>::quiet_NaN();
      return;
      }

    values[idx] *= jacobian;
    this->m_ReformattedFloatingChannels[flt][rindex] = static_cast<float>( values[idx] );
    }

  metricData += &(values[0]);
}

} // namespace cmtk
