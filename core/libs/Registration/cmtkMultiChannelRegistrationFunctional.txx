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

template<class TDataType,class TInterpolator>
void
MultiChannelRegistrationFunctional<TDataType,TInterpolator>
::ClearAllChannels()
{
  this->m_ReferenceChannels.resize( 0 );
  this->m_FloatingChannels.resize( 0 );
}

template<class TDataType,class TInterpolator>
void
MultiChannelRegistrationFunctional<TDataType,TInterpolator>
::AddReferenceChannel( UniformVolume::SmartPtr& channel )
{
  if ( this->m_ReferenceChannels.size() )
    {
    this->VerifyImageSize( this->m_ReferenceChannels[0], channel );
    }
  else
    {
    memcpy( this->m_ReferenceDims, channel->GetDims(), sizeof(this->m_ReferenceDims) );
    memcpy( this->m_ReferenceSize, channel->Size, sizeof(this->m_ReferenceSize) );
    channel->GetCropRegion( this->m_ReferenceCropFrom, this->m_ReferenceCropTo );
    }
  this->m_ReferenceChannels.push_back( channel );
  this->m_NumberOfChannels = this->m_ReferenceChannels.size() + this->m_FloatingChannels.size();

  if ( this->m_ReferenceChannels.size() == 1 )
    {
    this->NewReferenceChannelGeometry();
    }
}

template<class TDataType,class TInterpolator>
void
MultiChannelRegistrationFunctional<TDataType,TInterpolator>
::AddFloatingChannel( UniformVolume::SmartPtr& channel )
{
  if ( this->m_FloatingChannels.size() )
    {
    this->VerifyImageSize( this->m_FloatingChannels[0], channel );
    }
  else
    {
    memcpy( this->m_FloatingDims, channel->GetDims(), sizeof(this->m_FloatingDims) );
    memcpy( this->m_FloatingSize, channel->Size, sizeof(this->m_FloatingSize) );
    channel->GetCropRegion( this->m_FloatingCropFrom, this->m_FloatingCropTo );
    for ( int dim = 0; dim < 3; ++dim ) 
      {
      this->m_FloatingInverseDelta[dim] = 1.0 / channel->Delta[dim];
      }
    }
  this->m_FloatingChannels.push_back( channel );
  this->m_FloatingInterpolators.push_back( typename TInterpolator::SmartPtr( new TInterpolator( channel ) ) );

  this->m_NumberOfChannels = this->m_ReferenceChannels.size() + this->m_FloatingChannels.size();
}

template<class TDataType,class TInterpolator>
void
MultiChannelRegistrationFunctional<TDataType,TInterpolator>
::VerifyImageSize( const UniformVolume* imgA, const UniformVolume* imgB )
{
  for ( int dim  = 0; dim < 3; ++dim )
    {
    if ( imgA->GetDims(dim) != imgB->GetDims(dim) )
      {
      throw Exception( "MultiChannelRegistrationFunctional::VerifyImageSize(): Image dimension mismatch" );
      }
    if ( fabs( imgA->Size[dim] - imgB->Size[dim] ) > 1e-6 )
      {
      throw Exception( "MultiChannelRegistrationFunctional::VerifyImageSize(): Image size mismatch" );
      }
    }    
}

} // namespace cmtk
