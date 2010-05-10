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

#include <cmtkMultiChannelRegistrationFunctionalBase.h>

#include <cmtkException.h>
#include <cmtkMathUtil.h>
#include <cmtkTypes.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

void
MultiChannelRegistrationFunctionalBase
::ClearAllChannels()
{
  this->m_ReferenceChannels.resize( 0 );
  this->m_FloatingChannels.resize( 0 );
}

void
MultiChannelRegistrationFunctionalBase
::AddReferenceChannel( UniformVolume::SmartPtr& channel )
{
  if ( this->m_ReferenceChannels.size() )
    {
    this->VerifyImageSize( this->m_ReferenceChannels[0], channel );
    }
  else
    {
    this->m_ReferenceDims = channel->GetDims();
    this->m_ReferenceSize = channel->Size;
    this->m_ReferenceCropRegion = channel->CropRegion();
    }
  this->m_ReferenceChannels.push_back( channel );
  this->m_NumberOfChannels = this->m_ReferenceChannels.size() + this->m_FloatingChannels.size();

  if ( this->m_ReferenceChannels.size() == 1 )
    {
    this->NewReferenceChannelGeometry();
    }
}

void
MultiChannelRegistrationFunctionalBase
::AddFloatingChannel( UniformVolume::SmartPtr& channel )
{
  if ( this->m_FloatingChannels.size() )
    {
    this->VerifyImageSize( this->m_FloatingChannels[0], channel );
    }
  else
    {
    this->m_FloatingDims = channel->GetDims();
    this->m_FloatingSize = channel->Size;
    this->m_FloatingCropRegion = channel->GetCropRegionCoordinates();
    for ( int dim = 0; dim < 3; ++dim ) 
      {
      this->m_FloatingInverseDelta[dim] = 1.0 / channel->m_Delta[dim];
      }
    }
  this->m_FloatingChannels.push_back( channel );
  this->m_NumberOfChannels = this->m_ReferenceChannels.size() + this->m_FloatingChannels.size();
}

void
MultiChannelRegistrationFunctionalBase
::VerifyImageSize( const UniformVolume* imgA, const UniformVolume* imgB )
{
  for ( int dim  = 0; dim < 3; ++dim )
    {
    if ( imgA->GetDims()[dim] != imgB->GetDims()[dim] )
      {
      throw Exception( "MultiChannelRegistrationFunctionalBase::VerifyImageSize(): Image dimension mismatch" );
      }
    if ( fabs( imgA->Size[dim] - imgB->Size[dim] ) > 1e-6 )
      {
      throw Exception( "MultiChannelRegistrationFunctionalBase::VerifyImageSize(): Image size mismatch" );
      }
    }    
}

} // namespace cmtk
