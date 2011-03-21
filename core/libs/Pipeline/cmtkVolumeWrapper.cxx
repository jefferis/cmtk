/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <Pipeline/cmtkVolumeWrapper.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

VolumeWrapper*
VolumeWrapper::New() 
{
  return new VolumeWrapper; 
}

VolumeWrapper::VolumeWrapper()
{ 
  Volume = UniformVolume::SmartPtr::Null(); 
  this->m_AffineXform = AffineXform::SmartPtr::Null();
  this->m_WarpXform = WarpXform::SmartPtr::Null();
}

VolumeWrapper::~VolumeWrapper()
{
}


void 
VolumeWrapper::SetVolume( UniformVolume::SmartPtr& volume )
{
  if ( Volume != volume ) 
    {
    Volume = volume;
    this->UpdateModifiedTime();
    }
}

void
VolumeWrapper::SetAffineXform( AffineXform::SmartPtr& affineXform )
{
  if ( this->m_AffineXform != affineXform ) 
    {
    this->m_AffineXform = affineXform;
    this->UpdateModifiedTime();
    }
}

void
VolumeWrapper::SetWarpXform( WarpXform::SmartPtr& warpXform )
{
  if ( this->m_WarpXform != warpXform ) 
    {
    this->m_WarpXform = warpXform;
    this->UpdateModifiedTime();
    }
}

} // namespace cmtk
