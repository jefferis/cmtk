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

#include <cmtkStudyImageSet.h>

#include <cmtkVolumeFromStudy.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

bool
StudyImageSet::ReadVolume( const bool reRead, const char* )
{
  UniformVolume::SmartPtr oldVolume( NULL );

  if ( this->m_Volume && reRead ) 
    {
    oldVolume = this->m_Volume;
    this->m_Volume = UniformVolume::SmartPtr( NULL );
    }

  if ( !this->m_Volume ) 
    {
    this->m_Volume = UniformVolume::SmartPtr( VolumeFromStudy::Read( this ) );
    if ( this->m_Volume ) 
      {
      this->SetDims( this->m_Volume->GetDims()[0], this->m_Volume->GetDims()[1], this->m_Volume->GetDims()[2] );
      this->m_DisplayedImageIndex = this->m_Volume->GetDims( AXIS_Z ) / 2 ;
      this->m_ZoomFactor = 1;
      const TypedArray *dataArray = this->m_Volume->GetData();
      if ( dataArray ) 
	{
	dataArray->GetRange( this->m_MinimumValue, this->m_MaximumValue );
	this->m_Black = this->m_MinimumValue;
	this->m_White = this->m_MaximumValue;
	this->m_StandardColormap = 0;
	this->m_ReverseColormap = false;
	}
      }
    }

  if ( this->m_Volume && this->m_Volume->GetData() ) 
    {
    return true;
    }

  this->m_Volume = oldVolume;
  return false;
}

} // namespace cmtk
