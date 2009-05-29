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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
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

  if ( !Volume.IsNull() && reRead ) 
    {
    oldVolume = Volume;
    Volume = UniformVolume::SmartPtr( NULL );
    }

  if ( Volume.IsNull() ) 
    {
    Volume = UniformVolume::SmartPtr( VolumeFromStudy::Read( this ) );
    if ( Volume ) 
      {
      this->SetDims( Volume->GetDims()[0], Volume->GetDims()[1], Volume->GetDims()[2] );
      DisplayedImageIndex = Volume->GetDims( AXIS_Z ) / 2 ;
      ZoomFactor = 1;
      const TypedArray *dataArray = Volume->GetData();
      if ( dataArray ) 
	{
	dataArray->GetRange( MinimumValue, MaximumValue );
	Black = MinimumValue;
	White = MaximumValue;
	StandardColormap = 0;
	ReverseColormap = false;
	}
      }
    }

  if ( !Volume.IsNull() && !Volume->GetData().IsNull() ) 
    {
    return true;
    }

  Volume = oldVolume;
  return false;
}

} // namespace cmtk
