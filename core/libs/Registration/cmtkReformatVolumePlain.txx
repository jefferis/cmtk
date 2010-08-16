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

#include <Registration/cmtkReformatVolume.h>
#include <Base/cmtkUniformVolumeInterpolator.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template <class TInterpolator>
bool
ReformatVolume::Plain::operator()
  ( Types::DataItem& value, const Vector3D& inRef, const XformList& refToFloat, TInterpolator& interpolator )
{
  Vector3D inFlt( inRef );
  if ( ! refToFloat.ApplyInPlace( inFlt ) ) 
    return false;

  return interpolator->GetDataAt( inFlt, value ); 
}

} // namespace cmtk

