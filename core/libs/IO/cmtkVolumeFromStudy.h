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

#ifndef __cmtkVolumeFromStudy_h_included_
#define __cmtkVolumeFromStudy_h_included_

#include <cmtkconfig.h>

#include <cmtkVolumeFromSlices.h>

#include <cmtkStudy.h>
#include <cmtkStudyImageSet.h>
#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Class for building 3D volumes from an Study object.
class VolumeFromStudy : 
  /// Inherit basic slice-to-volume functions.
  public VolumeFromSlices 
{
public:
  /** Build volume from slice images.
   *@see VolumeFromSlices#AssembleVolume
   */
  virtual const UniformVolume::SmartPtr AssembleVolume ( const StudyImageSet* study, const bool verbose = false );

  /// Read from generic Study object.
  static const UniformVolume::SmartPtr Read( const Study* study, const bool verbose = false );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVolumeFromStudy_h_included_
