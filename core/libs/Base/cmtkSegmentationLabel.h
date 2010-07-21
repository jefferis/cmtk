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

#ifndef __cmtkSegmentationLabel_h_included_
#define __cmtkSegmentationLabel_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkMacros.h"
#include "Base/cmtkTypes.h"

#include <map>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// A label (class) in a segmentation.
class SegmentationLabel
{
public:
  /// Init constructor.
  SegmentationLabel() { this->m_Name = NULL; this->m_RGB[0] = this->m_RGB[1] = this->m_RGB[2] = 0; }

  /// Destructor.
  ~SegmentationLabel() { if ( this->m_Name ) free( this->m_Name ); }

  /// Name of this label.
  cmtkGetSetMacroString(Name);

  /// Color as RGB components for visualization.
  cmtkGetSetMacro3Array(byte,RGB);
};

/// Map from numerical IDs to labels.
typedef std::map<int,SegmentationLabel> SegmentationLabelMap;

/// Create system label map.
void CreateSystemLabelColorMap( SegmentationLabelMap& map );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSegmentationLabel_h_included_

