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

#ifndef __cmtkAnalyze_h_included_
#define __cmtkAnalyze_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// IDs for data types in Analyze image file.
typedef enum {
  ANALYZE_TYPE_NONE = 0,
  ANALYZE_TYPE_BINARY =  1,
  ANALYZE_TYPE_UNSIGNED_CHAR = 2,
  ANALYZE_TYPE_SIGNED_SHORT = 4,
  ANALYZE_TYPE_SIGNED_INT = 8,
  ANALYZE_TYPE_FLOAT = 16,
  ANALYZE_TYPE_COMPLEX = 32,
  ANALYZE_TYPE_DOUBLE = 64,
  ANALYZE_TYPE_RGB = 128,
  ANALYZE_TYPE_USHORT = 132, //SPM extension
  ANALYZE_TYPE_UINT = 136, // SPM extension
  ANALYZE_TYPE_ALL = 255
} AnalyzeDataType;

/// IDs for slice orientations in Analyze image file.
typedef enum {
  ANALYZE_AXIAL = 0,
  ANALYZE_CORONAL = 1,
  ANALYZE_SAGITTAL = 2,
  ANALYZE_AXIAL_FLIP = 3,
  ANALYZE_CORONAL_FLIP = 4,
  ANALYZE_SAGITTAL_FLIP = 5,
  ANALYZE_UNKNOWN = 255
} AnalyzeOrientation;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkAnalyze_h_included_
