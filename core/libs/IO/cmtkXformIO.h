/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkXformIO_h_included__
#define __cmtkXformIO_h_included__

#include <cmtkconfig.h>

#include "Base/cmtkXform.h"

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Utility class for one-stop transformation import.
 * When reading a transformation file using the Read() function, the file and transformation type
 * are automatically detected based on each file format's "magic number".
 *
 * When writing a transformation using the Write() function, the path or file name suffix determines
 * the output file format. Supported formats are: ITK Transformation file (".txt"; ".tfm"), Nrrd deformation
 * fields (".nrrd"; ".nhdr"), and legacy TypedStream (all other suffixes).
 */
class XformIO 
{
public:
  /// This class.
  typedef XformIO Self;

  /// Read transformation from filesystem.
  static Xform::SmartPtr Read( const char* path, const bool verbose = false );

  /// Write transformation to filesystem.
  static void Write( const Xform* xform, const char* path, const bool verbose = false );

protected:
#ifdef CMTK_BUILD_NRRD
  /// Read deformation field from Nrrd image file.
  static Xform::SmartPtr ReadNrrd( const char* path, const bool verbose = false );

  /// Write transformation to filesystem.
  static void WriteNrrd( const Xform* xform, const char* path, const bool verbose = false );
#endif // #ifdef CMTK_BUILD_NRRD
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkXformIO_h_included__
