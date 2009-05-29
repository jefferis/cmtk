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

#ifndef __cmtkXformIO_h_included__
#define __cmtkXformIO_h_included__

#include <cmtkconfig.h>

#include <cmtkXform.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Utility class for one-stop transformation import (later also export).
class XformIO 
{
public:
  /// This class.
  typedef XformIO Self;

  /// Read transformation from filesystem.
  static Xform* Read( const char* path, const bool verbose = false );

  /// Write transformation to filesystem.
  static void Write( const Xform* xform, const char* path, const bool verbose = false );

protected:
#ifdef CMTK_BUILD_NRRD
  /// Read deformation field from Nrrd image file.
  static Xform* ReadNrrd( const char* path, const bool verbose = false );

  /// Write transformation to filesystem.
  static void WriteNrrd( const Xform* xform, const char* path, const bool verbose = false );
#endif // #ifdef CMTK_BUILD_NRRD
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkXformIO_h_included__
