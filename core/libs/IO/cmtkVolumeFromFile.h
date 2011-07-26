/*
//
//  Copyright 2004-2011 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkVolumeFromFile_h_included_
#define __cmtkVolumeFromFile_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <System/cmtkException.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Class to read uniform volume from disk file(s).
 */
class VolumeFromFile 
{
public:
  /// Read volume from file automatically.
  static const UniformVolume::SmartPtr Read( const char *filename );

#ifdef CMTK_USE_DCMTK
  /// Read volume in multi-slice DICOM format.
  static const UniformVolume::SmartPtr ReadDICOM( const char *filename );
#else
  static const UniformVolume::SmartPtr ReadDICOM( const char* )
  {
    throw Exception( "Library was configured without DICOM support." );
  }
#endif

  /// Read volume in Vanderbilt format.
  static const UniformVolume::SmartPtr ReadVanderbilt( const char *filename );

  /// Read BioRad PIC file (confocal microscopy image).
  static const UniformVolume::SmartPtr ReadBioRad( const char* path );

  /// Read Analyze 7.5 file (separate header file).
  static const UniformVolume::SmartPtr ReadAnalyzeHdr( const char* pathHdr, const bool bigEndian = false, const bool readData = true );

  /// Read Nifti file.
  static const UniformVolume::SmartPtr ReadNifti( const char* pathHdr, const bool detached, const bool readData = true );

  /// Write Analyze 7.5 file (separate header file).
  static void WriteAnalyzeHdr( const char* pathHdr, const UniformVolume& volume );

  /// Write Nifti file.
  static void WriteNifti( const char* pathImg, const UniformVolume& volume );

  /// Write MetaImage file.
  static void WriteMetaImage( const char* pathHdr, const UniformVolume& volume );

#ifdef CMTK_BUILD_NRRD
  /// Read NRRD file.
  static const UniformVolume::SmartPtr ReadNRRD( const char* path );

  /// Write NRRD file.
  static void WriteNRRD( const char* pathHdr, const UniformVolume& volume );
#else
  /// Read NRRD file.
  static const UniformVolume::SmartPtr ReadNRRD( const char* ) 
  {
    throw Exception( "Library was configured without Nrrd support." );
  }

  /// Write NRRD file.
  static void WriteNRRD( const char*, const UniformVolume&, const bool = false )
  {
    throw Exception( "Library was configured without Nrrd support." );
  }
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVolumeFromFile_h_included_
