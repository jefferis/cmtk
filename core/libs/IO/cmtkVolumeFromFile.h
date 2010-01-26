/*
//
//  Copyright 2004-2010 SRI International
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

#include <cmtkUniformVolume.h>
#include <cmtkException.h>

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
  static UniformVolume* Read( const char *filename );

  /// Read volume in multi-slice DICOM format.
  static UniformVolume* ReadDICOM( const char *filename );

  /// Read volume in Vanderbilt format.
  static UniformVolume* ReadVanderbilt( const char *filename );

  /// Read BioRad PIC file (confocal microscopy image).
  static UniformVolume* ReadBioRad( const char* path );

  /// Read Analyze 7.5 file (separate header file).
  static UniformVolume* ReadAnalyzeHdr( const char* pathHdr, const bool bigEndian = false, const bool readData = true );

  /// Read Nifti file.
  static UniformVolume* ReadNifti( const char* pathHdr, const bool detached, const bool readData = true );

  /// Write Analyze 7.5 file (separate header file).
  static void WriteAnalyzeHdr( const char* pathHdr, const UniformVolume* volume, const bool verbose = false );

  /// Write Nifti file.
  static void WriteNifti( const char* pathImg, const UniformVolume* volume, const bool verbose = false );

  /// Write MetaImage file.
  static void WriteMetaImage( const char* pathHdr, const UniformVolume* volume );

#ifdef CMTK_BUILD_NRRD
  /// Read NRRD file.
  static UniformVolume* ReadNRRD( const char* path );

  /// Write NRRD file.
  static void WriteNRRD( const char* pathHdr, const UniformVolume* volume, const bool verbose = false );
#else
  /// Read NRRD file.
  static UniformVolume* ReadNRRD( const char* path ) 
  {
    throw Exception( "Library was configured without Nrrd support." );
  }

  /// Write NRRD file.
  static void WriteNRRD( const char* pathHdr, const UniformVolume* volume, const bool = false )
  {
    throw Exception( "Library was configured without Nrrd support." );
  }
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVolumeFromFile_h_included_
