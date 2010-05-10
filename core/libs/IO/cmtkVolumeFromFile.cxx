/*
//
//  Copyright 2004-2010 SRI International
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

#include <cmtkVolumeFromFile.h>

#include <stdio.h>

#include <cmtkFileFormat.h>
#include <cmtkTypedArray.h>

#include <cmtkCompressedStream.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

const UniformVolume::SmartPtr
VolumeFromFile::Read( const char *path )
{
  FileFormatID id = FileFormat::Identify( path );
  switch ( id )
    {
    case FILEFORMAT_DICOM:
      return VolumeFromFile::ReadDICOM( path );
    case FILEFORMAT_VANDERBILT:
      return VolumeFromFile::ReadVanderbilt( path );
    case FILEFORMAT_ANALYZE_HDR:
      return VolumeFromFile::ReadAnalyzeHdr( path, false /* bigendian */ );
    case FILEFORMAT_ANALYZE_HDR_BIGENDIAN:
      return VolumeFromFile::ReadAnalyzeHdr( path, true /* bigendian */ );
    default:
      ;
    }
  
  return UniformVolume::SmartPtr( NULL );
}

} // namespace cmtk
