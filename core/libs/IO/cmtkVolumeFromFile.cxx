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

#include <cmtkVolumeFromFile.h>

#include <stdio.h>

#include <cmtkFileFormat.h>
#include <cmtkTypedArray.h>

#include <cmtkCompressedStream.h>
#include <cmtkDICOM.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

UniformVolume* 
VolumeFromFile::Read( const char *path )
{
  UniformVolume *volume = NULL;

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
      // for all other file formats, we shouldn't be in this function at all.
      return NULL;
    }
  
  return volume;
}

UniformVolume* 
VolumeFromFile::ReadDICOM( const char *path )
{
  ImageInfo imageInfo;
  StudyInfo studyInfo;
  DICOM dicomIO;

  dicomIO.Read( path, imageInfo, studyInfo, 0 );

  Types::Coordinate size[3] = {
    (imageInfo.dims[0]-1) * imageInfo.calibrationx,
    (imageInfo.dims[1]-1) * imageInfo.calibrationy,
    (imageInfo.dims[2]-1) * imageInfo.slicedistance
  };

  ScalarDataType dtype = SelectDataTypeInteger( imageInfo.bytesperpixel, imageInfo.signbit );

  TypedArray::SmartPtr dataArray
    ( TypedArray::Create( dtype, dicomIO.GetReleaseDataPtr(), imageInfo.dims[0]*imageInfo.dims[1]*imageInfo.dims[2], true /*freeArray*/, imageInfo.Padding, &imageInfo.PaddingValue ) );
		       
  UniformVolume *volume = new UniformVolume( imageInfo.dims, size, dataArray );

  return volume;
}

} // namespace cmtk
