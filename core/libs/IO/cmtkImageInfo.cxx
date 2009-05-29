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

#include <cmtkImageInfo.h>

#include <stdlib.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

ImageInfo::ImageInfo ()
{
  imagepath = NULL;
  SetImagePosition( 0, 0, 0 );
  SetImageOrientation( 1, 0, 0, 0, 1, 0 );
  dims[0] = dims[1] = dims[2] = 0;
  signbit = 0;
  bytesperpixel = 0;
  Padding = false;
  custom = false;
}

ImageInfo::~ImageInfo () 
{
  if ( imagepath ) free(imagepath);
}

ImageInfo* ImageInfo::Clone() const
{
  ImageInfo *instance = new ImageInfo;

  instance->imagepath = strdup( imagepath );

  for ( int dim = 0; dim<3; ++dim )
    instance->dims[dim] = dims[dim];
  instance->offset = offset;
  instance->bytesperpixel = bytesperpixel;
  instance->swapbytes = swapbytes;
  instance->signbit = signbit;
  instance->datatype = datatype;

  instance->minimum = minimum;
  instance->maximum = maximum;
  instance->black = black;
  instance->white = white;

  instance->custom = custom;
  instance->slicedistance = slicedistance;
  instance->calibrationx = calibrationx;
  instance->calibrationy = calibrationy;
  instance->original_calibrationx = original_calibrationx;
  instance->original_calibrationy = original_calibrationy;
  instance->tablepos = tablepos;
  instance->original_tablepos = original_tablepos;

  instance->ImagePosition = ImagePosition;
  instance->ImageOrientation[0] = ImageOrientation[0];
  instance->ImageOrientation[1] = ImageOrientation[1];
  
  return instance;
}

const char*
ImageInfo::ImageSuffix ( const int format ) 
{
  switch (format) {
  case 0: return "ima"; // ACR-NEMA
  case 1: return "ima"; // DICOM 3.0
  case 2: return "pgm"; // portable greymap
  case 3: return "raw"; // raw data
  }
  return "";
}

void ImageInfo::SetBytesPerPixel ( const int bpp )
{
  bytesperpixel = bpp;
  datatype = SelectDataTypeInteger( bytesperpixel, signbit );
}

void ImageInfo::SetImagePosition ( const double x, const double y, const double z )
{
  ImagePosition.Set( x, y, z );
}

void 
ImageInfo::SetImageOrientation
( const double xx, const double xy, const double xz, const double yx, const double yy, const double yz )
{
  ImageOrientation[0].Set( xx, xy, xz );
  ImageOrientation[1].Set( yx, yy, yz );
}

} // namespace cmtk
