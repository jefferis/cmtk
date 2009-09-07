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

#include <cmtkVolumeFromFile.h>

#include <cmtkCompressedStream.h>
#include <cmtkTypes.h>
#include <cmtkConsole.h>

#include <stdio.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{ 

/// BioRad microscopy image file header
typedef struct
{
  unsigned short nx, ny;    //  0   2*2     image width and height in pixels
  short npic;               //  4   2       number of images in file
  short ramp1_min;          //  6   2*2     LUT1 ramp min. and max.
  short ramp1_max;
  int notes;                // 10   4       no notes=0; has notes=non zero
  short byte_format;        // 14   2       bytes=TRUE(1); words=FALSE(0)
  unsigned short n;         // 16   2       image number within file
  char name[32];            // 18   32      file name
  short merged;             // 50   2       merged format
  unsigned short color1;    // 52   2       LUT1 color status
  unsigned short file_id;   // 54   2       valid .PIC file=12345
  short ramp2_min;          // 56   2*2     LUT2 ramp min. and max.
  short ramp2_max;
  unsigned short color2;    // 60   2       LUT2 color status
  short edited;             // 62   2       image has been edited=TRUE(1)
  short lens;               // 64   2       Integer part of lens magnification
  float mag_factor;         // 66   4       4 byte real mag. factor (old ver.)
  unsigned short dummy[3];  // 70   6       NOT USED (old ver.=real lens mag.)
} FileHeaderBioRad;

UniformVolume*
VolumeFromFile::ReadBioRad( const char* path )
{
  CompressedStream stream( path );

  // Biorad header is 76 bytes
  char buffer[76];
  if ( 1 != stream.Read( &buffer, sizeof(buffer), 1 ) ) 
    {
    StdErr << "ERROR: cannot read header from BioRad file " << path << ". Bailing out.\n";
    return NULL;
    }
  
  FileHeaderBioRad header;
  memcpy( &header.nx, buffer+0, sizeof( header.nx ) );
  memcpy( &header.ny, buffer+2, sizeof( header.ny ) );
  memcpy( &header.npic, buffer+4, sizeof( header.npic ) );

  memcpy( &header.notes, buffer+10, sizeof( header.notes ) );
  memcpy( &header.byte_format, buffer+14, sizeof( header.byte_format ) );
  memcpy( &header.file_id, buffer+54, sizeof( header.file_id ) );

  // check MagicNumber
#ifdef WORDS_BIGENDIAN    
  if ( Memory::ByteSwap( header.file_id ) != 12345 ) 
    {
    StdErr << "ERROR: BioRad file " << path << " has invalid magic number. Bailing out.\n";
    return NULL;
    }
#else
  if ( header.file_id != 12345 ) 
    {
    StdErr << "ERROR: BioRad file " << path << " has invalid magic number. Bailing out.\n";
    return NULL;
    }
#endif

#ifdef WORDS_BIGENDIAN    
  int dims[3] = { Memory::ByteSwap( header.nx ), Memory::ByteSwap( header.ny ), Memory::ByteSwap( header.npic ) };
#else
  int dims[3] = { header.nx, header.ny, header.npic };
#endif
  int numPixels = dims[0] * dims[1] * dims[2];
  
  TypedArray::SmartPtr dataArray;
  if ( header.byte_format ) 
    {
    dataArray = TypedArray::SmartPtr( TypedArray::Create( TYPE_BYTE, numPixels ) );
    } 
  else
    {
    dataArray = TypedArray::SmartPtr( TypedArray::Create( TYPE_USHORT, numPixels ) );    
    }
  
  stream.Read( dataArray->GetDataPtr(), dataArray->GetItemSize(), dataArray->GetDataSize() );
  
  double pixelsizeX = 1, pixelsizeY = 1, pixelsizeZ = 1;
  bool flipX = false, flipY = false, flipZ = false;;
  
  while ( ! stream.Feof() ) 
    {
    char lineheader[16], line[80];
    stream.Read( lineheader, sizeof( lineheader ), 1 );
    stream.Read( line, sizeof( line ), 1 );
    
    //      StdErr.printf( "%s\n->%s\n", lineheader, line );
    
    double d1, d2, d3;
    if ( 3 == sscanf( line, "AXIS_2 %lf %lf %lf", &d1, &d2, &d3 ) ) 
      {
      pixelsizeX = fabs( d3 );
      flipX = (d3 < 0 );
      }
    if ( 3 == sscanf( line, "AXIS_3 %lf %lf %lf", &d1, &d2, &d3 ) ) 
      {
      pixelsizeY = fabs( d3 );
      flipY = (d3 < 0 );
      }
    if ( 3 == sscanf( line, "AXIS_4 %lf %lf %lf", &d1, &d2, &d3 ) ) 
      {
      pixelsizeZ = fabs( d3 );
      flipZ = (d3 < 0 );
      }
    }
  
#ifdef WORDS_BIGENDIAN
  Types::Coordinate lensScale = 1; //Memory::ByteSwap( header.lens );
#else
  Types::Coordinate lensScale = 1; //header.lens;
#endif
  
// GJ I think this was backwards - want to swap if we ARE on big endian
// since Biorad is always little endian
//#ifndef WORDS_BIGENDIAN
#ifdef WORDS_BIGENDIAN
  // change endianness from Sun to whatever we're currently on.
  dataArray->ChangeEndianness();
#endif
  
  const Types::Coordinate volSize[3] = { (dims[0] - 1) * lensScale * pixelsizeX, (dims[1] - 1) * lensScale * pixelsizeY, (dims[2] - 1) * pixelsizeZ };
  
  UniformVolume* volume = new UniformVolume( dims, volSize, dataArray );

  if ( flipX )
    {
    StdErr << "WARNING: x pixel spacing is negative. Resulting volume will be mirrored accordingly.\n";
    volume->ApplyMirrorPlane( AXIS_X );
    }
  if ( flipY )
    {
    StdErr << "WARNING: y pixel spacing is negative. Resulting volume will be mirrored accordingly.\n";
    volume->ApplyMirrorPlane( AXIS_Y );
    }
  if ( flipZ )
    {
    StdErr << "WARNING: z pixel spacing is negative. Resulting volume will be mirrored accordingly.\n";
    volume->ApplyMirrorPlane( AXIS_Z );
    }

  return volume;
}

} // namespace cmtk
