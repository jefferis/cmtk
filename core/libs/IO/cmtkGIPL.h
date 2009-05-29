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

#ifndef __cmtkGIPL_h_included_
#define __cmtkGIPL_h_included_

#include <cmtkconfig.h>

#include <string.h>

#define GIPL_HEADERSIZE 256

#define GIPL_MAGIC 0x2AE389B8

// GIPL primitive data types
#define GIPL_BINARY       1
#define GIPL_CHAR         7
#define GIPL_U_CHAR       8
#define GIPL_SHORT        15
#define GIPL_U_SHORT      16
#define GIPL_U_INT        31
#define GIPL_INT          32
#define GIPL_FLOAT        64
#define GIPL_DOUBLE       65
#define GIPL_C_SHORT      144
#define GIPL_C_INT        160
#define GIPL_C_FLOAT      192
#define GIPL_C_DOUBLE     193

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Interface to GIPL image file header data structure.
 */
class FileHeaderGIPL 
{
public:
  /// Get image dimensions (numbers of voxels per dimension).
  const unsigned short int* GetDims() const 
  {
    static unsigned short int dims[4];
    memcpy( dims, Data, sizeof( dims ) );
#ifndef WORDS_BIGENDIAN
    SwapBytes( dims, 4, sizeof( dims[0] ) );
#endif
    return dims;
  }

  /// Get primitive image data type ID.
  short GetType() const 
  {
    short type;
    memcpy( &type, Data+8, sizeof( type ) );
#ifndef WORDS_BIGENDIAN
    SwapBytes( &type, 1, sizeof( type ) );
#endif
    return type;
  }

  /// Get voxel size in mm.
  const float* GetSize() const 
  {
    static float size[4];
    memcpy( size, Data+10, sizeof( size ) );
#ifndef WORDS_BIGENDIAN
    SwapBytes( size, 4, sizeof( size[0] ) );
#endif
    return size;
  }

  /// Get magic number.
  unsigned int GetMagic() const 
  {
    unsigned int magic;
    memcpy( &magic, Data+252, sizeof( magic ) );
#ifndef WORDS_BIGENDIAN
    SwapBytes( &magic, 1, sizeof( magic ) );
#endif
    return magic;
  }

  /// Get pointer to internal data array.
  char *GetData() { return Data; }

  /// Set internal data array.
  void *SetData( void *data, const size_t size ) 
  { 
    memcpy( Data, data, (size > GIPL_HEADERSIZE) ? GIPL_HEADERSIZE : size );
    return data;
  }

  /// Get size of internal data array.
  static int GetDataSize() { return GIPL_HEADERSIZE; }
  
private:
  /// The 256 byte header of the current file.
  char Data[GIPL_HEADERSIZE];

  /// Swap bytes of a buffer.
  static void SwapBytes( void *const buffer, const size_t items, const size_t item_size ) 
  {
    char* c_buf = static_cast<char*>( buffer );
    for ( size_t i=0; i<items; ++i ) {
      for ( size_t j=0; j<item_size/2; ++j ) {
	char d = c_buf[item_size-j];
	c_buf[item_size-j] = c_buf[j];
	c_buf[j] = d;
      }
      c_buf += item_size;
    }
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkGIPL_h_included_
