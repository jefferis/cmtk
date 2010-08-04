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

#include "cmtkPGM.h"

#include "System/cmtkCompressedStream.h"
#include "System/cmtkConsole.h"
#include "IO/cmtkImageInfo.h"

#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

#ifndef isspace
#define isspace(c)      ((c=='\t') || (c==' ') || (c=='\n'))
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

ScalarImage* 
PGM::Read( const char* filename ) 
{
  CompressedStream stream(filename);
  if ( ! stream.IsValid() ) 
    {
    StdErr.printf( "File %s could not be opened.", filename );
    return NULL;
    }
  
  Types::Coordinate pixelSize[2] = { 1, 1 };
  
  int i;
  char c, fileID[3], line[1024];
  stream.Get(fileID[0]); 
  stream.Get(fileID[1]);
  stream.Get(c);
  
  do
    {
    stream.Get(line[0]);
    if (line[0]=='#') 
      {
      int idx = 1;
      c=0;
      while (c != '\n') { stream.Get(c); line[idx++] = c; }
      }
    float tmpf0, tmpf1;
    if ( sscanf( line, "# calibration %f %f", &tmpf0, &tmpf1 ) == 2 ) 
      {
      pixelSize[0] = tmpf0;
      pixelSize[1] = tmpf1;
      } 
    else
      {
      if ( sscanf( line, "# tablepos %f", &tmpf0 ) == 1 ) 
	{
	}
      }
    } while (line[0]=='#');
  
  i = 1;
  for ( int spaces=0; spaces<3; ++i ) 
    {
    stream.Get(line[i]);
    if (isspace(line[i]))
      ++spaces;
    }
  line[i]=0;
  
  unsigned int dimsx, dimsy, maxvalue;
  sscanf( line, "%d%d%d", &dimsx , &dimsy, &maxvalue );
  
  int bytesperpixel = 1;
  while (maxvalue > 255) 
    {
    ++bytesperpixel;
    maxvalue /= 256;
    }
  
  int dim = dimsx * dimsy;
  
  TypedArray::SmartPtr pixelData;
  switch ( bytesperpixel ) 
    {
    case 1:
      pixelData = TypedArray::Create( TYPE_BYTE, dim );
      break;
    case 2:
      pixelData = TypedArray::Create( TYPE_USHORT, dim );
      break;
    case 4:
      pixelData = TypedArray::Create( TYPE_INT, dim );
      break;
    default:
      return NULL;
    }
  stream.Read( pixelData->GetDataPtr(), bytesperpixel, dim);
  
  ScalarImage *image = new ScalarImage( dimsx, dimsy );
  image->SetPixelSize( pixelSize );
  image->SetPixelData( TypedArray::SmartPtr( pixelData ) );
  
  return image;
}

void
PGM::Write
( const char* filename, const ScalarImage *image, const Types::DataItem greyFrom, const Types::DataItem greyTo )
{
  const size_t numberOfPixels = image->GetNumberOfPixels();
  byte *pgmData = Memory::AllocateArray<byte>(  numberOfPixels  );

  const TypedArray *pixelData = image->GetPixelData();

  const Types::DataItem greyScale = 255.0 / (greyTo - greyFrom);
  
  for ( unsigned int i = 0; i < numberOfPixels; ++i ) 
    {
    Types::DataItem pixel;
    if ( pixelData->Get( pixel, i ) ) 
      {
      if ( pixel <= greyFrom )
	pgmData[i] = 0;
      else
	if ( pixel >= greyTo )
	  pgmData[i] = 255;
	else
	  pgmData[i] = static_cast<byte>( (pixel - greyFrom) * greyScale );
      }
    }
  
  FILE *fp = fopen( filename, "wb" );
  if ( fp ) 
    {
    fprintf( fp, "P5\n" );
    fprintf( fp, "# calibration %f %f\n", image->GetPixelSize()[0], image->GetPixelSize()[1] );
    fprintf( fp, "# tablepos %f \n", image->GetImageSlicePosition() );
    
    fprintf( fp, "%d %d %d\n", image->GetDims()[0], image->GetDims()[1], 255 );
    
    fwrite( pgmData, 1, numberOfPixels, fp );
    
    fclose(fp);
    }

  Memory::DeleteArray( pgmData );
}

void
PGM::Write16bit( const char* filename, const ScalarImage *image, const Types::DataItem greyFrom, const Types::DataItem greyTo )
{
  const size_t numberOfPixels = image->GetNumberOfPixels();
  
  const TypedArray *pixelData = image->GetPixelData();

  const Types::DataItem greyScale = 255.0 / (greyTo - greyFrom);
  
  unsigned short *pgmData = Memory::AllocateArray<unsigned short>(  numberOfPixels  );
  unsigned short maxData = 0;
  for ( size_t i = 0; i < numberOfPixels; ++i ) 
    {
    Types::DataItem pixel;
    if ( pixelData->Get( pixel, i ) ) 
      {
      if ( pixel <= greyFrom )
	pixel = 0;
      else
	if ( pixel >= greyTo )
	  pixel = 65535;
	else
	  pixel = (pixel - greyFrom) * greyScale;
      
      // gthumb and ImageMagick want 16bit pgm in little endian, so let's do it if necessary
#ifdef WORDS_BIGENDIAN
      const unsigned short tmp = static_cast<unsigned short>( pixel );
      pgmData[i] = ((tmp&255)<<8) + (tmp>>8);
#else
      pgmData[i] = static_cast<unsigned short>( pixel );
#endif
      }
    else
      {
      pgmData[i] = 0;
      }
    maxData = std::max( maxData, pgmData[i] );
    }
  
  FILE *fp = fopen( filename, "wb" );
  if ( fp ) 
    {
    fprintf( fp, "P5\n" );
    fprintf( fp, "# calibration %f %f\n", image->GetPixelSize()[0], image->GetPixelSize()[1] );
    fprintf( fp, "# tablepos %f \n", image->GetImageSlicePosition() );
    
    fprintf( fp, "%d %d %d\n", image->GetDims()[0], image->GetDims()[1], maxData );
    
    fwrite( pgmData, sizeof( *pgmData ), numberOfPixels, fp );
    fclose(fp);
    }
  
  Memory::DeleteArray( pgmData );
}

} // namespace cmtk
