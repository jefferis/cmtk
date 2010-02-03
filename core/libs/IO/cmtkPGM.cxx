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

#include <cmtkPGM.h>

#include <cmtkCompressedStream.h>
#include <cmtkImageInfo.h>

#include <stdlib.h>
#include <string.h>

#include <cmtkConsole.h>

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

#ifndef isspace
/// UNDOCUMENTED
#define isspace(c)      ((c=='\t') || (c==' ') || (c=='\n'))
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

ScalarImage* 
PGM::ReadPGM( const char* filename ) 
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
  
  TypedArray* pixelData;
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
  unsigned int numberOfPixels = image->GetNumberOfPixels();
  byte *pgmData = Memory::AllocateArray<byte>(  numberOfPixels  );

  const TypedArray *pixelData = image->GetPixelData();

  Types::DataItem greyScale = 255.0 / (greyTo - greyFrom);
  
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
  
  FILE *fp;
  if ( (fp = fopen( filename, "wb" ) ) == NULL ) 
    {
    return;
    }
  
  fprintf( fp, "P5\n" );
  fprintf( fp, "# calibration %f %f\n", image->GetPixelSize()[0], image->GetPixelSize()[1] );
  fprintf( fp, "# tablepos %f \n", image->GetImageSlicePosition() );

  fprintf( fp, "%d %d %d\n", image->GetDims()[0], image->GetDims()[1], 255 );
    
  fwrite( pgmData, 1, numberOfPixels, fp );

  fclose(fp);

  delete[] pgmData;
}

void
PGM::Write16bit( const char* filename, const ScalarImage *image )
{
  unsigned int numberOfPixels = image->GetNumberOfPixels();
  union
  {
    void *Void;
    byte *Byte;
    unsigned short *UShort;
  } pgmData;
  pgmData.Void = NULL;
  unsigned int bytesPerPixel = 0;

  const TypedArray *pixelData = image->GetPixelData();

  unsigned short maxData;
  Types::DataItem min, max;
  pixelData->GetRange( min, max );
  maxData = static_cast<unsigned short>( max );
  if ( maxData > 255 ) 
    {
    pgmData.UShort = Memory::AllocateArray<unsigned short>(  numberOfPixels  );
    bytesPerPixel = 2;
    for ( unsigned int i = 0; i < numberOfPixels; ++i ) 
      {
      Types::DataItem pixel;
      if ( pixelData->Get( pixel, i ) ) 
	{
	// pgm wants 16bit in little endian, so let's do it...
#ifdef WORDS_BIGENDIAN
	unsigned short tmp = static_cast<unsigned short>( pixel );
	pgmData.UShort[i] = ((tmp&255)<<8) + (tmp>>8);
#else
	pgmData.UShort[i] = static_cast<unsigned short>( pixel );
#endif
	}
      else
	{
	pgmData.UShort[i] = 0;
	}
      }
    } 
  else
    {
    pgmData.Byte = Memory::AllocateArray<byte>(  numberOfPixels  );
    bytesPerPixel = 1;
    for ( unsigned int i = 0; i < numberOfPixels; ++i ) 
      {
      Types::DataItem pixel;
      if ( pixelData->Get( pixel, i ) ) 
	{
	// pgm wants 16bit in little endian, so let's do it...
	pgmData.Byte[i] = static_cast<byte>( pixel );
	} 
      else
	{
	pgmData.Byte[i] = 0;
	}
      }
    }
  
  if ( pgmData.Void ) 
    {
    FILE *fp;
    if ( (fp = fopen( filename, "wb" ) ) == NULL ) 
      {
      return;
      }
    
    fprintf( fp, "P5\n" );
    fprintf( fp, "# calibration %f %f\n", image->GetPixelSize()[0], image->GetPixelSize()[1] );
    fprintf( fp, "# tablepos %f \n", image->GetImageSlicePosition() );
    
    fprintf( fp, "%d %d %d\n", image->GetDims()[0], image->GetDims()[1], maxData ? maxData : 1 );
    
    fwrite( pgmData.Void, bytesPerPixel, numberOfPixels, fp );
    
    fclose(fp);
    
    if ( bytesPerPixel == 1 ) 
      {
      delete[] pgmData.Byte;
      } 
    else
      {
      delete[] pgmData.UShort;
      }
    }
}

} // namespace cmtk
