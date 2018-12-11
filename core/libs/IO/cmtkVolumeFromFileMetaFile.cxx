/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkVolumeFromFile.h"

#include <System/cmtkCompressedStream.h>
#include <System/cmtkConsole.h>
#include <Base/cmtkUniformVolume.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

void
VolumeFromFile::WriteMetaImage
( const std::string& pathHdr, const UniformVolume& volume )
{
  const TypedArray* data = volume.GetData();
  if ( ! data ) return;

#ifdef _MSC_VER
  FILE *outfile = fopen( pathHdr.c_str(), "wb" );
#else
  FILE *outfile = fopen( pathHdr.c_str(), "w" );
#endif

  if ( ! outfile )
    {
    StdErr << "Could not open file " << pathHdr << " for writing.\n";
    return;
    }

  fprintf( outfile, "ObjectType = Image\n" );
  fprintf( outfile, "NDims = 3\n" );
  fprintf( outfile, "BinaryData = True\n" );
#ifndef WORDS_BIGENDIAN
  fprintf( outfile, "BinaryDataByteOrderMSB = False\n" );
  fprintf( outfile, "ElementByteOrderMSB = False\n" );
#else
  fprintf( outfile, "BinaryDataByteOrderMSB = True\n" );
  fprintf( outfile, "ElementByteOrderMSB = True\n" );
#endif
  const AffineXform::MatrixType matrix = volume.GetImageToPhysicalMatrix();
  fprintf( outfile, "TransformMatrix = %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
	   (double)matrix[0][0], (double)matrix[0][1], (double)matrix[0][2], 
	   (double)matrix[1][0], (double)matrix[1][1], (double)matrix[1][2], 
	   (double)matrix[2][0], (double)matrix[2][1], (double)matrix[2][2] );
  fprintf( outfile, "Offset = %lf %lf %lf\n", (double)matrix[3][0], (double)matrix[3][1], (double)matrix[3][2] );
  fprintf( outfile, "CenterOfRotation = 0 0 0\n" );
  fprintf( outfile, "ElementSpacing = %f %f %f\n", volume.m_Delta[AXIS_X], volume.m_Delta[AXIS_Y], volume.m_Delta[AXIS_Z] );
  fprintf( outfile, "DimSize = %d %d %d\n", volume.m_Dims[AXIS_X], volume.m_Dims[AXIS_Y], volume.m_Dims[AXIS_Z] );
  fprintf( outfile, "AnatomicalOrientation = %s\n", volume.GetMetaInfo( META_SPACE ).c_str() );
  fprintf( outfile, "ElementNumberOfChannels = 1\n" );
  
  fputs( "ElementType = ", outfile ) ;
  switch ( data->GetType() )
    {
    case TYPE_BYTE:
      fputs( "MET_UCHAR\n", outfile );
      break;
    case TYPE_CHAR:
      fputs( "MET_CHAR\n", outfile );
      break;
    case TYPE_SHORT:
      fputs( "MET_SHORT\n", outfile );
      break;
    case TYPE_USHORT:
      fputs( "MET_USHORT\n", outfile );
      break;
    case TYPE_INT:
      fputs( "MET_INT\n", outfile );
      break;
    case TYPE_UINT:
      fputs( "MET_UINT\n", outfile );
      break;
    case TYPE_FLOAT:
      fputs( "MET_FLOAT\n", outfile );
      break;
    case TYPE_DOUBLE:
      fputs( "MET_DOUBLE\n", outfile );
      break;
    default:
      fputs( "MET_UNKNOWN\n", outfile );
      break;      
    }
  fprintf( outfile, "ElementDataFile = LOCAL\n" );

  fwrite( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize(), outfile );
  fclose( outfile );
}

} // namespace cmtk
