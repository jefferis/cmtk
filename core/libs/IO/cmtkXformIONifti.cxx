/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
//  $Revision: 5129 $
//
//  $LastChangedDate: 2014-01-09 13:48:15 -0800 (Thu, 09 Jan 2014) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include "cmtkXformIO.h"

#include <Base/cmtkDeformationField.h>

#include <IO/cmtkVolumeIO.h>

#include <System/cmtkCompressedStream.h>

#include "nifti1.h"
#include "nifti1_io_math.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef HAVE_ZLIB
#  include <zlib.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

void 
XformIO::WriteNIFTI
( const Xform* xform, const std::string& path )
{
  const DeformationField* dfield = dynamic_cast<const DeformationField*>( xform );
  if ( ! dfield )
    {
    StdErr << "ERROR: XformIO::WriteNIFTI can only write DeformationField objects so far.\n"
	   << "       No data was written.\n";
    return;
    }
  
  const size_t dfieldRegionSize = dfield->m_Dims.Product();
  TypedArray::SmartPtr data = TypedArray::Create( TYPE_COORDINATE, 3 * dfieldRegionSize );
  for ( size_t ofs = 0; ofs < dfieldRegionSize; ++ofs )
    {
    for ( int dim = 0; dim < 3; ++dim )
      {
      data->Set( dfield->m_Parameters[3*ofs+dim], ofs + dim * dfieldRegionSize );
      }
    }
  
  //    if ( dfield->MetaKeyExists(META_SPACE_UNITS_STRING) )
//      nval->spaceUnits[0] = strdup( dfield->GetMetaInfo( META_SPACE_UNITS_STRING ).c_str() );
      
  bool detachedHeader = false;
  bool forceCompressed = false;

  std::string pathImg( path );

  // first, look for .gz
  size_t suffixPosGz = pathImg.rfind( std::string( ".gz" ) );
  if ( suffixPosGz != std::string::npos )
    {
    // found: set force compression flag and remove .gz from path
    forceCompressed = true;
    pathImg = pathImg.substr( 0, suffixPosGz );
    }
  
  std::string pathHdr( pathImg );
  size_t suffixPos = pathHdr.rfind( ".img" );
  if ( suffixPos != std::string::npos )
    {
    detachedHeader = true;
    pathHdr.replace( suffixPos, 4, ".hdr" );
    }
  
  nifti_1_header header;
  memset( &header, 0, sizeof( header ) );

  header.sizeof_hdr = 348; // header size
  header.dim_info = 0;

  // ndims
  header.dim[0] = 5;

  // dimensions
  header.dim[1] = dfield->m_Dims[AXIS_X];
  header.dim[2] = dfield->m_Dims[AXIS_Y];
  header.dim[3] = dfield->m_Dims[AXIS_Z];
  header.dim[4] = 1;
  header.dim[5] = 3;
  header.dim[6] = 0;
  header.dim[7] = 0;

  header.pixdim[0] = 1.0;
  header.pixdim[1] = static_cast<float>( dfield->m_Spacing[AXIS_X] );
  header.pixdim[2] = static_cast<float>( dfield->m_Spacing[AXIS_Y] );
  header.pixdim[3] = static_cast<float>( dfield->m_Spacing[AXIS_Z] );
  header.pixdim[4] = 0.0;
  header.pixdim[5] = 1.0;
  
  header.intent_code = NIFTI_INTENT_DISPVECT;

  header.qform_code = header.sform_code = 0;

  header.bitpix = 8 * sizeof( Types::Coordinate );
  if ( sizeof( Types::Coordinate ) == sizeof( float ) )
    {
    header.datatype = DT_FLOAT;
    }
  else
    {
    header.datatype = DT_DOUBLE;
    }  
  
  // determine data range;
  const Types::DataItemRange dataRange = data->GetRange();
  header.cal_max = static_cast<float>( dataRange.m_UpperBound );
  header.cal_min = static_cast<float>( dataRange.m_LowerBound );

#ifdef _MSC_VER
  const char *const modestr = "wb";
#else
  const char *const modestr = "w";
#endif
    
  if ( detachedHeader )
    {
    memcpy( &header.magic, "ni1\x00", 4 );
    header.vox_offset = 0;
    FILE *hdrFile = fopen( pathHdr.c_str(), modestr );
    if ( hdrFile ) 
      {
      fwrite( &header, 1, sizeof( header ), hdrFile );
      const int extension = 0;
      fwrite( &extension, 1, 4, hdrFile );
      fclose( hdrFile );
      }
    else
      {
      StdErr << "ERROR: NIFTI header file '" << pathHdr << "' could not be opened for writing!\n";
      }
    }
  else
    {
    memcpy( &header.magic, "n+1\x00", 4 );
    header.vox_offset = 352;
    }
  
  if ( VolumeIO::GetWriteCompressed() || forceCompressed )
    {
    struct stat buf;
    if ( ! stat( pathImg.c_str(), &buf ) )
      {
      StdErr << "WARNING: NIFTI file '" << path << "' will be written compressed, but uncompressed file exists!\n";
      }
    
    gzFile imgFile = gzopen( (pathImg+".gz").c_str(), modestr );
    if ( imgFile ) 
      {
      if ( ! detachedHeader )
	{
	gzwrite( imgFile, &header, sizeof( header ) );
	const int extension = 0;
	gzwrite( imgFile, &extension, 4 );
	}
      
      const size_t dataSize = data->GetItemSize() * data->GetDataSize();
      if ( dataSize != CompressedStream::Zlib::StaticSafeWrite( imgFile, data->GetDataPtr(), dataSize ) )
	{
	StdErr << "WARNING: gzwrite() returned error when writing to " << pathImg << "\n";
	}
      gzclose( imgFile );
      }
    else
      {
      StdErr << "ERROR: could not open file '" << pathImg << ".gz' for writing\n";
      }
    }
  else
    {
    FILE *imgFile = fopen( pathImg.c_str(), modestr );
    if ( imgFile ) 
      {
      if ( ! detachedHeader )
	{
	fwrite( &header, 1, sizeof( header ), imgFile );
	const int extension = 0;
	fwrite( &extension, 1, 4, imgFile );
	}

      fwrite( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize(), imgFile );
      fclose( imgFile );
      }
    else
      {
      StdErr << "ERROR: could not open file '" << pathImg << "' for writing\n";
      }
    }
}

//@}

} // namespace cmtk
