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

#include <cmtkImageIO.h>

#ifdef CMTK_HAVE_DCMTK
#  include <cmtkDICOM.h>
#endif

#include <cmtkPGM.h>
#include <cmtkRAW.h>

#ifdef HAVE_STDARG_H
#  include <stdarg.h>
#endif

namespace 
cmtk
{

/** \addtogroup IO */
//@{

ImageIO*
ImageIO::Create( const char* format )
{
  if ( !strcmp( "ACR-NEMA", format ) )
    // Format is obsolete; use DICOM library instead.
#ifdef CMTK_HAVE_DCMTK
    return new DICOM;
#else
    return NULL;
#endif
  else if ( !strcmp( "DICOM", format ) )
#ifdef CMTK_HAVE_DCMTK
    return new DICOM;
#else
    return NULL;
#endif
  else if ( !strcmp( "PGM", format ) )
    return new PGM;
  else if ( !strcmp( "RAW-DATA", format ) || !strcmp( "RAW3D", format ) )
    return new RAW;

  return NULL;
}

ImageIO*
ImageIO::Create( const FileFormatID format )
{
  switch ( format ) 
    {
    case FILEFORMAT_DICOM:
#ifdef CMTK_HAVE_DCMTK
      return new DICOM;
#else
      return NULL;
#endif
    case FILEFORMAT_PGM:
      return new PGM;
    case FILEFORMAT_RAW:
    case FILEFORMAT_RAW3D:
      return new RAW;
    default:
      return NULL;
    }    
  return NULL;
}

ImageIO::ImageIO()
{
  Error = 0;
  Size = 0;
  DataPtr = 0;
  FreeDataPtrFlag = false;
}

size_t ImageIO::GetSize( const size_t itemSize ) const
{
  return ( Size / itemSize );
}

void* ImageIO::GetReleaseDataPtr()
{
  void *dataPtr = DataPtr;
  DataPtr = NULL;
  FreeDataPtrFlag = false;
  Size = 0;
  return dataPtr;
}

void ImageIO::SetErrorMsg( const char* format, ... )
{
  va_list args;
  va_start(args, format);
  vsnprintf( ErrorMsg, sizeof( ErrorMsg ), format, args );
  va_end(args);

  this->SetError();
}

const char* ImageIO::GetErrorMsg() const 
{ 
  if ( Error )
    return ErrorMsg; 
  else
    return NULL;
}

void
ImageIO::SetDataPtr
( void *const dataPtr, const size_t size, const size_t itemSize )
{
  this->FreeDataPtr();
  DataPtr = dataPtr;
  FreeDataPtrFlag = true;
  Size = size * itemSize;
}

void
ImageIO::LinkDataPtr
( void *const dataPtr, const size_t size, const size_t itemSize )
{
  this->FreeDataPtr();
  DataPtr = dataPtr;
  FreeDataPtrFlag = false;
  Size = size * itemSize;
}

void
ImageIO::FreeDataPtr()
{
  if ( FreeDataPtrFlag && DataPtr ) free( DataPtr );
  DataPtr = NULL;
}

void* 
ImageIO::AllocDataPtr
( const size_t size, const size_t itemSize )
{
  if ( DataPtr && ((size*itemSize) != Size) ) 
    {
    free( DataPtr );
    DataPtr = NULL;
    }

  if ( ! DataPtr ) 
    {
    Size = size * itemSize;
    DataPtr = malloc( Size );
    FreeDataPtrFlag = true;
    }

  return DataPtr;
}

} // namespace cmtk
