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

#include <cmtkXformIO.h>

#include <cmtkConsole.h>

#include <cmtkMountPoints.h>
#include <cmtkFileFormat.h>
#include <cmtkClassStream.h>
#include <cmtkClassStreamAffineXform.h>
#include <cmtkTypedStreamStudylist.h>
#include <cmtkFileUtil.h>
#include <cmtkAffineXformITKIO.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

Xform* 
XformIO::Read( const char* path, const bool verbose )
{
  const char* realPath = MountPoints::Translate( path );
  
  switch ( FileFormat::Identify( realPath ) ) 
    {
    case FILEFORMAT_NRRD: 
    {
#ifdef CMTK_BUILD_NRRD
    return Self::ReadNrrd( realPath, verbose );
#else
    StdErr << "ERROR: " << realPath << " is a Nrrd file, but Nrrd support is not enabled.\n"
	      << "  Please re-configure software using either '--with-nrrd' or '--with-nrrd-teem' switch.\n";
    return NULL;
#endif
    }
    case FILEFORMAT_ITK_TFM:
      return AffineXformITKIO::Read( path );
      break;
    case FILEFORMAT_STUDYLIST: 
    {
    if ( verbose ) 
      {
      StdErr << "Reading transformation from studylist " << realPath << "\n";
      }
    
    TypedStreamStudylist studylist( realPath );
    if ( studylist.GetWarpXform() )
      return studylist.GetWarpXform().ReleasePtr();
    else
      return studylist.GetAffineXform()->MakeInverse();
    }
    case FILEFORMAT_TYPEDSTREAM: 
    {
    if ( verbose ) 
      {
      StdErr << "Reading transformation from typedstream file " << realPath << "\n";
      }
    
    ClassStream stream( realPath, ClassStream::READ );
    WarpXform* warpXform;
    stream >> warpXform;
    
    if ( warpXform ) return warpXform;
    
    stream.Open( realPath, ClassStream::READ );
    try
      {
      AffineXform affineXform;
      stream >> affineXform;
      return new AffineXform( affineXform );
      }
    catch (...)
      {
      return NULL;
      }
    }
    default:
      StdErr << "The file/directory " << realPath << " does not seem to be in a supported transformation format\n";
    }
  return NULL;
}

void 
XformIO::Write
( const Xform* xform, const char *path, const bool verbose )
{
  FileFormatID fileFormat = FILEFORMAT_TYPEDSTREAM;

  const char* suffix = strrchr( path, '.' );
  if ( suffix )
    {
    if ( ! strcmp( ".nrrd", suffix ) || ! strcmp( ".nhdr", suffix ) )
      {
      fileFormat = FILEFORMAT_NRRD;
      }
    else
      {
      if ( ! strcmp( ".tfm", suffix ) )
	{
	fileFormat = FILEFORMAT_ITK_TFM;
	}      
      }
    }
  
  char absolutePath[PATH_MAX];
  FileUtils::GetAbsolutePath( absolutePath, path );
  
  switch ( fileFormat )
    {
    case FILEFORMAT_NRRD:
#ifdef CMTK_BUILD_NRRD
      WriteNrrd( xform, absolutePath, verbose );
#else
      StdErr << "ERROR: XformIO::Write -- Nrrd support not configured.\n";
#endif
      break;
    case FILEFORMAT_ITK_TFM:
    {
    const AffineXform* affineXform = dynamic_cast<const AffineXform*>( xform );
    if ( affineXform )
      AffineXformITKIO::Write( path, affineXform );
    break;
    }
    case FILEFORMAT_TYPEDSTREAM:
    {
    ClassStream stream( absolutePath, ClassStream::WRITE );
    
    const AffineXform* affineXform = dynamic_cast<const AffineXform*>( xform );
    if ( affineXform )
      stream << *affineXform;
    
    const SplineWarpXform* splineWarpXform = dynamic_cast<const SplineWarpXform*>( xform );
    if ( splineWarpXform )
      stream << *splineWarpXform;
    }
    break;
    default:
      break;
    }
}

} // namespace cmtk
