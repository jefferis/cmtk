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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include "cmtkXformIO.h"

#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkFileUtils.h>
#include <System/cmtkMountPoints.h>
#include <System/cmtkExitException.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkClassStreamAffineXform.h>
#include <IO/cmtkClassStreamPolynomialXform.h>
#include <IO/cmtkTypedStreamStudylist.h>
#include <IO/cmtkAffineXformITKIO.h>

#include <string>

namespace
cmtk
{

/** \addtogroup IO */
//@{

Xform::SmartPtr 
XformIO::Read( const std::string& path )
{
  const std::string realPath = MountPoints::Translate( path );
  
  switch ( FileFormat::Identify( realPath ) ) 
    {
    case FILEFORMAT_NRRD: 
#ifdef CMTK_BUILD_NRRD
      return Self::ReadNrrd( realPath );
#else
      StdErr << "ERROR: " << realPath << " is a Nrrd file, but Nrrd support is not enabled.\n"
	     << "  Please re-configure software using either '--with-nrrd' or '--with-nrrd-teem' switch.\n";
      return Xform::SmartPtr( NULL );
#endif
    case FILEFORMAT_ITK_TFM:
      return AffineXformITKIO::Read( path );
    case FILEFORMAT_STUDYLIST: 
      DebugOutput( 1 ) << "Reading transformation from studylist " << realPath << "\n";
      {
      TypedStreamStudylist studylist( realPath );
      if ( studylist.GetWarpXform() )
	return studylist.GetWarpXform();
      else
	return studylist.GetAffineXform();
      }
    case FILEFORMAT_TYPEDSTREAM: 
      DebugOutput( 1 ) << "Reading transformation from typedstream file " << realPath << "\n";
      {
      ClassStreamInput stream( realPath );
      WarpXform* warpXform;
      stream >> warpXform;
      
      if ( warpXform ) 
	return Xform::SmartPtr( warpXform );
      
      stream.Open( realPath );
      PolynomialXform polyXform;
      try 
	{
	stream >> polyXform;
	return Xform::SmartPtr( new PolynomialXform( polyXform ) );
	}
      catch ( const cmtk::Exception& )
	{
	}
      
      stream.Open( realPath );
      AffineXform affineXform;
      try 
	{
	stream >> affineXform;
	}
      catch ( const cmtk::Exception& ex )
	{
	StdErr << "ERROR: " << ex.what() << "\n";
	return Xform::SmartPtr( NULL );
	}
      return Xform::SmartPtr( new AffineXform( affineXform ) );
      }
    case FILEFORMAT_NEXIST:
      StdErr << "The file/directory " << realPath << " does not exist or cannot be read\n";
      throw ExitException( 1 );
    default:
      StdErr << "The file/directory " << realPath << " does not seem to be in a supported transformation format\n";
      throw ExitException( 1 );
    }
  return Xform::SmartPtr( NULL );
}

void 
XformIO::Write
( const Xform* xform, const std::string& path )
{
  FileFormatID fileFormat = FILEFORMAT_TYPEDSTREAM;

  const size_t period = path.rfind( '.' );
  if ( period != std::string::npos )
    {
    const std::string suffix = path.substr( period );
    if ( (suffix == ".nrrd") || (suffix == ".nhdr") )
      {
      fileFormat = FILEFORMAT_NRRD;
      }
    else if ( suffix == ".nii")
      {
      fileFormat = FILEFORMAT_NIFTI_SINGLEFILE;
      }
    else if ( suffix == ".img" )
      {
      fileFormat = FILEFORMAT_NIFTI_DETACHED;
      }
    else
      {
      if ( (suffix == ".tfm") || (suffix == ".txt") )
	{
	fileFormat = FILEFORMAT_ITK_TFM;
	}      
      }
    }
  
  const std::string absolutePath = FileUtils::GetAbsolutePath( path );
  
  switch ( fileFormat )
    {
    case FILEFORMAT_NRRD:
#ifdef CMTK_BUILD_NRRD
      WriteNrrd( xform, absolutePath );
#else
      StdErr << "ERROR: XformIO::Write -- Nrrd support not configured.\n";
#endif
      break;
    case FILEFORMAT_NIFTI_DETACHED:
    case FILEFORMAT_NIFTI_SINGLEFILE:
      Self::WriteNIFTI( xform, absolutePath );
      break;
    case FILEFORMAT_ITK_TFM:
    {
    const AffineXform* affineXform = dynamic_cast<const AffineXform*>( xform );
    if ( affineXform )
      AffineXformITKIO::Write( path, *affineXform );
    break;
    }
    case FILEFORMAT_TYPEDSTREAM:
    {
    ClassStreamOutput stream( absolutePath, ClassStreamOutput::MODE_WRITE );
    
    const AffineXform* affineXform = dynamic_cast<const AffineXform*>( xform );
    if ( affineXform )
      stream << *affineXform;
    
    const PolynomialXform* polyXform = dynamic_cast<const PolynomialXform*>( xform );
    if ( polyXform )
      stream << *polyXform;
    
    const SplineWarpXform* splineWarpXform = dynamic_cast<const SplineWarpXform*>( xform );
    if ( splineWarpXform )
      stream << *splineWarpXform;
    }
    break;
    default:
      // cannot really get here, but gcc doesn't know that
      break;
    }
}

} // namespace cmtk
