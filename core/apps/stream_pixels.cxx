/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypedArray.h>

#include <IO/cmtkVolumeIO.h>

#include <iostream>
#include <vector>
#include <string>

int
doMain( int argc, const char *argv[] )
{
  const char* readOrientation = NULL;
  std::vector<std::string> imagePaths;
  cmtk::ScalarDataType convertToType = cmtk::TYPE_NONE;
  bool changeEndian = false;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Stream pixel data from one or more images to Standard Output." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool reads one or more images and writes all their pixels to standard output in binary form. Optionally, each image can be reoriented to a specified anatomical orientation "
		       "and/or converted to a different data type. This is useful for piping image data through a pipeline, e.g., the Camino DTI toolkit." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( "reorient" ), &readOrientation, "Reorient all images to anatomy-based standard orientation, e.g., RAS or LPS for axial. Default: no reorientation." );
    cl.AddSwitch( Key( "change-endian" ), &changeEndian, true, "Change endianness. Default: native endianness." );

    cmtk::CommandLine::EnumGroup<cmtk::ScalarDataType>::SmartPtr typeGroup = cl.AddEnum( "convert", &convertToType, "Convert data to new scalar type. Default: no conversion." );
    typeGroup->AddSwitch( Key( "char" ), cmtk::TYPE_CHAR, "8 bits, signed integer" );
    typeGroup->AddSwitch( Key( "byte" ), cmtk::TYPE_BYTE, "8 bits, unsigned integer" );
    typeGroup->AddSwitch( Key( "short" ), cmtk::TYPE_SHORT, "16 bits, signed integer" );
    typeGroup->AddSwitch( Key( "ushort" ), cmtk::TYPE_USHORT, "16 bits, unsigned integer" );
    typeGroup->AddSwitch( Key( "int" ), cmtk::TYPE_INT, "32 bits, signed integer" );
    typeGroup->AddSwitch( Key( "uint" ), cmtk::TYPE_UINT, "32 bits, unsigned integer" );
    typeGroup->AddSwitch( Key( "float" ), cmtk::TYPE_FLOAT, "32-bits, single-precision floating point" );
    typeGroup->AddSwitch( Key( "double" ), cmtk::TYPE_DOUBLE, "64 bits, double-precision floating point" );

    cl.AddParameterVector( &imagePaths, "ImagePaths", "List of paths to one or more image files. These will be read and processed in the given order." );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr volume;
  for ( size_t idx = 0; idx < imagePaths.size(); ++idx )
    {
    if ( readOrientation )
      volume = cmtk::VolumeIO::ReadOriented( imagePaths[idx], readOrientation );
    else
      volume = cmtk::VolumeIO::Read( imagePaths[idx] );

    if ( ! volume )
      {
      cmtk::StdErr << "ERROR: failed to read image " << imagePaths[idx] << "\n";
      throw cmtk::ExitException( 1 );
      }

    cmtk::TypedArray::SmartPtr dataArray = volume->GetData();
    if ( ! dataArray )
      {
      cmtk::StdErr << "ERROR: image " << imagePaths[idx] << " does not have any valid pixel data\n";
      throw cmtk::ExitException( 1 );
      }
    
    if ( convertToType != cmtk::TYPE_NONE )
      {
      dataArray = dataArray->Convert( convertToType );
      }

    if ( changeEndian )
      {
      dataArray->ChangeEndianness();
      }

    std::cout.write( reinterpret_cast<const char*>( dataArray->GetDataPtr() ), dataArray->GetDataSizeBytes() );
    }
  
  return 0;
}

#include "cmtkSafeMain"

