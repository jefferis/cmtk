/*
//
//  Copyright 2004-2011, 2013 SRI International
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

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <IO/cmtkAnalyze.h>
#include <IO/cmtkFileHeader.h>

#include <Base/cmtkTypes.h>

#include <stdio.h>

int DimsX = 0;
int DimsY = 0;
int DimsZ = 0;
bool PutDims = false;

int LegacyMode = 0;

void
SetDims( const char* arg )
{
  sscanf( arg, "%6d,%6d,%6d", &DimsX, &DimsY, &DimsZ );
  PutDims = true;
}

float DeltaX = 1.0;
float DeltaY = 1.0;
float DeltaZ = 1.0;
bool PutDeltas = false;

void
SetDeltas( const char* arg )
{
  sscanf( arg, "%15f,%15f,%15f", &DeltaX, &DeltaY, &DeltaZ );
  PutDeltas = true;
}

bool LittleEndian = false;
cmtk::ScalarDataType DataType = cmtk::TYPE_NONE;
cmtk::AnalyzeOrientation Orientation = cmtk::ANALYZE_UNKNOWN;

float Offset = 0;
bool PutOffset = false;

const char* HdrFileName = "header.nii";
const char* ImportHdrFile = NULL;

const char* Description = NULL;

int 
doMain( const int argc, const char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Make Analyze header file" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Make header file according Analzye 7.5 format based on user-supplied parameters for geometry, data type, orientation, etc." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "mk_analyze_hdr [options] [output.nii]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddCallback( Key( 'D', "dims" ), SetDims, "Set dimensions in voxels" );
    cl.AddCallback( Key( 'V', "voxel" ), SetDeltas, "Set voxel size in [mm]" );

    cl.AddOption( Key( 'O', "offset" ), &Offset, "Binary data file offset", &PutOffset );

    cl.AddSwitch( Key( 'c', "char" ), &DataType, cmtk::TYPE_CHAR, "8 bits, signed" );
    cl.AddSwitch( Key( 'b', "byte" ), &DataType, cmtk::TYPE_BYTE, "8 bits, unsigned" );
    cl.AddSwitch( Key( 's', "short" ), &DataType, cmtk::TYPE_SHORT, "16 bits, signed" );
    cl.AddSwitch( Key( 'u', "ushort" ), &DataType, cmtk::TYPE_USHORT, "16 bits, unsigned" );
    cl.AddSwitch( Key( 'i', "int" ), &DataType, cmtk::TYPE_INT, "32 bits signed" );
    cl.AddSwitch( Key( 'f', "float" ), &DataType, cmtk::TYPE_FLOAT, "32 bits floating point" );
    cl.AddSwitch( Key( 'd', "double" ), &DataType, cmtk::TYPE_DOUBLE, "64 bits floating point\n" );
    
    cl.AddSwitch( Key( 'B', "big-endian" ), &LittleEndian, false, "Big endian data" );
    cl.AddSwitch( Key( 'L', "little-endian" ), &LittleEndian, true, "Little endian data" );

    cl.AddSwitch( Key( "axial" ), &Orientation, cmtk::ANALYZE_AXIAL, "Set slice orientation to axial" );
    cl.AddSwitch( Key( "axial-flip" ), &Orientation, cmtk::ANALYZE_AXIAL_FLIP, "Set slice orientation to reverse axial" );
    cl.AddSwitch( Key( "sagittal" ), &Orientation, cmtk::ANALYZE_SAGITTAL, "Set slice orientation to sagittal" );
    cl.AddSwitch( Key( "sagittal-flip" ), &Orientation, cmtk::ANALYZE_SAGITTAL_FLIP, "Set slice orientation to reverse sagittal" );
    cl.AddSwitch( Key( "coronal" ), &Orientation, cmtk::ANALYZE_CORONAL, "Set slice orientation to corona;" );
    cl.AddSwitch( Key( "coronal-flip" ), &Orientation, cmtk::ANALYZE_CORONAL_FLIP, "Set slice orientation to reverse coronal" );

    cl.AddSwitch( Key( "legacy" ), &LegacyMode, 1, "Legacy mode; do not write 'SRI1' tag into header" );
    cl.AddSwitch( Key( "SRI1" ), &LegacyMode, -1, "Force 'SRI1' tag even if not present in imported header" );

    cl.AddOption( Key( 'I', "import" ), &ImportHdrFile, "Import data from given header file." );
    cl.AddOption( Key( "description" ), &Description, "Set description string [max. 80 characters]" );

    cl.Parse( argc, argv );

    HdrFileName = cl.GetNext();
  }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  char buffer[348];
  memset( buffer, 0, sizeof( buffer ) );

  if ( ImportHdrFile )
    {
    FILE *hdrIn = fopen( ImportHdrFile, "r" );
    if ( hdrIn )
      {
      fread( buffer, sizeof( buffer ), sizeof( *buffer ), hdrIn );
      fclose( hdrIn );

      cmtk::DebugOutput( 1 ) << "Imported header file " << ImportHdrFile << "\n";
      }
    else 
      {
      cmtk::StdErr << "ERROR: Could not open file " << ImportHdrFile << " for import.\n";
      throw cmtk::ExitException( 1 );
      }
    
    LittleEndian = (buffer[0] != 0x00);
    }
  
  cmtk::FileHeader header( buffer, !LittleEndian );
  // ID
  header.StoreField<int>( 0, 0x0000015C );

  // ndims
  if ( ! ImportHdrFile )
    {
    header.StoreField<short>( 40, 3 );
    }

  // dimensions
  if ( !ImportHdrFile || PutDims )
    {
    cmtk::DebugOutput( 1 ) << "Setting image dimensions\n";

    header.StoreField<short>( 42, DimsX );
    header.StoreField<short>( 44, DimsY );
    header.StoreField<short>( 46, DimsZ );
    header.StoreField<short>( 48, 0 ); // write dims 3 and 4
    header.StoreField<short>( 50, 0 ); // just for safety
    }

  if ( !ImportHdrFile && DataType == cmtk::TYPE_NONE )
    {
    cmtk::StdErr << "ERROR: you must either select a data type (e.g., byte, float) or import an existing header file.";
    return 1;
    }

  if ( !ImportHdrFile || DataType != cmtk::TYPE_NONE )
    {
    cmtk::DebugOutput( 1 ) << "Setting data type\n";
    
    switch ( DataType ) 
      {
      default:
	header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_NONE );
	header.StoreField<short>( 72, 0 );
	break;
      case cmtk::TYPE_BYTE:
	header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_UNSIGNED_CHAR );
	header.StoreField<short>( 72, 8 * sizeof( byte ) );
	break;
      case cmtk::TYPE_SHORT:
      case cmtk::TYPE_USHORT:
	header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_SIGNED_SHORT );
	header.StoreField<short>( 72, 8 * sizeof( short ) );
	break;
      case cmtk::TYPE_INT:
	header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_SIGNED_INT );
	header.StoreField<short>( 72, 8 * sizeof( signed int ) );
	break;
      case cmtk::TYPE_FLOAT:
	header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_FLOAT );
	header.StoreField<short>( 72, 8 * sizeof( float ) );
	break;
      case cmtk::TYPE_DOUBLE:
	header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_DOUBLE );
	header.StoreField<short>( 72, 8 * sizeof( double ) );
	break;
      }
    }

  if ( !ImportHdrFile || PutDeltas )
    {
    cmtk::DebugOutput( 1 ) << "Setting pixel size\n";

    header.StoreField<float>( 80, (float)DeltaX );
    header.StoreField<float>( 84, (float)DeltaY );
    header.StoreField<float>( 88, (float)DeltaZ );
    }
  
  if ( !ImportHdrFile )
    {
    header.StoreField<float>( 92, 1.0f ); // write sizes in dims 3 and
    header.StoreField<float>( 96, 1.0f ); // 4 just to be safe
    }

  // set offset for binary file.
  if ( !ImportHdrFile || PutOffset )
    {
    cmtk::DebugOutput( 1 ) << "Setting image file offset\n";
    header.StoreField<float>( 108, Offset );
    }

  // slice orientation
  if ( !ImportHdrFile && Orientation == cmtk::ANALYZE_UNKNOWN )
    {
    cmtk::StdErr << "ERROR: you must either select a slice orientation (e.g., axial, sagittal) or import an existing header file.";
    return 1;
    }

  if ( !ImportHdrFile || Orientation != cmtk::ANALYZE_UNKNOWN )
    {
    cmtk::DebugOutput( 1 ) << "Setting image orientation\n";
    header.StoreField<byte>( 252, Orientation );
    }

  if ( (!ImportHdrFile && LegacyMode == 0) || (LegacyMode < 0) )
    header.StoreFieldString( 344, "SRI1", 4 );
  else
    {
    // in "legacy" mode, remove "SRI1" tag even if present in imported header
    if ( LegacyMode > 0 )
      {
      header.StoreFieldString( 344, "\0x\0x\0x\0x", 4 );
      }
    }
  
  if ( Description )
    {
    cmtk::DebugOutput( 1 ) << "Setting image description\n";
    header.StoreFieldString( 148, Description, 80 );
    }

  // write header info
#ifdef _MSC_VER
  FILE *hdrFile = fopen( HdrFileName, "wb" );
#else
  FILE *hdrFile = fopen( HdrFileName, "w" );
#endif
  if ( hdrFile )
    {
    if ( 348 != fwrite( buffer, 1, 348, hdrFile ) ) 
      {
      cmtk::StdErr.printf( "ERROR: could not write 348 bytes to header file %s\n", HdrFileName );
      }
    fclose( hdrFile );
    }

  return 0;
}

#include "cmtkSafeMain"
