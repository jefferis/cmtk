/*
//
//  Copyright 2004-2011, 2013 SRI International
//
//  Copyright 1997-2011 Torsten Rohlfing
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
#include <System/cmtkCompressedStream.h>

#include <IO/cmtkAnalyze.h>
#include <IO/cmtkFileHeader.h>
#include <IO/nifti1.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkAffineXform.h>

int DimsX = 0;
int DimsY = 0;
int DimsZ = 0;
bool PutDims = false;

bool DetachedHeader = true;

void
SetDims( const char* arg )
{
  if ( 3 == sscanf( arg, "%6d,%6d,%6d", &DimsX, &DimsY, &DimsZ ) )
    {
    PutDims = true;
    }
  else
    {
    cmtk::StdErr << "ERROR: could not parse X,Y,Z dimensions from argument '" << arg << "'\n";
    throw cmtk::ExitException( 1 );
    }
}

float DeltaX = 1.0;
float DeltaY = 1.0;
float DeltaZ = 1.0;
bool PutDeltas = false;

void
SetDeltas( const char* arg )
{
  if ( 3 == sscanf( arg, "%15f,%15f,%15f", &DeltaX, &DeltaY, &DeltaZ ) )
    {
    PutDeltas = true;
    }
  else
    {
    cmtk::StdErr << "ERROR: could not parse X,Y,Z pixel size from argument '" << arg << "'\n";
    throw cmtk::ExitException( 1 );
    }
}

cmtk::ScalarDataType DataType = cmtk::TYPE_NONE;

float Offset = 0;
bool PutOffset = false;

bool ResetCoordinates = false;

const char* HdrFileName = "header.nii";
const char* ImportHdrFile = NULL;

const char* Description = NULL;

int 
doMain( const int argc, const char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Make NIFTI header file" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Make header file according to NIFTI file format based on user-supplied parameters for geometry, data type, orientation, etc." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "mk_nifti_hdr [options] [output.nii]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddCallback( Key( 'D', "dims" ), SetDims, "Set dimensions in voxels. Provided as 'DimsX,DimsY,DimsZ'" );
    cl.AddCallback( Key( 'V', "voxel" ), SetDeltas, "Set voxel size. Provided as 'PixX,PixY,PixZ' [in mm]" );

    cl.AddOption( Key( 'O', "offset" ), &Offset, "Binary data file offset", &PutOffset );

    cl.AddSwitch( Key( 'b', "byte" ), &DataType, cmtk::TYPE_BYTE, "8 bits, unsigned (this is the default)" );
    cl.AddSwitch( Key( 'c', "char" ), &DataType, cmtk::TYPE_CHAR, "8 bits, signed" );
    cl.AddSwitch( Key( 's', "short" ), &DataType, cmtk::TYPE_SHORT, "16 bits, signed" );
    cl.AddSwitch( Key( 'u', "ushort" ), &DataType, cmtk::TYPE_USHORT, "16 bits, unsigned" );
    cl.AddSwitch( Key( 'i', "int" ), &DataType, cmtk::TYPE_INT, "32 bits signed" );
    cl.AddSwitch( Key( 'f', "float" ), &DataType, cmtk::TYPE_FLOAT, "32 bits floating point" );
    cl.AddSwitch( Key( 'd', "double" ), &DataType, cmtk::TYPE_DOUBLE, "64 bits floating point\n" );
    
    cl.AddSwitch( Key( "reset-coordinates" ), &ResetCoordinates, true, "In an imported header, reset the image-to-physical coordinate mapping to an identity matrix." );
    cl.AddSwitch( Key( "attached" ), &DetachedHeader, false, "Create attached header for use in single-file images. Image data must be concatenated directly onto the created header file." );
    
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
  
  nifti_1_header header;
  memset( &header, 0, sizeof( header ) );

  if ( ImportHdrFile )
    {
    cmtk::CompressedStream hdrStream( ImportHdrFile );
    if ( !hdrStream.IsValid() ) 
      {
      cmtk::StdErr.printf( "ERROR: could not file %s\n", ImportHdrFile );
      throw cmtk::ExitException( 1 );
      }
    
    if ( sizeof(header) != hdrStream.Read( &header, 1, sizeof(header) ) ) 
      {
      cmtk::StdErr.printf( "ERROR: could not read %d bytes from header file %s\n", int( sizeof( header ) ), ImportHdrFile );
      throw cmtk::ExitException( 1 );
      }
    hdrStream.Close();
    }
  
  // set offset for binary file.
  if ( !ImportHdrFile || PutOffset )
    {
    cmtk::DebugOutput( 1 ) << "Setting image file offset\n";
    header.vox_offset = Offset;
    }

  header.sizeof_hdr = 348; // header size
  header.dim_info = 0;

  if ( ! ImportHdrFile || PutDims )
    {
    // ndims
    header.dim[0] = 4;
    
    // dimensions
    header.dim[1] = DimsX;
    header.dim[2] = DimsY;
    header.dim[3] = DimsZ;
    header.dim[4] = 1;
    header.dim[5] = 0;
    header.dim[6] = 0;
    header.dim[7] = 0;
    }

  if ( ! ImportHdrFile || PutDeltas )
    {
    header.pixdim[0] = 1;
    header.pixdim[1] = DeltaX;
    header.pixdim[2] = DeltaY;
    header.pixdim[3] = DeltaZ;
    header.pixdim[4] = 0.0;
    header.pixdim[5] = 0.0;
    }
  
  if ( ! ImportHdrFile || ResetCoordinates )
    {
    header.qform_code = 0;
    header.sform_code = 1;
    
    cmtk::AffineXform::MatrixType m4 = cmtk::AffineXform::MatrixType::Identity();
    for ( int i = 0; i < 4; ++i )
      {
      header.srow_x[i] = header.pixdim[1+i] * static_cast<float>( m4[i][0] );
      header.srow_y[i] = header.pixdim[1+i] * static_cast<float>( m4[i][1] );
      header.srow_z[i] = header.pixdim[1+i] * static_cast<float>( m4[i][2] );
      }
    }
  
  if ( !ImportHdrFile || DataType != cmtk::TYPE_NONE )
    {
    cmtk::DebugOutput( 1 ) << "Setting data type\n";

    switch ( DataType ) 
      {
      default:
      case cmtk::TYPE_BYTE:
	header.datatype = DT_UNSIGNED_CHAR;
	header.bitpix = 8 * sizeof(byte);
	break;
      case cmtk::TYPE_CHAR:
	header.datatype = DT_INT8;
	header.bitpix = 8 * sizeof(char);
	break;
      case cmtk::TYPE_SHORT:
	header.datatype = DT_INT16;
	header.bitpix = 8 * sizeof(short);
	break;
      case cmtk::TYPE_USHORT:
	header.datatype = DT_UINT16;
	header.bitpix = 8 * sizeof(unsigned short);
	break;
      case cmtk::TYPE_INT:
	header.datatype = DT_INT32;
	header.bitpix = 8 * sizeof(int);
	break;
      case cmtk::TYPE_FLOAT:
	header.datatype = DT_FLOAT;
	header.bitpix = 8 * sizeof(float);
	break;
      case cmtk::TYPE_DOUBLE:
	header.datatype = DT_DOUBLE;
	header.bitpix = 8 * sizeof(double);
	break;
      }
    }

  if ( Description )
    {
    cmtk::DebugOutput( 1 ) << "Setting image description\n";    
    memcpy( header.descrip, Description, std::min<size_t>( strlen( Description )+1, 80 ) );
    }

#ifdef _MSC_VER
  const char *const modestr = "wb";
#else
  const char *const modestr = "w";
#endif
    
  if ( DetachedHeader )
    {
    memcpy( &header.magic, "ni1\x00", 4 );
    header.vox_offset = 0;
    FILE *hdrFile = fopen( HdrFileName, modestr );
    if ( hdrFile ) 
      {
      fwrite( &header, 1, sizeof( header ), hdrFile );
      const int extension = 0;
      fwrite( &extension, 1, 4, hdrFile );
      fclose( hdrFile );
      }
    else
      {
      cmtk::StdErr << "ERROR: file '" << HdrFileName << "' could not be opened for writing!\n";
      }
    }
  else
    {
    memcpy( &header.magic, "n+1\x00", 4 );
    header.vox_offset = 352;
    FILE *hdrFile = fopen( HdrFileName, modestr );
    if ( hdrFile ) 
      {
      fwrite( &header, 1, sizeof( header ), hdrFile );
      const int extension = 0;
      fwrite( &extension, 1, 4, hdrFile );
      fclose( hdrFile );
      }
    }

  return 0;
}

#include "cmtkSafeMain"
