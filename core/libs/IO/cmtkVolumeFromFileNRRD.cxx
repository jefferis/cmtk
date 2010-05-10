/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifdef CMTK_BUILD_NRRD

#include <cmtkVolumeIO.h>
#include <cmtkVolumeFromFile.h>
#include <cmtkTypes.h>

#include <cmtkConsole.h>

#include <cmtkUniformVolume.h>
#include <cmtkAnatomicalOrientation.h>

#ifdef CMTK_BUILD_NRRD_TEEM
#  include <teem/nrrd.h>
#else
#  include <NrrdIO.h>
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

const UniformVolume::SmartPtr
VolumeFromFile::ReadNRRD( const char* pathHdr )
{
  UniformVolume::SmartPtr volume( NULL );
  try 
    {
    Nrrd *nrrd = nrrdNew();
    if ( nrrdLoad( nrrd, pathHdr, NULL ) )
      throw biffGetDone(NRRD);

    if ( nrrd->dim > 3 )
      {
      StdErr << "ERROR: for now, nrrd input can only handle data with dimension 3 or less.\n";
      return UniformVolume::SmartPtr( NULL );
      }

    const int dims[3] = 
      { 
	(nrrd->dim > 0) ? nrrd->axis[0].size : 1,
	(nrrd->dim > 1) ? nrrd->axis[1].size : 1,
	(nrrd->dim > 2) ? nrrd->axis[2].size : 1 
      };

    // for each axis, if spacing is NaN, use direction vector to compute spacing.
    double spacing[3] = { 1, 1, 1 };
    for ( size_t ax = 0; ax < nrrd->dim; ++ax )
      {
      switch ( nrrdSpacingCalculate( nrrd, ax, spacing+ax, nrrd->axis[ax].spaceDirection ) )
	{
	case nrrdSpacingStatusScalarNoSpace:
	  break;
	case nrrdSpacingStatusDirection:
	  break;
	case nrrdSpacingStatusScalarWithSpace:
	  StdErr << "WARNING: nrrdSpacingCalculate returned nrrdSpacingStatusScalarWithSpace\n";
	  spacing[ax] = nrrd->axis[ax].spacing;
	  break;
	case nrrdSpacingStatusNone:
	default:
	  StdErr << "WARNING: no pixel spacings in Nrrd for axis " << ax << "; setting to 1.0\n";
	  spacing[ax] = 1.0;
	  break;
	}
      }
    const Types::Coordinate size[3] = { (dims[0]-1) * spacing[0], (dims[1]-1) * spacing[1], (dims[2]-1) * spacing[2] };
    volume = UniformVolume::SmartPtr( new UniformVolume( DataGrid::IndexType( dims ),UniformVolume::CoordinateVectorType( size ) ) );

    ScalarDataType type = TYPE_NONE;
    switch ( nrrd->type )
      {
      case nrrdTypeUChar:  type = TYPE_BYTE; break;
      case nrrdTypeChar:   type = TYPE_CHAR; break;
      case nrrdTypeUShort: type = TYPE_USHORT; break;
      case nrrdTypeShort:  type = TYPE_SHORT; break;
      case nrrdTypeInt:    type = TYPE_INT; break;
      case nrrdTypeFloat:  type = TYPE_FLOAT; break;
      case nrrdTypeDouble: type = TYPE_DOUBLE; break;
      default: break;
      }

    if ( type != TYPE_NONE )
      {
      TypedArray::SmartPtr data( TypedArray::Create( type, nrrd->data, volume->GetNumberOfPixels() ) );
      volume->SetData( data );
      }
    else
      {
      StdErr << "ERROR: unsupported data type in nrrd file.\n";
      return volume;
      }

    const char* orientationSpaceAnatomical = NULL;
    switch ( nrrd->space )
      {
      case nrrdSpaceRightAnteriorSuperior:
      case nrrdSpaceRightAnteriorSuperiorTime:
	orientationSpaceAnatomical = "RAS";
	volume->m_MetaInformation[META_SPACE] = orientationSpaceAnatomical;
	break;
      case nrrdSpaceLeftAnteriorSuperior:
      case nrrdSpaceLeftAnteriorSuperiorTime:
	orientationSpaceAnatomical = "LAS";
	volume->m_MetaInformation[META_SPACE] = orientationSpaceAnatomical;
	break;
      case nrrdSpaceLeftPosteriorSuperior:
      case nrrdSpaceLeftPosteriorSuperiorTime:
	orientationSpaceAnatomical = "LPS";
	volume->m_MetaInformation[META_SPACE] = orientationSpaceAnatomical;
	break;
      case nrrdSpace3DRightHanded:
	volume->m_MetaInformation[META_SPACE] = "3DRH";
	break;
      case nrrdSpace3DLeftHanded:
	volume->m_MetaInformation[META_SPACE] = "3DLH";
	break;
      case nrrdSpace3DRightHandedTime:
	volume->m_MetaInformation[META_SPACE] = "3DRHT";
	break;
      case nrrdSpace3DLeftHandedTime:
	volume->m_MetaInformation[META_SPACE] = "3DLHT";
	break;
      default:
	break;
      }

    volume->m_MetaInformation[META_SPACE_ORIGINAL] = volume->m_MetaInformation[META_SPACE];
    
    const Types::Coordinate directions[3][3] = 
      {
	{ nrrd->axis[0].spaceDirection[0] * spacing[0],
	  nrrd->axis[0].spaceDirection[1] * spacing[0],
	  nrrd->axis[0].spaceDirection[2] * spacing[0] },
	{ nrrd->axis[1].spaceDirection[0] * spacing[1], 
	  nrrd->axis[1].spaceDirection[1] * spacing[1],
	  nrrd->axis[1].spaceDirection[2] * spacing[1] },
	{ nrrd->axis[2].spaceDirection[0] * spacing[2], 
	  nrrd->axis[2].spaceDirection[1] * spacing[2], 
	  nrrd->axis[2].spaceDirection[2] * spacing[2] }
      };
    
    const Matrix3x3<Types::Coordinate> m3( directions );
    Matrix4x4<Types::Coordinate> m4( m3 );
    for ( int i = 0; i < 3; ++i )
      m4[3][i] = nrrd->spaceOrigin[i];
    volume->m_IndexToPhysicalMatrix = m4;
    
    if ( orientationSpaceAnatomical )
      {        
      char orientationImage[4];
      AnatomicalOrientation::GetOrientationFromDirections( orientationImage, m4, orientationSpaceAnatomical );
      volume->m_MetaInformation[META_IMAGE_ORIENTATION] = volume->m_MetaInformation[META_IMAGE_ORIENTATION_ORIGINAL] = orientationImage;
      volume->ChangeCoordinateSpace( AnatomicalOrientation::ORIENTATION_STANDARD );
      }
    else
      {
      }

    if ( nrrd->spaceUnits[0] )
      volume->m_MetaInformation[META_SPACE_UNITS_STRING] = nrrd->spaceUnits[0];

    nrrdNix( nrrd );
    }
  catch ( char* err )
    {
    StdErr << "ERROR: nrrd library returned error '" << err << "'\n";
    free( err );
    }

  return volume;
}


void
VolumeFromFile::WriteNRRD
( const char* pathHdr, const UniformVolume& volume, const bool )
{
  UniformVolume::SmartPtr writeVolume( volume.Clone() );
  if ( writeVolume->MetaKeyExists( META_SPACE_ORIGINAL ) )
    writeVolume->ChangeCoordinateSpace( writeVolume->m_MetaInformation[META_SPACE_ORIGINAL] );
  
  void* val = const_cast<void*>( writeVolume->GetData()->GetDataPtr() );
  int type = nrrdTypeUnknown;
  switch ( writeVolume->GetData()->GetType() )
    {
    case TYPE_BYTE:   type = nrrdTypeUChar; break;
    case TYPE_CHAR:   type = nrrdTypeChar; break;
    case TYPE_USHORT: type = nrrdTypeUShort; break;
    case TYPE_SHORT:  type = nrrdTypeShort; break;
    case TYPE_INT:    type = nrrdTypeInt; break;
    case TYPE_FLOAT:  type = nrrdTypeFloat; break;
    case TYPE_DOUBLE: type = nrrdTypeDouble; break;
    default: break;
    }
  
  Nrrd *nval = nrrdNew();
  NrrdIoState *nios = nrrdIoStateNew();

  if ( VolumeIO::GetWriteCompressed() && nrrdEncodingGzip->available() )
    {
    nrrdIoStateEncodingSet( nios, nrrdEncodingGzip );
    nrrdIoStateSet( nios, nrrdIoStateZlibLevel, 9 );
    }
  else
    {
    StdErr << "WARNING: Nrrd library does not support Gzip compression encoding.\n"
	      << " Please add -DTEEM_ZLIB to compiler options when building Nrrd library.\n";
    }
  
  try
    {
    if ( nrrdWrap_va( nval, val, type, (size_t)3, (size_t)writeVolume->m_Dims[0], (size_t)writeVolume->m_Dims[1], (size_t)writeVolume->m_Dims[2] ) )
      {
      throw( biffGetDone(NRRD) );
      }

    nrrdSpaceDimensionSet( nval, 3 );

    if ( writeVolume->MetaKeyExists(META_SPACE_UNITS_STRING) )
      {
      nval->spaceUnits[0] = strdup( writeVolume->m_MetaInformation[META_SPACE_UNITS_STRING].c_str() );
      nval->spaceUnits[1] = strdup( writeVolume->m_MetaInformation[META_SPACE_UNITS_STRING].c_str() );
      nval->spaceUnits[2] = strdup( writeVolume->m_MetaInformation[META_SPACE_UNITS_STRING].c_str() );
      }
      
    int kind[NRRD_DIM_MAX] = { nrrdKindDomain, nrrdKindDomain, nrrdKindDomain };
    nrrdAxisInfoSet_nva( nval, nrrdAxisInfoKind, kind );

    const std::string space = writeVolume->m_MetaInformation[META_SPACE];
    
    // if the volume has a direction table, write it to the Nrrd
    if ( space == "RAS" )
      {
      nval->space = nrrdSpaceRightAnteriorSuperior;
      }
    else if ( space == "LAS" )
      {
      nval->space = nrrdSpaceLeftAnteriorSuperior;
      }
    else if ( space == "LPS" )
      {
      nval->space = nrrdSpaceLeftPosteriorSuperior;
      }
    else if ( space == "3DRH" )
      {
      nval->space = nrrdSpace3DRightHanded;
      }
    else if ( space == "3DLH" )
      {
      nval->space = nrrdSpace3DLeftHanded;
      }
    else if ( space == "3DRHT" )
      {
      nval->space = nrrdSpace3DRightHandedTime;
      }
    else if ( space == "3DLHT" )
      {
      nval->space = nrrdSpace3DLeftHandedTime;
      }
    else
      {
      if ( space.length() == 3 )
	{
	writeVolume->ChangeCoordinateSpace( "RAS" );
	nval->space = nrrdSpaceRightAnteriorSuperior;
	}
      }
    
    const AffineXform::MatrixType& matrix = writeVolume->m_IndexToPhysicalMatrix;
    double spaceDir[NRRD_DIM_MAX][NRRD_SPACE_DIM_MAX];
    for ( int i = 0; i < 3; ++i ) 
      {
      for ( int j = 0; j < 3; ++j )
	spaceDir[i][j] = matrix[i][j];
      }
    nrrdAxisInfoSet_nva( nval, nrrdAxisInfoSpaceDirection, spaceDir );
    
    double origin[NRRD_DIM_MAX] = { matrix[3][0], matrix[3][1], matrix[3][2] };
    if ( nrrdSpaceOriginSet( nval, origin ) )
      {
      throw( biffGetDone(NRRD) );
      }	
  
    nrrdAxisInfoSet_va( nval, nrrdAxisInfoLabel, "x", "y", "z" );

    if ( nrrdSave( pathHdr, nval, nios ) )
      {
      throw( biffGetDone(NRRD) );
      }
    }
  catch ( char* err )
    {
    StdErr << "ERROR: NrrdIO library returned error '" << err << "'\n";
    free( err );
    }
  
  nrrdIoStateNix( nios );
  nrrdNix(nval);    
}

//@}

} // namespace cmtk

#endif // #ifdef CMTK_BUILD_NRRD
