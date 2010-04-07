/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkXformIO.h>

#include <cmtkDeformationField.h>

#ifdef CMTK_BUILD_NRRD_TEEM
#  include <teem/nrrd.h>
#else
#  ifdef CMTK_BUILD_NRRD
#    include <NrrdIO.h>
#  endif
#endif

#ifdef CMTK_BUILD_NRRD

namespace
cmtk
{

/** \addtogroup IO */
//@{

Xform::SmartPtr
XformIO::ReadNrrd( const char* path, const bool )
{
  DeformationField::SmartPtr dfield( NULL );
  try 
    {
    Nrrd *nrrd = nrrdNew();
    if ( nrrdLoad( nrrd, path, NULL ) )
      throw biffGetDone(NRRD);

    if ( nrrd->dim != 4 )
      {
      StdErr << "ERROR: deformation field must be stored as 4-dimensional Nrrd.\n";
      return dfield;
      }

    if ( nrrd->axis[0].kind != nrrdKindVector )
      {
      StdErr << "ERROR: deformation field vectors in Nrrd must be stored together.\n";
      return dfield;
      }
    
    if ( nrrd->axis[0].size != 3 )
      {
      StdErr << "ERROR: deformation field vectors in Nrrd must be three dimensional.\n";
      return dfield;
      }

    NrrdAxisInfo* nrrdSpaceAxes = nrrd->axis+1;
    const int dims[3] = { nrrdSpaceAxes[0].size, nrrdSpaceAxes[1].size, nrrdSpaceAxes[2].size };

    // for each axis, if spacing is NaN, use direction vector to compute spacing.
    double spacing[3] = { 1, 1, 1 };
    for ( size_t ax = 0; ax < 3; ++ax )
      {
      switch ( nrrdSpacingCalculate( nrrd, ax+1, spacing+ax, nrrd->axis[ax+1].spaceDirection ) )
	{
	case nrrdSpacingStatusScalarNoSpace:
	  break;
	case nrrdSpacingStatusDirection:
	  break;
	case nrrdSpacingStatusScalarWithSpace:
	  StdErr << "WARNING: nrrdSpacingCalculate returned nrrdSpacingStatusScalarWithSpace\n";
	  spacing[ax] = nrrdSpaceAxes[ax].spacing;
	  break;
	case nrrdSpacingStatusNone:
	default:
	  StdErr << "WARNING: no pixel spacings in Nrrd for axis " << ax << "; setting to 1.0\n";
	  spacing[ax] = 1.0;
	  break;
	}
      }

    const Types::Coordinate size[3] = { (dims[0]-1) * spacing[0], (dims[1]-1) * spacing[1], (dims[2]-1) * spacing[2] };
    const Types::Coordinate origin[3] = { nrrd->spaceOrigin[0], nrrd->spaceOrigin[1], nrrd->spaceOrigin[2] };
    dfield = DeformationField::SmartPtr( new DeformationField( size, dims, origin ) );
    
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
      TypedArray::SmartPtr data( TypedArray::Create( type, nrrd->data, 3 * dims[0] * dims[1] * dims[2] ) );
      data->ConvertSubArray( dfield->m_Parameters, TYPE_COORDINATE, 0, data->GetDataSize() );
      }
    else
      {
      StdErr << "ERROR: unsupported data type in nrrd file.\n";
      return dfield;
      }

    const char* orientationSpace = NULL;
    switch ( nrrd->space )
      {
      case nrrdSpaceRightAnteriorSuperior:
      case nrrdSpaceRightAnteriorSuperiorTime:
	orientationSpace = "RAS";
	break;
      case nrrdSpaceLeftAnteriorSuperior:
      case nrrdSpaceLeftAnteriorSuperiorTime:
	orientationSpace = "LAS";
	break;
      case nrrdSpaceLeftPosteriorSuperior:
      case nrrdSpaceLeftPosteriorSuperiorTime:
	orientationSpace = "LPS";
	break;
      default:
	break;
      }

    if ( orientationSpace )
      {
      dfield->m_MetaInformation[CMTK_META_SPACE] = dfield->m_MetaInformation[CMTK_META_SPACE_ORIGINAL] = orientationSpace;
      
      const Types::Coordinate directions[3][3] = 
	{
	  { nrrdSpaceAxes[0].spaceDirection[0] * spacing[0], 
	    nrrdSpaceAxes[0].spaceDirection[1] * spacing[0], 
	    nrrdSpaceAxes[0].spaceDirection[2] * spacing[0] },
	  { nrrdSpaceAxes[1].spaceDirection[0] * spacing[1], 
	    nrrdSpaceAxes[1].spaceDirection[1] * spacing[1], 
	    nrrdSpaceAxes[1].spaceDirection[2] * spacing[1] },
	  { nrrdSpaceAxes[2].spaceDirection[0] * spacing[2], 
	    nrrdSpaceAxes[2].spaceDirection[1] * spacing[2], 
	    nrrdSpaceAxes[2].spaceDirection[2] * spacing[2] }
	};

      const Matrix3x3<Types::Coordinate> m3( directions );
      Matrix4x4<Types::Coordinate> m4( m3 );
      for ( int i = 0; i < 3; ++i )
	m4[3][i] = nrrd->spaceOrigin[i];

      AffineXform::SmartPtr xform( new AffineXform( m4 ) ) ;
      dfield->SetInitialAffineXform( xform );
      
      char orientationImage[4];
      AnatomicalOrientation::GetOrientationFromDirections( orientationImage, m4, orientationSpace );
      dfield->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION] = orientationImage;
      dfield->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION_ORIGINAL] = orientationImage;
      }
    
    nrrdNix( nrrd );
    }
  catch ( char* err )
    {
    StdErr << "ERROR: nrrd library returned error '" << err << "'\n";
    free( err );
    }

  return dfield;
}

void 
XformIO::WriteNrrd
( const Xform* xform, const char *path, const bool )
{
  const DeformationField* dfield = dynamic_cast<const DeformationField*>( xform );
  if ( ! dfield )
    {
    StdErr << "ERROR: XformIO::WriteNrrd can only write DeformationField objects so far.\n"
	      << "       No data was written.\n";
    return;
    }

  void* val = static_cast<void*>( dfield->m_Parameters );
  const int type = (sizeof(Types::Coordinate) == sizeof(float)) ? nrrdTypeFloat : nrrdTypeDouble;
  
  Nrrd *nval = nrrdNew();
  NrrdIoState *nios = nrrdIoStateNew();
  
  if ( nrrdEncodingGzip->available() )
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
    if ( nrrdWrap_va( nval, val, type, 4, (size_t)3, (size_t)dfield->m_Dims[0], (size_t)dfield->m_Dims[1], (size_t)dfield->m_Dims[2] ) )
      {
      throw( biffGetDone(NRRD) );
      }
    
    nrrdSpaceDimensionSet( nval, 3 );
    
    if ( dfield->MetaKeyExists(CMTK_META_SPACE_UNITS_STRING) )
      {
      nval->spaceUnits[0] = strdup( dfield->m_MetaInformation[CMTK_META_SPACE_UNITS_STRING].c_str() );
      }
      
    int kind[NRRD_DIM_MAX] = { nrrdKindVector, nrrdKindDomain, nrrdKindDomain, nrrdKindDomain };
    nrrdAxisInfoSet_nva( nval, nrrdAxisInfoKind, kind );
    nrrdAxisInfoSet_va( nval, nrrdAxisInfoLabel, "Vx;Vy;Vz", "x", "y", "z" );
    
    double origin[NRRD_DIM_MAX] = { dfield->m_Offset.XYZ[0], dfield->m_Offset.XYZ[1], dfield->m_Offset.XYZ[2] };
    if ( nrrdSpaceOriginSet( nval, origin ) )
      {
      throw( biffGetDone(NRRD) );
      }
    
    nval->space = nrrdSpaceRightAnteriorSuperior;
    double spaceDir[NRRD_DIM_MAX][NRRD_SPACE_DIM_MAX];
    for ( int i = 0; i < 4; ++i )
      {
      for ( int j = 0; j < 3; ++j )
	{
	if ( i )
	  {
	  if ( i-1 == j )
	    spaceDir[i][j] = dfield->Spacing[ i-1 ];
	  else
	    spaceDir[i][j] = 0.0;
	  }
	else
	  {
	  spaceDir[i][j] = AIR_NAN;
	  }
	}
      }
    nrrdAxisInfoSet_nva( nval, nrrdAxisInfoSpaceDirection, spaceDir );
    
    if ( nrrdSave( path, nval, nios ) )
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
