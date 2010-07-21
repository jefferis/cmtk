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

#include "cmtkClassStreamAffineXform.h"

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
operator << ( ClassStream& stream, const AffineXform& affineXform )
{
  stream.Begin( "affine_xform" );
  stream.WriteCoordinateArray( "xlate", affineXform.RetXlate(), 3 );
  stream.WriteCoordinateArray( "rotate", affineXform.RetAngles(), 3 );
  if ( affineXform.GetUseLogScaleFactors() )
    stream.WriteCoordinateArray( "log_scale", affineXform.RetScales(), 3 );
  else
    stream.WriteCoordinateArray( "scale", affineXform.RetScales(), 3 );
  stream.WriteCoordinateArray( "shear", affineXform.RetShears(), 3 );
  stream.WriteCoordinateArray( "center", affineXform.RetCenter(), 3 );
  stream.End();

  return stream;
}
 
ClassStream& 
operator >> ( ClassStream& stream, AffineXform::SmartPtr& affineXform )
{
  try
    {
    affineXform = AffineXform::SmartPtr( new AffineXform );
    stream >> (*affineXform);
    }
  catch (...)
    {
    affineXform = AffineXform::SmartPtr::Null;
    }
  return stream;
}

ClassStream& 
operator >> ( ClassStream& stream, AffineXform& affineXform )
{
  CoordinateVector pVector( 15 );
  Types::Coordinate* parameters = pVector.Elements;

  if ( stream.Seek( "affine_xform", true /*forward*/ ) != TYPEDSTREAM_OK )
    {
    stream.Rewind();
    if ( stream.Seek( "registration", true /*forward*/ ) != TYPEDSTREAM_OK )
      {
      throw Exception( "Did not find 'registration' section in affine xform archive" );
      }
    if ( stream.Seek( "affine_xform", true /*forward*/ ) != TYPEDSTREAM_OK )
      {
      throw Exception( "Did not find 'affine_xform' section in affine xform archive" );
      }
    }

  if ( stream.ReadCoordinateArray( "xlate", parameters, 3 ) != TYPEDSTREAM_OK )
    {
    parameters[0] = parameters[1] = parameters[2] = 0;
    }
  if ( stream.ReadCoordinateArray( "rotate", parameters+3, 3 ) != TYPEDSTREAM_OK )
    {
    parameters[3] = parameters[4] = parameters[5] = 0;
    }
  bool logScaleFactors = false;
  if ( stream.ReadCoordinateArray( "scale", parameters+6, 3 ) != TYPEDSTREAM_OK )
    {
    if ( stream.ReadCoordinateArray( "log_scale", parameters+6, 3 ) == TYPEDSTREAM_OK )
      {
      logScaleFactors = true;
      }
    else
      {
      parameters[6] = parameters[7] = parameters[8] = 1;
      }
    }
  if ( stream.ReadCoordinateArray( "shear", parameters+9, 3 ) != TYPEDSTREAM_OK )
    {
    parameters[9] = parameters[10] = parameters[11] = 0;
    }
  if ( stream.ReadCoordinateArray( "center", parameters+12, 3 ) != TYPEDSTREAM_OK )
    {
    parameters[12] = parameters[13] = parameters[14] = 0;
    }
  stream.End();

  affineXform.SetUseLogScaleFactors( logScaleFactors );
  affineXform.SetParamVector( pVector );
  affineXform.m_MetaInformation[META_SPACE] = AnatomicalOrientation::ORIENTATION_STANDARD;

  return stream;
}

} // namespace cmtk
