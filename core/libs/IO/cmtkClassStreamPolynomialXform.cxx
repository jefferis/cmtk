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
//  $Revision: 5046 $
//
//  $LastChangedDate: 2013-11-27 14:59:04 -0800 (Wed, 27 Nov 2013) $
//
//  $LastChangedBy: torsten_at_home $
//
*/

#include "cmtkClassStreamPolynomialXform.h"

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStreamOutput& 
operator << ( ClassStreamOutput& stream, const PolynomialXform& xform )
{
  stream.Begin( "polynomial_xform" );
  stream.WriteInt( "degree", xform.Degree() );
  stream.WriteCoordinateArray( "center", xform.Center().begin(), 3 );
  stream.WriteCoordinateArray( "coefficients", xform.m_Parameters, xform.m_NumberOfParameters );
  stream.End();

  return stream;
}
 
ClassStreamInput& 
operator >> ( ClassStreamInput& stream, PolynomialXform& xform )
{
  const char *fixedImage = NULL;
  const char *movingImage = NULL;

  if ( stream.Seek( "polynomial_xform", true /*forward*/ ) != TypedStream::CONDITION_OK )
    {
    stream.Rewind();
    if ( stream.Seek( "registration", true /*forward*/ ) != TypedStream::CONDITION_OK )
      {
      throw Exception( "Did not find 'registration' section in archive" );
      }
    
    fixedImage = stream.ReadString( "reference_study", NULL );
    movingImage = stream.ReadString( "floating_study", NULL );

    if ( stream.Seek( "polynomial_xform", false /*forward*/ ) != TypedStream::CONDITION_OK )
      {
      throw Exception( "Did not find 'polynomial_xform' section in archive" );
      }
    }

  int degree = stream.ReadInt( "degree", -1 );
  if ( degree == -1 )
    {
    throw Exception( "Did not find 'degree' value in polynomial xform archive" );
    }

  xform = PolynomialXform( degree );

  Types::Coordinate center[3];
  if ( stream.ReadCoordinateArray( "center", center, 3 ) != TypedStream::CONDITION_OK )
    {
    throw Exception( "Could not read 'center' array from polynomial xform archive" );
    }
  xform.SetCenter( PolynomialXform::SpaceVectorType::FromPointer( center ) );

  if ( stream.ReadCoordinateArray( "coefficients", xform.m_Parameters, xform.m_NumberOfParameters ) != TypedStream::CONDITION_OK )
    {
    throw Exception( "Could not read 'coeffients' array from polynomial xform archive" );
    }
  stream.End();

  xform.SetMetaInfo( META_SPACE, AnatomicalOrientation::ORIENTATION_STANDARD );
  if ( fixedImage )
    xform.SetMetaInfo( META_XFORM_FIXED_IMAGE_PATH, fixedImage );

  if ( movingImage )
    xform.SetMetaInfo( META_XFORM_MOVING_IMAGE_PATH, movingImage );

  return stream;
}

} // namespace cmtk
