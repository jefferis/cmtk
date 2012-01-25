/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkClassStream.h"

#include <IO/cmtkClassStreamAffineXform.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
ClassStream::operator << ( const WarpXform *warpXform )
{
  return this->PutWarp( warpXform );
}

ClassStream& 
ClassStream::PutWarp
( const WarpXform *warpXform )
{
  const Types::Coordinate *nCoeff = warpXform->m_Parameters;
  
  if ( dynamic_cast<const SplineWarpXform*>( warpXform ) )
    this->Begin( "spline_warp" );
  
  if ( warpXform->GetInitialAffineXform() )
    *this << (*warpXform->GetInitialAffineXform());
  
  this->WriteBool ( "absolute", true );
  this->WriteIntArray( "dims", warpXform->m_Dims.begin(), 3 );
  
  this->WriteCoordinateArray( "domain", warpXform->m_Domain.begin(), 3 );
  this->WriteCoordinateArray( "origin", warpXform->m_Offset.begin(), 3 );
  this->WriteCoordinateArray( "coefficients", nCoeff, warpXform->m_NumberOfParameters, 3 );
  
  const BitVector* activeFlags = warpXform->GetActiveFlags();
  if ( activeFlags ) 
    {
    this->WriteBoolArray( "active", activeFlags->GetBitVector(), warpXform->m_NumberOfParameters, 30 );
    }			 
  
  this->End();
  
  return *this;
}

ClassStream&
ClassStream::Get
( WarpXform::SmartPtr& warpXform, const AffineXform* initialXform )
{
  WarpXform* warp;
  this->Get( warp, initialXform );
  warpXform = WarpXform::SmartPtr( warp );
  return *this;
}

ClassStream&
ClassStream::Get
( WarpXform*& warpXform, const AffineXform* initialXform )
{
  warpXform = NULL;

  int WarpType = -1;
  if ( this->Seek( "spline_warp" ) == TypedStream::CONDITION_OK ) 
    WarpType = 1;
  else
    if ( this->Seek( "linear_warp" ) == TypedStream::CONDITION_OK )
      WarpType = 0;
    else 
      {
      this->Rewind();
      if ( this->Seek( "registration", true /*forward*/ ) != TypedStream::CONDITION_OK )
	{
	return *this;
	}
      if ( this->Seek( "spline_warp" ) == TypedStream::CONDITION_OK ) 
	WarpType = 1;
      else
	if ( this->Seek( "linear_warp" ) == TypedStream::CONDITION_OK )
	  WarpType = 0;
	else
	  return *this;
      }
  
  AffineXform::SmartPtr initialInverse( NULL );
  if ( initialXform == NULL ) 
    {
    AffineXform::SmartPtr newInitialXform;
    *this >> newInitialXform;
    initialInverse = AffineXform::SmartPtr( newInitialXform );
    } 
  else 
    {
    initialInverse = AffineXform::SmartPtr( initialXform->MakeInverse() );
    }
  
  int absolute = this->ReadBool( "absolute", 0 );
  
  int dims[3];
  if ( TypedStream::CONDITION_OK != this->ReadIntArray( "dims", dims, 3 ) ) 
    {
    return *this;
    }
  
  int numControlPoints = dims[0] * dims[1] * dims[2];
  int numberOfParameters = 3 * numControlPoints;
  CoordinateVector::SmartPtr parameters( new CoordinateVector( numberOfParameters ) );
  Types::Coordinate *Coefficients = parameters->Elements;
  
  Vector3D domain;
  Vector3D origin;

  if ( this->ReadCoordinateArray( "domain", domain.begin(), 3 ) != TypedStream::CONDITION_OK )
    this->ReadCoordinateArray( "extent", domain.begin(), 3 );
  
  int readOrigin = this->ReadCoordinateArray( "origin", origin.begin(), 3 );
  this->ReadCoordinateArray( "coefficients", Coefficients, numberOfParameters );
  if ( !absolute && (readOrigin == TypedStream::CONDITION_OK) ) 
    {
    Types::Coordinate *p = Coefficients;
    for ( int z=0; z<dims[2]; ++z )
      for ( int y=0; y<dims[1]; ++y )
	for ( int x=0; x<dims[0]; ++x, p+=3 ) 
	  {
	  if ( WarpType == 0 ) 
	    {
	    p[0] += (origin[0] + x * domain[0]/(dims[0]-1));
	    p[1] += (origin[1] + y * domain[1]/(dims[1]-1));
	    p[2] += (origin[2] + z * domain[2]/(dims[2]-1));
	    } 
	  else
	    {
	    p[0] += (origin[0] + x * domain[0]/(dims[0]-3));
	    p[1] += (origin[1] + y * domain[1]/(dims[1]-3));
	    p[2] += (origin[2] + z * domain[2]/(dims[2]-3));
	    }
	  }
    }
  
  switch ( WarpType ) 
    {
    case 0: 
      warpXform = NULL; // linear warp no longer supported
      break;
    case 1: 
      warpXform = new SplineWarpXform( domain, SplineWarpXform::IndexType( dims ), parameters, initialInverse );
      break;
    };
  
  byte *active = Memory::ArrayC::Allocate<byte>( (numberOfParameters / 8)+1 );
  if ( this->ReadBoolArray( "active", active, numberOfParameters ) == TypedStream::CONDITION_OK ) 
    {
    BitVector::SmartPtr bitSet( new BitVector( numberOfParameters, active ) );
    warpXform->SetActiveFlags( bitSet );
    } 
  else 
    {
    Memory::ArrayC::Delete( active );
    }
  
  this->End();

  if ( warpXform )
    {
    warpXform->SetMetaInfo( META_SPACE, AnatomicalOrientation::ORIENTATION_STANDARD );
    }
  
  return *this;
}

ClassStream& 
ClassStream::operator >> ( WarpXform::SmartPtr& warpXform )
{
  this->Get( warpXform );
  return *this; 
}

ClassStream& 
ClassStream::operator >> ( WarpXform*& warpXform )
{
  this->Get( warpXform );
  return *this; 
}

} // namespace cmtk
