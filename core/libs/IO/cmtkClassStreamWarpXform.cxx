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

#include <cmtkClassStream.h>
#include <cmtkClassStreamAffineXform.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
ClassStream::operator << ( const WarpXform *warpXform )
{
  const SplineWarpXform *splineWarpXform = dynamic_cast<const SplineWarpXform*>( warpXform );
  if ( splineWarpXform ) return (*this << splineWarpXform );

  return *this;
}

ClassStream& 
ClassStream::Put
( const WarpXform *warpXform,
  const AffineXform* initialXform )
{
  const SplineWarpXform* splineWarpXform = dynamic_cast<const SplineWarpXform*>( warpXform );
  if ( splineWarpXform )
    return this->PutWarp( splineWarpXform, initialXform );
  
  return *this;
}

ClassStream& 
ClassStream::Put
( const SplineWarpXform *splineWarpXform,
  const AffineXform* initialXform )
{
  this->Begin( "spline_warp" );
  return this->PutWarp( splineWarpXform, initialXform );
}

ClassStream& 
ClassStream::PutWarp
( const WarpXform *warpXform, 
  const AffineXform* initialXform )
{
  Types::Coordinate *nCoeff;
  // Undo initial transformation.
  if ( initialXform ) 
    {
    Types::Coordinate *oCoeff = warpXform->m_Parameters;
    nCoeff = Memory::AllocateArray<Types::Coordinate>( warpXform->m_NumberOfParameters );
    Types::Coordinate *p = nCoeff;
    for ( int z=0; z<warpXform->m_Dims[2]; ++z )
      for ( int y=0; y<warpXform->m_Dims[1]; ++y )
	for ( int x=0; x<warpXform->m_Dims[0]; ++x, oCoeff+=3, p+=3 ) 
	  {
	  Vector3D P( oCoeff );
	  initialXform->ApplyInPlace( P );
	  p[0] = P[0]; p[1] = P[1]; p[2] = P[2];
	  
	  // Undo offset tranformation; regain deltas.
	  p[0] -= (warpXform->m_Offset[0] + x * warpXform->Spacing[0]);
	  p[1] -= (warpXform->m_Offset[1] + y * warpXform->Spacing[1]);
	  p[2] -= (warpXform->m_Offset[2] + z * warpXform->Spacing[2]);
	  
	  // Is this correct, anyway? a) we may have to use the INVERSE
	  // affine transformation, as it is with respect to the model while
	  // the warp refers to the reference. b) we should do deformed-affine
	  // to get the relative offset.
	  
	  // Filter out unreasonably small values induced by rounding errors.
	  for ( int i=0; i<3; ++i ) 
	    {
	    if ( fabs( p[i] ) < 1e-5 )
	      p[i] = 0;
	    }
	  }
    } 
  else
    {
    nCoeff = warpXform->m_Parameters;
    }
  
  if ( warpXform->GetInitialAffineXform() )
    *this << (*warpXform->GetInitialAffineXform());
  
  this->WriteBool ( "absolute", (initialXform == NULL) );
  this->WriteIntArray( "dims", warpXform->m_Dims.begin(), 3 );

  this->WriteCoordinateArray( "domain", warpXform->Domain, 3 );
  this->WriteCoordinateArray( "origin", warpXform->m_Offset.XYZ, 3 );
  this->WriteCoordinateArray( "coefficients", nCoeff, warpXform->m_NumberOfParameters, 3 );

  const BitVector* activeFlags = warpXform->GetActiveFlags();
  if ( activeFlags ) 
    {
    this->WriteBoolArray( "active", activeFlags->GetBitVector(), warpXform->m_NumberOfParameters, 30 );
    }			 
  
  this->End();
  
  if ( initialXform ) delete[] nCoeff;
  
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
  if ( this->Seek( "spline_warp" ) == TYPEDSTREAM_OK ) 
    WarpType = 1;
  else
    if ( this->Seek( "linear_warp" ) == TYPEDSTREAM_OK )
      WarpType = 0;
    else 
      {
      this->Rewind();
      if ( this->Seek( "registration", true /*forward*/ ) != TYPEDSTREAM_OK )
	{
	return *this;
	}
      if ( this->Seek( "spline_warp" ) == TYPEDSTREAM_OK ) 
	WarpType = 1;
      else
	if ( this->Seek( "linear_warp" ) == TYPEDSTREAM_OK )
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
  if ( TYPEDSTREAM_OK != this->ReadIntArray( "dims", dims, 3 ) ) 
    {
    return *this;
    }
  
  int numControlPoints = dims[0] * dims[1] * dims[2];
  int numberOfParameters = 3 * numControlPoints;
  CoordinateVector::SmartPtr parameters( new CoordinateVector( numberOfParameters ) );
  Types::Coordinate *Coefficients = parameters->Elements;
  
  Types::Coordinate domain[3];
  Types::Coordinate origin[3];

  if ( this->ReadCoordinateArray( "domain", domain, 3 ) != TYPEDSTREAM_OK )
    this->ReadCoordinateArray( "extent", domain, 3 );
  
  int readOrigin = this->ReadCoordinateArray( "origin", origin, 3 );
  this->ReadCoordinateArray( "coefficients", Coefficients, numberOfParameters );
  if ( !absolute && (readOrigin == TYPEDSTREAM_OK) ) 
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
      warpXform = new SplineWarpXform( domain, dims, parameters, initialInverse );
      break;
    };
  
  byte *active = Memory::AllocateArray<byte>( (numberOfParameters / 8)+1 );
  if ( this->ReadBoolArray( "active", active, numberOfParameters ) == TYPEDSTREAM_OK ) 
    {
    BitVector::SmartPtr bitSet( new BitVector( numberOfParameters, active ) );
    warpXform->SetActiveFlags( bitSet );
    } 
  else 
    {
    delete[] active;
    }
  
  this->End();

  if ( warpXform )
    {
    warpXform->m_MetaInformation[META_SPACE] = AnatomicalOrientation::ORIENTATION_STANDARD;
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
