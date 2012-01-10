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

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
ClassStream::operator>>( ParametricPlane*& parametricPlane )
{
  parametricPlane = NULL;

  if ( this->Seek( "plane" ) != TypedStream::CONDITION_OK )
    return *this;
  
  parametricPlane = new ParametricPlane();

  Types::Coordinate planeOrigin[3];
  this->ReadCoordinateArray( "origin", planeOrigin, 3 );
  parametricPlane->SetOrigin( FixedVector<3,Types::Coordinate>( planeOrigin ) );

  parametricPlane->SetRho( this->ReadCoordinate( "rho" ) );
  parametricPlane->SetTheta( Units::Degrees( this->ReadCoordinate( "theta" ) ) );
  parametricPlane->SetPhi( Units::Degrees( this->ReadCoordinate( "phi" ) ) );

  return *this;
}

ClassStream&
ClassStream::operator<<( const ParametricPlane* parametricPlane )
{  
  this->Begin( "plane" );
  this->WriteCoordinateArray( "origin", parametricPlane->GetOrigin().begin(), 3 );
  this->WriteDouble( "rho", parametricPlane->GetRho() );
  this->WriteDouble( "theta", Units::Degrees( parametricPlane->GetTheta() ).Value() );
  this->WriteDouble( "phi", Units::Degrees( parametricPlane->GetPhi() ).Value() );

  this->WriteCoordinateArray( "normal", parametricPlane->GetNormal().begin(), 3 );
  return *this;
}

} // namespace cmtk
