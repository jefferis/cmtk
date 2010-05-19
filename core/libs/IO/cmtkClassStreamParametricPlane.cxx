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

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
ClassStream::operator>>( ParametricPlane*& infinitePlane )
{
  infinitePlane = NULL;

  if ( this->Seek( "plane" ) != TYPEDSTREAM_OK )
    return *this;
  
  infinitePlane = new ParametricPlane();

  Types::Coordinate planeOrigin[3];
  this->ReadCoordinateArray( "origin", planeOrigin, 3 );
  infinitePlane->SetOrigin( Vector3D( planeOrigin ) );

  infinitePlane->SetRho( this->ReadCoordinate( "rho" ) );
  infinitePlane->SetTheta( Units::Degrees( this->ReadCoordinate( "theta" ) ) );
  infinitePlane->SetPhi( Units::Degrees( this->ReadCoordinate( "phi" ) ) );

  return *this;
}

ClassStream&
ClassStream::operator<<( const ParametricPlane* infinitePlane )
{  
  this->Begin( "plane" );
  this->WriteCoordinateArray( "origin", infinitePlane->GetOrigin().begin(), 3 );
  this->WriteDouble( "rho", infinitePlane->GetRho() );
  this->WriteDouble( "theta", Units::Degrees( infinitePlane->GetTheta() ).Value() );
  this->WriteDouble( "phi", Units::Degrees( infinitePlane->GetPhi() ).Value() );

  this->WriteCoordinateArray( "normal", infinitePlane->GetNormal().begin(), 3 );
  return *this;
}

} // namespace cmtk
