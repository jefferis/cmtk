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

#include "cmtkPlaneSource.h"

#include "Base/cmtkMathUtil.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

PlaneSource::PlaneSource()
{
  Direction = 0;
  Position = 0;
  Resolution = 1;
  Input[0] = Input[1] = NULL;
  this->RegisterInput( &Input[0] );
  this->RegisterInput( &Input[1] );
  ReferenceVolumeIndex = 0;
}

int
PlaneSource::HasValidInputs() const
{
  if ( (Input[0] == NULL) || (Input[1] == NULL) ) return 0;
  Input[0]->Update();
  Input[1]->Update();
  if ( ! Input[0]->GetVolume() ) return 0;
  if ( ! Input[1]->GetVolume() ) return 0;
  if ( ! (Input[0]->GetVolume()->GetData()) ) return 0;
  if ( ! (Input[1]->GetVolume()->GetData()) ) return 0;
  return 1;
} 

void
PlaneSource::SetInput
( const int index, VolumeWrapper *const input )
{
  if ( (index < 0) || (index > 1) ) return;
 
  this->ReplaceObject( Input[index], input );
}

Types::Coordinate
PlaneSource::GetMinPosition()
{
  return 0;
}

Types::Coordinate
PlaneSource::GetMaxPosition()
{
  if ( !Input[ReferenceVolumeIndex] ) return 0;

  this->Update();

  const Volume *volume = Input[ReferenceVolumeIndex]->GetVolume();
  if ( volume == NULL ) return 0;

  switch ( Direction ) 
    {
    case SCANDIRECTION_CAUDAL_CRANIAL:
    case SCANDIRECTION_CRANIAL_CAUDAL:
      return volume->Size[2];
    case SCANDIRECTION_RIGHT_LEFT:
    case SCANDIRECTION_LEFT_RIGHT:
      return volume->Size[0];
    case SCANDIRECTION_VENTRAL_DORSAL:
    case SCANDIRECTION_DORSAL_VENTRAL:
      return volume->Size[1];
    }
  // default:
  return 0;
}

Types::Coordinate
PlaneSource::GetMaxResolution()
{
  this->Update();
  
  double maxResolution = 1000;
  for ( int i=0; i<2; ++i ) 
    {
    if ( Input[i] ) 
      {
      const UniformVolume *volume = Input[i]->GetVolume();
      if ( volume ) 
	maxResolution = std::min<Types::Coordinate>( maxResolution, volume->GetMinDelta() );
      }
    }
  
  // enforce a somehow reasonable minimum pixel size.
  return std::max<Types::Coordinate>( 0.05, maxResolution );
}

void
PlaneSource::Execute()
{
  if ( !Input[ReferenceVolumeIndex] ) return;
  const Volume *volume = Input[ReferenceVolumeIndex]->GetVolume();
  if ( volume == NULL ) return;

  Plane *output = this->GetOutput();

  if ( Position < 0 ) Position = 0;
  output->SetSpacing( Resolution, Resolution );

  switch ( Direction ) 
    {
    case SCANDIRECTION_CAUDAL_CRANIAL:
    case SCANDIRECTION_CRANIAL_CAUDAL:
      if ( Position > volume->Size[2] ) Position = volume->Size[2];
      
      output->SetDims( 1+(int)(volume->Size[0]/Resolution), 1+(int)(volume->Size[1]/Resolution) );
      
      output->SetDirectionY( 0, 1, 0 );
      
      if ( Direction == SCANDIRECTION_CAUDAL_CRANIAL ) 
	{
	output->SetOrigin( 0, 0, Position );
	output->SetDirectionX( 1, 0, 0 );
	} 
      else
	{
	output->SetOrigin( volume->Size[0], 0, Position );
	output->SetDirectionX( -1, 0, 0 );
	}
      break;
    case SCANDIRECTION_RIGHT_LEFT:
    case SCANDIRECTION_LEFT_RIGHT:
      if ( Position > volume->Size[0] ) Position = volume->Size[0];
      
      output->SetDims( 1+(int)(volume->Size[1]/Resolution), 1+(int)(volume->Size[2]/Resolution) );
      
      output->SetDirectionY( 0, 0, -1 );
      
      if ( Direction == SCANDIRECTION_RIGHT_LEFT ) 
	{
	output->SetOrigin( Position, 0, volume->Size[2] );
	output->SetDirectionX( 0, 1, 0 );
	} 
      else
	{
	output->SetOrigin( Position, volume->Size[1], volume->Size[2] );
	output->SetDirectionX( 0, -1, 0 );
	}
      break;
    case SCANDIRECTION_VENTRAL_DORSAL:
    case SCANDIRECTION_DORSAL_VENTRAL:
      if ( Position > volume->Size[1] ) Position = volume->Size[1];
      
      output->SetDims( 1+(int)(volume->Size[0]/Resolution), 1+(int)(volume->Size[2]/Resolution) );
      
      output->SetDirectionY( 0, 0, -1 );
      
      if ( Direction == SCANDIRECTION_DORSAL_VENTRAL ) 
	{
	output->SetOrigin( 0, Position, volume->Size[2] );
	output->SetDirectionX( 1, 0, 0 );
	} 
      else
	{
	output->SetOrigin( volume->Size[0], Position, volume->Size[2] );
	output->SetDirectionX( -1, 0, 0 );      
	}
      
      break;
    }
}

} // namespace cmtk
