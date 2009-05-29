/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkVtkLookupTable.h>

void
cmtk::VtkLookupTable::SetFromStudy( const Study* study )
{
  Types::DataItem range = study->GetMaximumValue() - study->GetMinimumValue();
  if ( range < 1024 ) {
    this->SetNumberOfTableValues( (int)range + 1 );
  } else {
    this->SetNumberOfTableValues( 1024 );
  }

  switch ( study->GetStandardColormap() ) {
  case PALETTE_GRAY : 
    this->SetHueRange( 0, 0 );
    this->SetSaturationRange( 0, 0 );
    if ( study->GetReverseColormap() )
      this->SetValueRange( 1, 0 );
    else
      this->SetValueRange( 0, 1 );
    break;
  case PALETTE_RED : 
    this->SetHueRange( 0, 0 );
    this->SetSaturationRange( 1, 1 );
    if ( study->GetReverseColormap() )
      this->SetValueRange( 1, 0 );
    else
      this->SetValueRange( 0, 1 );
    break;
  case PALETTE_GREEN : 
    this->SetHueRange( 0.33, 0.33 );
    this->SetSaturationRange( 1, 1 );
    if ( study->GetReverseColormap() )
      this->SetValueRange( 1, 0 );
    else
      this->SetValueRange( 0, 1 );
    break;
  case PALETTE_BLUE : 
    this->SetHueRange( 0.66, 0.66 );
    this->SetSaturationRange( 1, 1 );
    if ( study->GetReverseColormap() )
      this->SetValueRange( 1, 0 );
    else
      this->SetValueRange( 0, 1 );
    break;
  case PALETTE_RAINBOW : 
    if ( study->GetReverseColormap() )
      this->SetHueRange( 0, 0.66 );
    else
      this->SetHueRange( 0.66, 0 );
    this->SetSaturationRange( 1, 1 );
    this->SetValueRange( 1, 1 );
    break;
  }
}
