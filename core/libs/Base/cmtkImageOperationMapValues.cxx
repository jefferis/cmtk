/*
//
//  Copyright 2010 Torsten Rohlfing
//
//  Copyright 2011, 2012 SRI International
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

#include "cmtkImageOperationMapValues.h"

#include <limits>

cmtk::ImageOperationMapValues::ImageOperationMapValues( const char* mapping, const bool exclusive )
  : m_Exclusive( exclusive )
{
  const char* rptr = mapping;

  while ( rptr )
    {
    const char* comma = strchr( rptr, ',' );
    const char* plus = strchr( rptr, '+' );
    
    double value;
    std::vector<Types::DataItem> fromValues;
    while ( comma && ( !plus || comma < plus) )
      {
      if ( 1 == sscanf( rptr, "%20lf", &value ) )
	fromValues.push_back( value );
      rptr = comma+1;
      comma = strchr( rptr, ',' );
      }
    
    double newValue;
    if ( 2 == sscanf( rptr, "%20lf:%20lf", &value, &newValue ) )
      {
      fromValues.push_back( value );
      
      for ( size_t i = 0; i < fromValues.size(); ++i )
	{
	this->m_Mapping[fromValues[i]] = newValue;
	}
      }
    else
      {
      if ( 1 == sscanf( rptr, "%20lf", &value ) )
	{
	fromValues.push_back( value );
	
	for ( size_t i = 0; i < fromValues.size(); ++i )
	  {
	  this->m_Mapping[fromValues[i]] = std::numeric_limits<double>::quiet_NaN();
	  }
	}
      else
	{
	StdErr << "ERROR: could not parse mapping\n\t" << mapping << "\nwhich is supposed to be VAL0[,VAL1,...][:NEWVAL]\n";
	}
      }

    // if more rules are following, separated by ";", then move to next one, otherwise terminate loop.
    if ( plus )
      {
      rptr = plus + 1;
      plus = strchr( rptr, '+' );
      }
    else
      rptr = NULL;
    }
}

cmtk::UniformVolume::SmartPtr
cmtk::ImageOperationMapValues::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  TypedArray& volumeData = *(volume->GetData());
#pragma omp parallel for
  for ( int i = 0; i < static_cast<int>( volumeData.GetDataSize() ); ++i )
    {
    Types::DataItem value = 0;
    if ( volumeData.Get( value, i ) )
      {
      std::map<Types::DataItem,Types::DataItem>::const_iterator it = this->m_Mapping.find( value );
      if ( it != this->m_Mapping.end() )
	{
	const Types::DataItem newValue = it->second;
	if ( finite( newValue ) )
	  volumeData.Set( newValue, i );
	else
	  volumeData.SetPaddingAt( i );
	}
      else
	{
	// value not explicitly mapped; see if we're in "exclusive" mode.
	if ( this->m_Exclusive )
	  {
	  volumeData.SetPaddingAt( i );
	  }
	}
      }
    }
  
  return volume;
}
