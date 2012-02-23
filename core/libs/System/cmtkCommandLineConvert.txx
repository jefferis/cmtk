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

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Convert string to bool.
template<> 
inline bool CommandLine::Item::Convert<bool>( const char* str ) 
{ 
  if ( ! str ) 
    return false; 
  else
    return !strcmp( str, "yes" ); 
}

/// Convert string to float.
template<> 
inline float CommandLine::Item::Convert<float>( const char* str ) 
{
  return static_cast<float>( this->ConvertStrToDouble( str ) );
}

/// Convert string to double.
template<> 
inline double CommandLine::Item::Convert<double>( const char* str ) 
{
  return static_cast<double>( this->ConvertStrToDouble( str ) );
}

/// Convert string to long int.
template<> 
inline long int
CommandLine::Item::Convert<long int>( const char* str ) 
{
  return static_cast<long int>( this->ConvertStrToLong( str ) );
}

/// Convert string to int.
template<> 
inline int CommandLine::Item::Convert<int>( const char* str ) 
{
  return static_cast<int>( this->ConvertStrToLong( str ) );
}

/// Convert string to unsigned int.
template<> 
inline unsigned int 
CommandLine::Item::Convert<unsigned int>( const char* str ) 
{
  return static_cast<unsigned int>( this->ConvertStrToLong( str ) );
}

/// Convert string to unsigned long
template<> 
inline unsigned long int
CommandLine::Item::Convert<unsigned long int>( const char* str ) 
{
  return static_cast<unsigned long int>( this->ConvertStrToLong( str ) );
}

/// Convert string to char.
template<> 
inline char
CommandLine::Item::Convert<char>( const char* str ) 
{
  return static_cast<char>( this->ConvertStrToLong( str ) );
}

/// Convert string to unsigned char.
template<> 
inline unsigned char
CommandLine::Item::Convert<unsigned char>( const char* str ) 
{ 
  return static_cast<unsigned char>( this->ConvertStrToLong( str ) );
}

/// Convert string to signed char.
template<> 
inline signed char
CommandLine::Item::Convert<signed char>( const char* str ) 
{ 
  return static_cast<signed char>( this->ConvertStrToLong( str ) );
}

/// Convert string to short integer.
template<> 
inline short 
CommandLine::Item::Convert<short>( const char* str ) 
{
  return static_cast<short>( this->ConvertStrToLong( str ) );
}

/// Convert string to unsigned short integer.
template<> 
inline unsigned short 
CommandLine::Item::Convert<unsigned short>( const char* str ) 
{
  return static_cast<unsigned short>( this->ConvertStrToLong( str ) );
}

/// Convert string to string.
template<> 
inline
const char* 
CommandLine::Item::Convert<const char*>( const char* str ) 
{
  return str; 
}

/// Convert string to vector.
template<> 
inline
std::vector<std::string> 
CommandLine::Item::Convert< std::vector<std::string> >( const char* str ) 
{
  std::vector<std::string> v;
  v.push_back( std::string( str ) );
  return v; 
}

} // namespace cmtk
