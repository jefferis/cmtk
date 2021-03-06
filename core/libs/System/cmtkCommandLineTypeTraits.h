/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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
#include <cmtkconfig.h>

#include <string>
#include <sstream>

namespace
cmtk
{

/** \addtogroup System */
//@{

template<class T>
class CommandLineTypeTraitsBase
{
public:
  /// Convert a value of this type to string.
  static std::string ValueToString( const T& value )
  {
    std::ostringstream stream;
    stream << value;
    return stream.str();
  }

  /// Convert a value of this type to string with minimal added markup (for XML output).
  static std::string ValueToStringMinimal( const T& value )
  {
    std::ostringstream stream;
    stream << value;
    return stream.str();
  }
};

/// Template for traits to handle command line arguments of different types.
template<class T>
class
CommandLineTypeTraits
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<T>
{
public:
  /// Return name of the parameter type (for XML).
  static const char* GetName() 
  { 
    return "none"; 
  }
};

template<>
class 
CommandLineTypeTraits<const char*>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<const char*>
{
public:
  static const char* GetName() 
  { 
    return "string";
  }

  static std::string ValueToString( const char*& value )
  {
    std::ostringstream stream;
    if ( value )
      stream << "\"" << value << "\"";
    else
      stream << "NONE";
    return stream.str();
  }

  static std::string ValueToStringMinimal( const char*& value )
  {
    std::ostringstream stream;
    if ( value )
      stream << value;
    return stream.str();
  }
};

template<>
class 
CommandLineTypeTraits<std::string>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<std::string>
{
public:
  static const char* GetName() 
  { 
    return "string";
  }

  static std::string ValueToString( const std::string& value )
  {
    std::ostringstream stream;
    if ( value.length() )
      stream << "\"" << value << "\"";
    else
      stream << "NONE";
    return stream.str();
  }

  static std::string ValueToStringMinimal( const std::string& value )
  {
    std::ostringstream stream;
    stream << value;
    return stream.str();
  }
};

template<>
class 
CommandLineTypeTraits< std::vector<std::string> >
/// Inherit generic template members
  : public CommandLineTypeTraitsBase< std::vector<std::string> >
{
public:
  static const char* GetName() 
  { 
    return "vector<string>";
  }

  static std::string ValueToString( const std::vector<std::string>& value )
  {
    std::ostringstream stream;
    for ( size_t i = 0; i < value.size(); ++i )
      stream << value[i] << " ";
    return stream.str();
  }

  /// Convert a value of this type to string with minimal added markup (for XML output).
  static std::string ValueToStringMinimal( const std::vector<std::string>& value )
  {
    return ValueToString( value );
  }
};

template<>
class 
CommandLineTypeTraits<int>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<int>
{
public:
  static const char* GetName() 
  { 
    return "integer";
  }
};

template<>
class 
CommandLineTypeTraits<unsigned int>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<unsigned int>
{
public:
  static const char* GetName() 
  { 
    return "integer";
  }
};

template<>
class 
CommandLineTypeTraits<short>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<short>
{
public:
  static const char* GetName() 
  { 
    return "integer";
  }
};

template<>
class 
CommandLineTypeTraits<unsigned short>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<unsigned short>
{
public:
  static const char* GetName() 
  { 
    return "integer";
  }
};

template<>
class 
CommandLineTypeTraits<signed char>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<signed char>
{
public:
  /// Convert a value of this type to numerical string.
  static std::string ValueToString( const signed char& value )
  {
    std::ostringstream stream;
    stream << static_cast<int>( value );
    return stream.str();
  }

  /// Convert a value of this type to numerical string with minimal added markup (for XML output).
  static std::string ValueToStringMinimal( const signed char& value )
  {
    std::ostringstream stream;
    stream << static_cast<int>( value );
    return stream.str();
  }

  static const char* GetName() 
  { 
    return "integer";
  }
};

template<>
class 
CommandLineTypeTraits<unsigned char> :
  /// Inherit generic template members
  public CommandLineTypeTraitsBase<unsigned char>
{
public:
  /// Convert a value of this type to numerical string.
  static std::string ValueToString( const unsigned char& value )
  {
    std::ostringstream stream;
    stream << static_cast<int>( value );
    return stream.str();
  }

  /// Convert a value of this type to numerical string with minimal added markup (for XML output).
  static std::string ValueToStringMinimal( const unsigned char& value )
  {
    std::ostringstream stream;
    stream << static_cast<int>( value );
    return stream.str();
  }

  static const char* GetName() 
  { 
    return "integer";
  }
};

template<>
class 
CommandLineTypeTraits<float>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<float>
{
public:
  static const char* GetName()
  { 
    return "float"; 
  }
};

template<>
class 
CommandLineTypeTraits<double>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<double>
{
public:
  static const char* GetName() 
  { 
    return "double"; 
  }
};

template<>
class 
CommandLineTypeTraits<bool>
/// Inherit generic template members
  : public CommandLineTypeTraitsBase<bool>
{
public:
  static const char* GetName() 
  { 
    return "boolean"; 
  }
};

} // namespace cmtk
