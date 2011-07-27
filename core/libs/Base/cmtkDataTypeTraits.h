/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkDataTypeTraits_h_included_
#define __cmtkDataTypeTraits_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkMathUtil.h>

#include <math.h>
#include <stdlib.h>
#include <limits>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Data type traits */
template<class TType>
class DataTypeTraits
{
public:
  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_NONE; }

  static const ScalarDataType DataTypeID = TYPE_NONE;
};

/** Data type traits for single-precision floating point. */
template<>
class DataTypeTraits<float>
{
public:
  /** This class. */
  typedef DataTypeTraits<float> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_FLOAT; }

  static const ScalarDataType DataTypeID = TYPE_FLOAT;

  /** Get absolute value. */
  static inline float Abs( const float value ) 
  {
    return fabsf( value );
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline float Convert ( const T value, const bool = false, const float = 0 ) 
  { 
    return static_cast<float>( value );
  }

  /** Return padding data value (i.e. Inf) for the given type. */
  static inline float ChoosePaddingValue () 
  { 
    return std::numeric_limits<float>::infinity();
  }
};

/** Data type traits for double-precision floating point. */
template<>
class DataTypeTraits<double>
{
public:
  /** This class. */
  typedef DataTypeTraits<double> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_DOUBLE; }

  static const ScalarDataType DataTypeID = TYPE_DOUBLE;

  /** Get absolute value. */
  static inline double Abs( const double value ) 
  {
    return fabs( value );
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline double Convert ( const T value, const bool = false, const double = 0 ) 
  { 
    return static_cast<double>( value );
  }
  
  /** Return padding data value (i.e. Inf) for the given type. */
  static inline double ChoosePaddingValue () 
  { 
    return std::numeric_limits<double>::infinity();
  }
};

/** Data type traits for unsigned char (byte). */
template<>
class DataTypeTraits<byte>
{
public:
  /** This class. */
  typedef DataTypeTraits<byte> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_BYTE; }

  static const ScalarDataType DataTypeID = TYPE_BYTE;

  /** Minimum value in this type. */
  static const Types::DataItem Min;

  /** Maximum value in this type. */
  static const Types::DataItem Max;

  /** Get absolute value. */
  static inline byte Abs( const byte value ) 
  {
    return abs( value );
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline byte Convert ( const T value, const bool paddingFlag = false, const byte paddingData = 0 ) 
  { 
    using namespace std;
    if ( MathUtil::IsFinite( value ) )
      {
      return (byte) ((value<Self::Min) ? Self::Min : (value+0.5>Self::Max) ? Self::Max : floor(value+0.5));
      }
    else
      {
      if ( paddingFlag )
	return paddingData;
      else
	return ChoosePaddingValue();
      }
  }
  
  /** Return padding data value for the given type. */
  static inline byte ChoosePaddingValue () 
  { 
    return 255;
  }
};

/** Data type traits for signed char. */
template<>
class DataTypeTraits<char>
{
public:
  /** This class. */
  typedef DataTypeTraits<char> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_CHAR; }

  static const ScalarDataType DataTypeID = TYPE_CHAR;

  /** Minimum value in this type. */
  static const Types::DataItem Min;

  /** Maximum value in this type. */
  static const Types::DataItem Max;

  /** Get absolute value. */
  static inline char Abs( const char value ) 
  {
    return abs( value );
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline char Convert ( const T value, const bool paddingFlag = false, const char paddingData = 0 ) 
  { 
    using namespace std;
    if ( MathUtil::IsFinite( value ) )
      {
      return (char) ((value<Self::Min) ? Self::Min : (value+0.5>Self::Max) ? Self::Max : floor(value+0.5));
      }
    else
      {
      if ( paddingFlag )
	return paddingData;
      else
	return ChoosePaddingValue();
      }
  }
  
  /** Return padding data value for the given type. */
  static inline char ChoosePaddingValue () 
  { 
    return -1;
  }
};

/** Data type traits for signed short. */
template<>
class DataTypeTraits<signed short>
{
public:
  /** This class. */
  typedef DataTypeTraits<signed short> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_SHORT; }

  static const ScalarDataType DataTypeID = TYPE_SHORT;

  /** Minimum value in this type. */
  static const Types::DataItem Min;

  /** Maximum value in this type. */
  static const Types::DataItem Max;

  /** Get absolute value. */
  static inline signed short Abs( const signed short value ) 
  {
    return abs( value );
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline signed short Convert ( const T value, const bool paddingFlag = false, const signed short paddingData = 0 ) 
  { 
    using namespace std;
    if ( MathUtil::IsFinite( value ) )
      {
      return (signed short) (((value<Self::Min) ? Self::Min : (value+0.5>Self::Max) ? Self::Max : floor(value+0.5)));
      }
    else
      {
      if ( paddingFlag )
	return paddingData;
      else
	return ChoosePaddingValue();
      }
  }
  
  /** Return padding data value for the given type. */
  static inline signed short ChoosePaddingValue () 
  { 
    return -1;
  }
};

/** Data type traits for unsigned short. */
template<>
class DataTypeTraits<unsigned short>
{
public:
  /** This class. */
  typedef DataTypeTraits<unsigned short> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_USHORT; }

  static const ScalarDataType DataTypeID = TYPE_USHORT;

  /** Minimum value in this type. */
  static const Types::DataItem Min;

  /** Maximum value in this type. */
  static const Types::DataItem Max;

  /** Get absolute value. */
  static inline unsigned short Abs( const unsigned short value ) 
  {
    return abs( value );
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline unsigned short Convert ( const T value, const bool paddingFlag = false, const unsigned short paddingData = 0 ) 
  { 
    using namespace std;
    if ( MathUtil::IsFinite( value ) )
      {
      return (unsigned short) ((value<Self::Min) ? Self::Min : (value+0.5>Self::Max) ? Self::Max : floor(value+0.5));
      }
    else
      {
      if ( paddingFlag )
	return paddingData;
      else
	return ChoosePaddingValue();
      }
  }
  
  /** Return padding data value for the given type. */
  static inline unsigned short ChoosePaddingValue () 
  { 
    return USHRT_MAX;
  }
};

/** Data type traits for int. */
template<>
class DataTypeTraits<int>
{
public:
  /** This class. */
  typedef DataTypeTraits<int> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_INT; }

  static const ScalarDataType DataTypeID = TYPE_INT;

  /** Minimum value in this type. */
  static const Types::DataItem Min;

  /** Maximum value in this type. */
  static const Types::DataItem Max;

  /** Get absolute value. */
  static inline int Abs( const int value ) 
  {
    return abs( value );
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline int Convert ( const T value, const bool paddingFlag = false, const int paddingData = 0 ) 
  { 
    using namespace std;
    if ( MathUtil::IsFinite( value ) )
      {
      return (int) ((value<Self::Min) ? Self::Min : (value+0.5>Self::Max) ? Self::Max : floor(value+0.5));
      }
    else
      {
      if ( paddingFlag )
	return paddingData;
      else
	return ChoosePaddingValue();
      }
  }
  
  /** Return padding data value for the given type. */
  static inline int ChoosePaddingValue () 
  { 
    return -1;
  }
};

/** Data type traits for unsigned int. */
template<>
class DataTypeTraits<unsigned int>
{
public:
  /** This class. */
  typedef DataTypeTraits<unsigned int> Self;

  /** Get the scalar type ID constant for this type. */
  static ScalarDataType GetScalarDataType() { return TYPE_UINT; }

  static const ScalarDataType DataTypeID = TYPE_UINT;

  /** Minimum value in this type. */
  static const Types::DataItem Min;

  /** Maximum value in this type. */
  static const Types::DataItem Max;

  /** Get absolute value. */
  static inline unsigned int Abs( const unsigned int value ) 
  {
    // usigned is always positive
    return value;
  }

  /** Return given value converted (and rounded) to discrete type. */
  template<class T>
  static inline unsigned int Convert ( const T value, const bool paddingFlag = false, const unsigned int paddingData = 0 ) 
  { 
    using namespace std;
    if ( MathUtil::IsFinite( value ) )
      {
      return (unsigned int) ((value<Self::Min) ? Self::Min : (value+0.5>Self::Max) ? Self::Max : floor(value+0.5));
      }
    else
      {
      if ( paddingFlag )
	return paddingData;
      else
	return ChoosePaddingValue();
      }
  }
  
  /** Return padding data value for the given type. */
  static inline unsigned int ChoosePaddingValue () 
  { 
    return static_cast<unsigned int>( -1 );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDataTypeTraits_h_included_
