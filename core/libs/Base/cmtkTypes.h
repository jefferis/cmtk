/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkTypes_h_included_
#define __cmtkTypes_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkMathUtil.h"

#include <stdlib.h>
#include <limits.h>
#include <cfloat>
#include <cstddef>

#ifdef HAVE_VALUES_H
#  include <values.h>
#endif

#ifdef _MSC_VER
typedef unsigned short ushort;
#endif

#ifndef NULL
#define NULL 0
#endif

#ifndef _HAVE_BYTE_
typedef unsigned char byte;
#define _HAVE_BYTE_
#endif // #ifndef _HAVE_BYTE_

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Identifiers for coordinate axes.
 *\note It is very important that the enumeration constants remain defined
 * equal to the integer constants 0 through 2, since quite a bit of code
 * depends on this for array indexing.
 */
enum 
{
  /// x-axis.
  AXIS_X = 0,
  /// y-axis.
  AXIS_Y = 1,
  /// z-axis.
  AXIS_Z = 2
};

/// Class of image data.
typedef enum 
{
  /// Grey-level data.
  DATACLASS_GREY,
  /// (Segmented) label data.
  DATACLASS_LABEL,
  /// Data type unknown.
  DATACLASS_UNKNOWN
} DataClass;

/// Convert string to data class identifier.
DataClass StringToDataClass( const char *dataClassStr );

/// Convert data class identifier to string.
const char* DataClassToString( const DataClass dataClass );

/// Scalar data type identifiers.
typedef enum 
{
  /// Unsigned byte data (8 bit, range 0-255).
  TYPE_BYTE = 0,
  /// Signed byte data (8 bit, range -128-127).
  TYPE_CHAR = 1,
  /// Signed byte data (16 bit, range -32768-32767).
  TYPE_SHORT = 2,
  /// Unsigned byte data (16 bit, range 0-65535).
  TYPE_USHORT = 3,
  /// Signed integer data (32 bit).
  TYPE_INT = 4,
  /// Unsigned integer data (32 bits).
  TYPE_UINT = 5,
  /// Single precision float data (32 bits).
  TYPE_FLOAT = 6,
  /// Double precision float data (64 bits).
  TYPE_DOUBLE = 7,
  /// No data type defined.
  TYPE_NONE = -1
} ScalarDataType;

/// Names of scalar data types.
extern const char* DataTypeName[];

#ifdef CMTK_DATA_FLOAT
namespace Types
{ 
/** @memo Definition of the data exchange item type
 * All data retrievals, stores and conversions are done using this type.
 */
typedef float DataItem; 
}
const ScalarDataType TYPE_ITEM = TYPE_FLOAT;
#define CMTK_ITEM_MAX FLT_MAX
#define CMTK_ITEM_MIN FLT_MIN
#define CMTK_ITEM_NAN CMTK_FLOAT_NAN
#else
namespace Types 
{
typedef double DataItem; 
}
const ScalarDataType TYPE_ITEM = TYPE_DOUBLE;
#define CMTK_ITEM_MAX DBL_MAX
#define CMTK_ITEM_MIN DBL_MIN
#define CMTK_ITEM_NAN CMTK_DOUBLE_NAN
#endif // #ifdef CMTK_DATA_FLOAT

namespace Types
{

/// Range of DataItem values specified as lower and upper bound.
template<class T>
class Range
{
public:
  /// Default constructor: do nothing at all.
  Range() {}

  /// Constructor.
  Range( const T& lowerBound, const T& upperBound ) : m_LowerBound( lowerBound ), m_UpperBound( upperBound ) {}

  /// Conversion constructor.
  template<class T2>
  explicit Range( const Range<T2>& range ) : m_LowerBound( range.m_LowerBound ), m_UpperBound( range.m_UpperBound ) {}

  /// Compute "width" of range, i.e., upper minus lower bound.
  T Width() const
  {
    return this->m_UpperBound - this->m_LowerBound;
  }

  /// Lower bound.
  T m_LowerBound;

  /// Upper bound.
  T m_UpperBound;
};

/// Convenience declaration: range of DataItem values.
typedef Range<DataItem> DataItemRange;

}

#ifdef CMTK_COORDINATES_FLOAT
/** @memo Definition of the coordinate data type
 * All spatial locations, distances, etc. are stored using this type.
 */
namespace Types 
{ 
typedef float Coordinate; 
}
const ScalarDataType TYPE_COORDINATE = TYPE_FLOAT;
#else
/// Define float type used for coordinates.
namespace Types 
{ 
typedef double Coordinate; 
}
const ScalarDataType TYPE_COORDINATE = TYPE_DOUBLE;
#endif

/// Return item size for given scalar data type.
size_t TypeItemSize ( const ScalarDataType dtype );

/// Select integer data type based on item size and sign bit.
ScalarDataType SelectDataTypeInteger( const byte itemSize,const bool signBit );

/// Return signed datatype ID corresponding to given datatype.
ScalarDataType GetSignedDataType( const ScalarDataType dtype );

/// Return difference datatype ID for given pair of datatypes.
ScalarDataType GetDifferenceDataType( const ScalarDataType dtype1, const ScalarDataType dtype2 );

/// Return unsigned datatype ID corresponding to given datatype.
ScalarDataType GetUnsignedDataType( const ScalarDataType dtype );

namespace
Types
{

/// Template for traits class to combine two different real (floating point) data types.
template<class T1,class T2>
class Combined
{
};

/// Combination of two single-precision floating point values.
template<>
class Combined<float,float>
{
public:
  /// Single-precision floating point.
  typedef float Type;
};

/// Combination of single- and double-precision floating point values.
template<>
class Combined<float,double>
{
public:
  /// Double-precision floating point. 
 typedef double Type;
};

/// Combination of double- and single-precision floating point values.
template<>
class Combined<double,float>
{
public:
  /// Double-precision floating point. 
  typedef double Type;
};

/// Combination of two double-precision floating point values.
template<>
class Combined<double,double>
{
public:
  /// Double-precision floating point. 
  typedef double Type;
};

} // namespace Types

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypes_h_included_
