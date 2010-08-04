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

#ifndef __cmtkNumericTraits_h_included_
#define __cmtkNumericTraits_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkTypes.h"
#include "Base/cmtkMathUtil.h"

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Base class for numeric traits definition for primitive data types..
template<class TType>
class NumericTraits
{
public:
  /// Default value to use for padding data with a given type.
  static const TType DefaultPaddingValue;

  /// Convert from float to given data type.
  static TType ConvertFromDataItem( const Types::DataItem value );
};

template<>
class NumericTraits<char>
{
public:
  static const char DefaultPaddingValue = -1;

  static char ConvertFromDataItem( const Types::DataItem value )
  {
    using namespace std;
    if ( MathUtil::IsNaN( value ) )
      {
      return DefaultPaddingValue;
      }
    else
      {
      return static_cast<char>( value + 0.5 );
      }
  }
};

template<>
class NumericTraits<unsigned char>
{
public:
  static const unsigned char DefaultPaddingValue = 255;

  static unsigned char ConvertFromDataItem( const Types::DataItem value )
  {
    using namespace std;
    if ( MathUtil::IsNaN( value ) )
      {
      return DefaultPaddingValue;
      }
    else
      {
      return static_cast<unsigned char>( value + 0.5 );
      }
  }
};

template<>
class NumericTraits<short>
{
public:
  static const short DefaultPaddingValue = -32768;

  static short ConvertFromDataItem( const Types::DataItem value )
  {
    using namespace std;
    if ( MathUtil::IsNaN( value ) )
      {
      return DefaultPaddingValue;
      }
    else
      {
      return static_cast<short>( value + 0.5 );
      }
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkNumericTraits_h_included_
