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

#include <Base/cmtkTypes.h>

#include <string.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

const char *DataClassString[] = 
{ "grey", "label", "unknown",
  NULL };

const char* DataTypeName[] = 
{ "byte (8bit unsigned)",
  "char (8bit signed)",
  "short (16bit signed)",
  "ushort (16bit unsigned)",
  "int (32bit signed)",
  "uint (32bit unsigned)",
  "float (32bit)",
  "double (64bit)"
};

DataClass
StringToDataClass( const char *dataClassStr )
{
  if ( dataClassStr ) 
    {
    for ( int idx=0; DataClassString[idx]; ++idx ) 
      {
      if ( !strcmp( dataClassStr, DataClassString[idx] ) )
	return (DataClass) idx;
      }
    }
  
  return DATACLASS_UNKNOWN;
}

const char*
DataClassToString( const DataClass dataClass )
{
  return DataClassString[ static_cast<int>( dataClass ) ];
}

size_t
TypeItemSize ( const ScalarDataType dtype ) 
{
  switch (dtype)
    {
    case TYPE_BYTE: return sizeof(byte);      
    case TYPE_CHAR: return sizeof(char);
    case TYPE_SHORT: return sizeof(short);
    case TYPE_USHORT: return sizeof(unsigned short);
    case TYPE_INT: return sizeof(int);
    case TYPE_UINT: return sizeof(int);
    case TYPE_FLOAT: return sizeof(float);
    case TYPE_DOUBLE: return sizeof(double);
    default : return 0;
    }
}

ScalarDataType 
SelectDataTypeInteger( const byte itemSize, const bool signBit )
{
  if ( signBit ) 
    {
    switch ( itemSize ) 
      {
      case 1 : return TYPE_CHAR;
      case 2 : return TYPE_SHORT;
      case 4 : return TYPE_INT;
      default: return TYPE_NONE;
      }
    } 
  else
    {
    switch ( itemSize ) 
      {
      case 1 : return TYPE_BYTE;
      case 2 : return TYPE_USHORT;
      case 4 : return TYPE_INT;
      default: return TYPE_NONE;
      }
    }
}

ScalarDataType
GetSignedDataType( const ScalarDataType dtype )
{
  switch ( dtype ) 
    {
    case TYPE_BYTE:
      return TYPE_CHAR;
    case TYPE_USHORT:
      return TYPE_SHORT;
    case TYPE_UINT:
      return TYPE_INT;
    default:
      return dtype;
    }
}

ScalarDataType
GetUnsignedDataType( const ScalarDataType dtype )
{
  switch ( dtype ) 
    {
    case TYPE_CHAR:
      return TYPE_BYTE;
    case TYPE_SHORT:
      return TYPE_USHORT;
    case TYPE_INT:
      return TYPE_UINT;
    default:
      return dtype;
    }
}

} // namespace cmtk
