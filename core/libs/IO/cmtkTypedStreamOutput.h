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

#ifndef __cmtkTypedStreamOutput_h_included_
#define __cmtkTypedStreamOutput_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>

#include <IO/cmtkTypedStream.h>

#include <stack>
#include <stdio.h>

#include <zlib.h>

#ifndef NULL
#define NULL 0
#endif

#include <string>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Class for writing "typedstream" archives.
 */
class TypedStreamOutput : public TypedStream
{
public:
  /// This class.
  typedef TypedStreamOutput Self;

  /// Base class.
  typedef TypedStream Superclass;

  /// Access modes for archives.
  typedef enum 
  {
    /// Currently unset.
    MODE_UNSET,
    /// Write-only access.
    MODE_WRITE,
    /// Write-only access piped through zlib/gzip compression.
    MODE_WRITE_ZLIB,
    /// Open existing archive and append to it.
    MODE_APPEND
  } Mode;
  
  /// Default constructor.
  TypedStreamOutput() : TypedStream() {}

  /** Open constructor.
   *\param filename Name of the archive to open.
   *\param mode Access mode, ie. read-only, write-only, etc.
   */
  TypedStreamOutput( const std::string& filename, const Self::Mode mode );

  /** Open constructor for separate path and archive names.
   *\param dir Directory to open archive in.
   *\param archive Name of the archive to open.
   *\param mode Access mode, ie. read-only, write-only, etc.
   */
  TypedStreamOutput( const std::string& dir, const std::string& archive, const Self::Mode mode );

  /** Destructor.
   * Close() is called to close a possibly open archive.
   */
  virtual ~TypedStreamOutput();

  /** Open another archive without constructing a new object.
   */
  void Open( const std::string& filename, const Self::Mode mode );

  /** Open another archive in explicit directory.
   */
  void Open( const std::string& dir, const std::string& archive, const Self::Mode mode );

  /** Close an open archive.
   */
  void Close();

  /** Begin a section.
   * This function will start a new section and increase the indentation level by one.
   *\param section Name of the new section.
   *\return Error condition.
   */
  Self::Condition Begin( const std::string& section = NULL );

  /** End a section.
   * In the open archive, this function will close the last section and 
   * decrease the nesting level by one.
   *\param flush If this flag is set, the output file buffer will be flushed
   * after closing the section.
   *\return Error condition.
   */
  Self::Condition End( const bool flush = false );

  /// Write a boolean value to an open archive.
  Self::Condition WriteBool( const char* key /*!< The name of the field under which to write this value in the archive.*/, 
				  const bool value /*!< Value to write to the archive under the given key. */ );

  /// Write an integer value to an open archive.
  Self::Condition WriteInt( const char* key /*!< The name of the field under which to write this value in the archive.*/, 
				 const int value /*!< Value to write to the archive under the given key. */ );

  /// Write a float value to an open archive.
  Self::Condition WriteFloat( const char* key /*!< The name of the field under which to write this value in the archive.*/, 
				   const float value /*!< Value to write to the archive under the given key. */ );

  /// Write a double precision float value to an open archive.
  Self::Condition WriteDouble( const char* key /*!< The name of the field under which to write this value in the archive.*/, 
				    const double value /*!< Value to write to the archive under the given key. */ );

  /// Write an Types::Coordinate value to an open archive.
  Self::Condition WriteCoordinate( const char* key /*!< The name of the field under which to write this value in the archive.*/, 
					const Types::Coordinate value /*!< Value to write to the archive under the given key. */ ) 
  {
#ifdef CMTK_COORDINATES_FLOAT
    return this->WriteFloat( key, value );
#else
    return this->WriteDouble( key, value );
#endif
  }
  
  /// Write an Types::DataItem value to an open archive.
  Self::Condition WriteItem( const char* key /*!< The name of the field under which to write this value in the archive.*/, 
				  const Types::DataItem value /*!< Value to write to the archive under the given key. */ ) 
  {
#ifdef CMTK_DATA_FLOAT
    return this->WriteFloat( key, value );
#else
    return this->WriteDouble( key, value );
#endif
  }
  
  /// Write a string to an open archive.
  Self::Condition WriteString( const char* key /*!< The name of the field under which to write this string in the archive.*/, 
				    const char* value /*!< String to write to the archive under the given key. */ );

  /// Write a string to an open archive.
  Self::Condition WriteString( const char* key /*!< The name of the field under which to write this string in the archive.*/, 
				    const std::string& value /*!< String to write to the archive under the given key. */ );

  /** Write array of integer values to an open archive.
   */
  Self::Condition WriteIntArray( const char* key /*!< The name of the field under which to write this array in the archive.*/, 
				      const int* array /*!< Pointer to the array to be written.*/, 
				      const int size /*!< Number of values in the array. This is the number of values written to the archive. */,
				      const int valuesPerLine = 10 /*!< Optional number of values per line of text written to the archive. This improves readability of the resulting archive as a text. */ );

  /** Write array of binay encoded boole values to an open archive.
   */
  Self::Condition WriteBoolArray( const char* key /*!< The name of the field under which to write this array in the archive.*/, 
				       const byte* array /*!< Pointer to the array to be written.*/, 
				       const int size /*!< Number of values in the array. This is the number of values written to the archive. */, 
				       const int valuesPerLine = 10 /*!< Optional number of values per line of text written to the archive. This improves readability of the resulting archive as a text. */ );

  /** Write array of single-precision values to an open archive.
   */
  Self::Condition WriteFloatArray( const char* key/*!< The name of the field under which to write this array in the archive.*/, 
					const float* array /*!< Pointer to the array to be written.*/, 
					const int size /*!< Number of values in the array. This is the number of values written to the archive. */,
					const int valuesPerLine = 10 /*!< Optional number of values per line of text written to the archive. This improves readability of the resulting archive as a text. */ );

  /** Write array of double-precision values to an open archive.
   */
  Self::Condition WriteDoubleArray( const char* key /*!< The name of the field under which to write this array in the archive.*/, 
					 const double* array /*!< Pointer to the array to be written.*/, 
					 const int size /*!< Number of values in the array. This is the number of values written to the archive. */,
					 const int valuesPerLine = 10 /*!< Optional number of values per line of text written to the archive. This improves readability of the resulting archive as a text. */ );

  /** Write array of double- or single precision values to an open archive.
   * Whether double- or single-precision data is written depends on the 
   * definition of the CMTK_COORDINATES_DOUBLE preprocessor symbol. This function
   * is thus guaranteed to always match the Types::Coordinate type.
   *\see CMTK_COORDINATES_DOUBLE
   *\see Types::Coordinate
   */
  Self::Condition WriteCoordinateArray( const char* key/*!< The name of the field under which to write this array in the archive.*/, 
					     const Types::Coordinate* array /*!< Pointer to the array to be written.*/, 
					     const int size /*!< Number of values in the array. This is the number of values written to the archive. */,
					     const int valuesPerLine = 10 /*!< Optional number of values per line of text written to the archive. This improves readability of the resulting archive as a text. */ )
  { 
#ifdef CMTK_COORDINATES_DOUBLE
    return this->WriteDoubleArray( key, array, size, valuesPerLine );
#else
    return this->WriteFloatArray( key, array, size, valuesPerLine );
#endif
  }
  
  /** Write array of double- or single precision values to an open archive.
   * Whether double- or single-precision data is written depends on the 
   * definition of the CMTK_DATA_DOUBLE preprocessor symbol. This function
   * is thus guaranteed to always match the Types::DataItem type.
   *\see CMTK_DATA_DOUBLE
   *\see Types::DataItem
   */
  Self::Condition WriteItemArray( const char* key /*!< The name of the field under which to write this array in the archive.*/, 
				  const Types::DataItem* array /*!< Pointer to the array to be written.*/, 
				  const int size /*!< Number of values in the array. This is the number of values written to the archive. */, 
				  const int valuesPerLine = 10 /*!< Optional number of values per line of text written to the archive. This improves readability of the resulting archive as a text. */ )
  { 
#ifdef CMTK_DATA_DOUBLE
    return this->WriteDoubleArray( key, array, size, valuesPerLine );
#else
    return this->WriteFloatArray( key, array, size, valuesPerLine );
#endif
  }

private:
  /// Mode the current archive was opened with.
  Self::Mode m_Mode;

};

//@}

} // namespace cmtk

//@}

#endif // #ifndef __cmtkTypedStreamOutput_h_included_
