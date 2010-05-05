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

#ifndef __cmtkTypedStream_h_included_
#define __cmtkTypedStream_h_included_

#include <cmtkconfig.h>

#include <cmtkTypes.h>

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

/**@name TypedStream.h
 */
//@{

/// Access modes for archives.
typedef enum 
{
  /// Read-only access.
  TYPEDSTREAM_READ,
  /// Write-only access.
  TYPEDSTREAM_WRITE,
  /// Write-only access piped through zlib/gzip compression.
  TYPEDSTREAM_WRITE_ZLIB,
  /// Open existing archive and append to it.
  TYPEDSTREAM_APPEND
} TypedStreamMode;

/// Condition upon function return.
typedef enum 
{
  /// There was an error.
  TYPEDSTREAM_ERROR,
  /// No error encountered; operation completed successfully.
  TYPEDSTREAM_OK
} TypedStreamCondition;

/// Classes of error conditions
typedef enum 
{
  /// No error.
  TYPEDSTREAM_ERROR_NONE,
  /// Unknown error.
  TYPEDSTREAM_ERROR_UNKNOWN,
  /** A call to a system function returned an error condition.
   * To find out more details, "errno" may be consulted.
   */
  TYPEDSTREAM_ERROR_SYSTEM,
  /** Error in the format of the open archive.
   */
  TYPEDSTREAM_ERROR_FORMAT,
  /** Wrong or invalid arguments given to a function.
   */
  TYPEDSTREAM_ERROR_ARG,
  /** The requested operation is not available in the current stream mode.
   * This usually means that a write access was requested on a read-only
   * archive or a Seek() operation on a write-only archive.
   */
  TYPEDSTREAM_ERROR_MODE,
  /** Error in a primitive data object.
   * A value in the archive does not have the correct syntax for the expected
   * type.
   */
  TYPEDSTREAM_ERROR_TYPE,
  /** An internal limit was exhausted.
   * As we are now using a proper STL stack for keeping track of open levels,
   * this condition should not occur any more.
   */
  TYPEDSTREAM_ERROR_LIMIT,
  /** Close of a level was requested when none was open.
   */
  TYPEDSTREAM_ERROR_LEVEL,
  /** The current stream is invalid.
   * This condition is set when an access is tried without opening a file 
   * first.
   */
  TYPEDSTREAM_ERROR_INVALID,
  TYPEDSTREAM_ERROR_MAX
} TypedStreamStatus;

/// Identifiers for supported primitive data types.
typedef enum 
{
  /// Interger.
  TYPEDSTREAM_TYPE_INT,
  /// Boolean (Yes/No).
  TYPEDSTREAM_TYPE_BOOL,
  /// Binary boolean (0/1).
  TYPEDSTREAM_TYPE_BINARYBOOL,
  /// Single-precision float.
  TYPEDSTREAM_TYPE_FLOAT,
  /// Double-precision float
  TYPEDSTREAM_TYPE_DOUBLE,
  /// String (char*).
  TYPEDSTREAM_TYPE_STRING
} TypedStreamType;

/// Identifiers for tokens in archives.
typedef enum 
{
  /// End-of-file.
  TYPEDSTREAM_EOF,
  /// Section beginning "{".
  TYPEDSTREAM_BEGIN,
  /// Section end "}".
  TYPEDSTREAM_END,
  /// Key (field name).
  TYPEDSTREAM_KEY,
  /// Field value.
  TYPEDSTREAM_VALUE,
  /// Comment.
  TYPEDSTREAM_COMMENT
} TypedStreamToken;

/// Debug flag values.
typedef enum 
{
  /// There was an error.
  TYPEDSTREAM_DEBUG_OFF,
  /// No error encountered; operation completed successfully.
  TYPEDSTREAM_DEBUG_ON
} TypedStreamDebugFlag;

/// Internal: Length of the read buffer for one archive line.  
#define TYPEDSTREAM_LIMIT_BUFFER	1024

/// Constant for use with flush parameter of End() member function.
#define TYPEDSTREAM_FLUSH true

/// Constant for use with flush parameter of End() member function.
#define TYPEDSTREAM_NOFLUSH false

/** Class for reading and writing og "typedstream" archives.
 * This class provides the same functions as DHZB's old, C-based "typedstream"
 * library. The interface has been remodelled, especially primitive read
 * functions have been made more convenient to use. Also, the implementation of
 * certain features has been made clearer, eg. using an explicit stack to keep
 * track of the currently open archive levels.
 */
class TypedStream 
{
public:
  /// Default constructor.
  TypedStream();

  /** Open constructor.
   *@param name Name of the archive to open.
   *@param mode Access mode, ie. read-only, write-only, etc.
   */
  TypedStream( const char *filename, const TypedStreamMode mode );

  /** Open constructor for separate path and archive names.
   *@param dir Directory to open archive in.
   *@param archive Name of the archive to open.
   *@param mode Access mode, ie. read-only, write-only, etc.
   */
  TypedStream( const char *dir, const char *archive, const TypedStreamMode mode );

  /** Destructor.
   * Close() is called to close a possibly open archive.
   */
  virtual ~TypedStream();

  /** Open another archive without constructing a new object.
   */
  void Open( const char *filename, const TypedStreamMode mode );

  /** Open another archive in explicit directory.
   */
  void Open( const char *dir, const char *archive, const TypedStreamMode mode );

  /** Close an open archive.
   */
  void Close();

  /** Move to a particular section in the open archive.
   * The named section is found if it is either inside the currently open
   * section or after it on the same level.
   *
   * This function may only be called for read-only archive, ie. for such that
   * were opened in TYPEDSTREAM_READONLY mode. For writeable archive, it 
   * will return an error.
   */
  TypedStreamCondition Seek( const char *section, const bool forward = false );

  /** Rewind archive.
   * This function resets filepointer of an open archive to the beginning of
   * the current section.
   */
  TypedStreamCondition Rewind();

  /** Return validity of archive.
   *@return 1 if an archive is currently open, 0 if not.
   */
  int IsValid() { return (File != NULL) || (GzFile != NULL); }

  /** Return status of last operation.
   */
  TypedStreamStatus GetStatus() const { return Status; }

  /** Begin a section.
   * In an archive opened for writing, this function will start a new section
   * and increase the indentation level by one. For a read-only archive, this
   * function will generate an error condition.
   *@param section Name of the new section.
   *@return Error condition.
   */
  TypedStreamCondition Begin( const char* section = NULL );

  /** End a section.
   * In the open archive, this function will close the last section and 
   * decrease the nesting level by one.
   *@param flush If this flag is set, the output file buffer will be flushed
   * after closing the section.
   *@return Error condition.
   */
  TypedStreamCondition End( const bool flush = false );

  /** Read boolean value from an open archive.
   * This function recognizes both yes/no and 0/1 entries in the archive.
   * First, "yes" and "no" is tried, if that doesn't work the function reads
   * an integer value from the same key.
   *@param key The name of the boolean entry in the archive.
   *@param defaultValue Default value returned if no valid entry can be read. 
   * This parameter can be omitted and defaults to 0.
   *@return If reading was succesful, the value from the archive is returned.
   * Otherwise the value given as the "defaultValue" parameter is returned.
   */
  bool ReadBool( const char* key, const bool defaultValue = false, const bool forward = false );

  /** Read array of boole values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  TypedStreamCondition ReadBoolArray( const char* key, byte *const array, const int size = 0, const bool forward = false );
  
  /** Read integer value from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  int ReadInt( const char* key, const int defaultValue = 0, const bool forward = false );

  /** Read array of integer values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  TypedStreamCondition ReadIntArray( const char* key, int *const array, const int size = 0, const bool forward = false );

  /** Read single-precision value from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  float ReadFloat( const char* key, const float defaultValue = 0, const bool forward = false );

  /** Read array of single-precision values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  TypedStreamCondition ReadFloatArray( const char* key, float *const array, const int size = 0, const bool forward = false );

  /** Read double-precision value from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  double ReadDouble( const char* key, const double defaultValue = 0, const bool forward = false );

  /** Read array of double-precision values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  TypedStreamCondition ReadDoubleArray( const char* key, double *const array, const int size = 0, const bool forward = false );

  /** Read double- or single precision value from an open archive.
   * Whether double- or single-precision data is read depends on the definition
   * of the CMTK_COORDINATES_DOUBLE preprocessor symbol. This function is thus
   * guaranteed to always match the Types::Coordinate type.
   *@see CMTK_COORDINATES_DOUBLE
   *@see Types::Coordinate
   */
  Types::Coordinate ReadCoordinate( const char* key, const Types::Coordinate defaultValue = 0, const bool forward = false ) 
  {
#ifdef CMTK_COORDINATES_DOUBLE
    return this->ReadDouble( key, defaultValue, forward );
#else
    return this->ReadFloat( key, defaultValue, forward );
#endif
  }
  
  /** Read double- or single precision value from an open archive.
   * Whether double- or single-precision data is read depends on the definition
   * of the CMTK_DATA_DOUBLE preprocessor symbol. This function is thus
   * guaranteed to always match the Types::DataItem type.
   *@see CMTK_DATA_DOUBLE
   *@see Types::DataItem
   */
  Types::DataItem ReadItem( const char* key, const Types::DataItem defaultValue = 0, const bool forward = false ) 
  {
#ifdef CMTK_DATA_DOUBLE
    return this->ReadDouble( key, defaultValue, forward );
#else
    return this->ReadFloat( key, defaultValue, forward );
#endif
  }
  
  /** Read array of double- or single precision values from an open archive.
   * Whether double- or single-precision data is read depends on the definition
   * of the CMTK_COORDINATES_DOUBLE preprocessor symbol. This function is thus
   * guaranteed to always match the Types::Coordinate type.
   *@see CMTK_COORDINATES_DOUBLE
   *@see Types::Coordinate
   */
  TypedStreamCondition ReadCoordinateArray
( const char* key, Types::Coordinate *const array, const int size = 0, const bool forward = false ) 
  {
#ifdef CMTK_COORDINATES_DOUBLE
    return this->ReadDoubleArray( key, array, size, forward );
#else
    return this->ReadFloatArray( key, array, size, forward );
#endif
  }

  /** Read array of double- or single precision values from an open archive.
   * Whether double- or single-precision data is read depends on the definition
   * of the CMTK_DATA_DOUBLE preprocessor symbol. This function is thus
   * guaranteed to always match the Types::DataItem type.
   *@see CMTK_DATA_DOUBLE
   *@see Types::DataItem
   */
  TypedStreamCondition ReadItemArray( const char* key, Types::DataItem *const array, const int size = 0, const bool forward = false ) 
  {
#ifdef CMTK_DATA_DOUBLE
    return this->ReadDoubleArray( key, array, size, forward );
#else
    return this->ReadFloatArray( key, array, size, forward );
#endif
  }

  /** Read null-terminated string from an open archive.
   * The string returned is newly allocated by this function. So unless NULL
   * is returner, the string must later be freed by the caller in order to
   * avoid memory leaks.
   *@return A pointer to a newly allocated string is returned if reading was
   * succesful. If no valid entry could be read from the archive, a copy of
   * the string given as "defaultValue" parameter is returned. If that 
   * parameter was NULL, the same value is also returned.
   */
  char* ReadString( const char* key, const char* defaultValue = NULL, const bool forward = false );

  /// Write a boolean value to an open archive.
  TypedStreamCondition WriteBool( const char* key, const bool value );

  /// Write an integer value to an open archive.
  TypedStreamCondition WriteInt( const char* key, const int value );

  /// Write a float value to an open archive.
  TypedStreamCondition WriteFloat( const char* key, const float value );

  /// Write a double precision float value to an open archive.
  TypedStreamCondition WriteDouble( const char* key, const double value );

  /// Write an Types::Coordinate value to an open archive.
  TypedStreamCondition WriteCoordinate( const char* key, const Types::Coordinate value ) 
  {
#ifdef CMTK_COORDINATES_FLOAT
    return this->WriteFloat( key, value );
#else
    return this->WriteDouble( key, value );
#endif
  }
  
  /// Write an Types::DataItem value to an open archive.
  TypedStreamCondition WriteItem( const char* key, const Types::DataItem value ) 
  {
#ifdef CMTK_DATA_FLOAT
    return this->WriteFloat( key, value );
#else
    return this->WriteDouble( key, value );
#endif
  }
  
  /// Write a string to an open archive.
  TypedStreamCondition WriteString( const char* key, const char* value );

  /// Write a string to an open archive.
  TypedStreamCondition WriteString( const char* key, const std::string& value );

  /// Write a formated comment to an open archive.
  TypedStreamCondition WriteComment( const char* fmt, ... );

  /// Write string array as comment to an open archive.
  TypedStreamCondition WriteComment( const int argc, const char* argv[] );

  /// Write string array as comment to an open archive.
  TypedStreamCondition WriteComment( int argc, char* argv[] );

  /** Write array of integer values to an open archive.
   */
  TypedStreamCondition WriteIntArray( const char* key, const int* array, const int size, const int valuesPerLine = 10  );

  /** Write array of binay encoded boole values to an open archive.
   */
  TypedStreamCondition WriteBoolArray( const char* key, const byte* array, const int size, const int valuesPerLine = 10 );

  /** Write array of single-precision values to an open archive.
   */
  TypedStreamCondition WriteFloatArray( const char* key, const float* array, const int size, const int valuesPerLine = 10 );

  /** Write array of double-precision values to an open archive.
   */
  TypedStreamCondition WriteDoubleArray( const char* key, const double* array, const int size, const int valuesPerLine = 10 );

  /** Write array of double- or single precision values to an open archive.
   * Whether double- or single-precision data is written depends on the 
   * definition of the CMTK_COORDINATES_DOUBLE preprocessor symbol. This function
   * is thus guaranteed to always match the Types::Coordinate type.
   *@see CMTK_COORDINATES_DOUBLE
   *@see Types::Coordinate
   */
  TypedStreamCondition WriteCoordinateArray
  ( const char* key, const Types::Coordinate* array, const int size, const int valuesPerLine = 10  )
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
   *@see CMTK_DATA_DOUBLE
   *@see Types::DataItem
   */
  TypedStreamCondition WriteItemArray( const char* key, const Types::DataItem* array, const int size, const int valuesPerLine = 10  )
  { 
#ifdef CMTK_DATA_DOUBLE
    return this->WriteDoubleArray( key, array, size, valuesPerLine );
#else
    return this->WriteFloatArray( key, array, size, valuesPerLine );
#endif
  }

  /// Set debugging flag.
  void SetDebugFlag( const TypedStreamDebugFlag debugFlag = TYPEDSTREAM_DEBUG_ON ) { DebugFlag = debugFlag; }
  
private:
  /** Initialize internal data structures.
   * This function is called from both constructors to initialize the internal
   * data structures of this object.
   */
  void InitInternals();

  /// Pointer to the actual file.
  FILE *File;

  /// Pointer to the compressed file in decompression mode.
  gzFile GzFile;

  /// Mode the current archive was opened with.
  TypedStreamMode Mode;

  /** Holds the status of the last operation.
   */
  TypedStreamStatus Status;

  /** Number of significant digits for "float" fields.
   * As all float numbers are written to the archive as strings, part of the
   * native resolution is lost. This field determines, how many significant
   * digits are preserved when converting single precision float numbers to
   * strings.
   */
  int PrecisionFloat;

  /** Number of significant digits for "double" fields.
   * As all float numbers are written to the archive as strings, part of the
   * native resolution is lost. This field determines, how many significant
   * digits are preserved when converting double precision float numbers to
   * strings.
   */
  int PrecisionDouble;

  /// Buffer for the current line read from the archive.
  char Buffer[TYPEDSTREAM_LIMIT_BUFFER];

  /// Pointer to the "key" part of the line currently in Buffer.
  char *BufferKey;

  /// Pointer to the "value" part of the line currently in Buffer.
  char *BufferValue;

  /** Stack of open section levels.
   * This stack holds the starting positions of all currently open sections.
   * The entries are byte positions relative to the beginning of the file.
   */
  std::stack<int> LevelStack;

  /** Utility function: Read an array of arbitrary type.
   * This function is called by all reader functions. Internally, a "switch"
   * statement selects the correct code for the effective data type to be read.
   * Besides, common functions such as the skipping of inserted sections are
   * implemented as shared code for all data types.
   *\param forward If this flag is set, then the search will be forward only,
   * i.e., the file pointer will NOT be reset to the beginning of the current 
   * section.
   */
  TypedStreamCondition GenericReadArray
  ( const char * key, const int type, void *const array, const int arraySize, const bool forward = false );

  /// Read the next archive line to the buffer.
  TypedStreamToken ReadLineToken();

  /** Compare two strings.
   * Other than the standard library's strcmp() function, this implementation
   * ignores upper and lowercase. Also, strings are terminated by either NULL
   * characters or any white space or newline.
   *@return 0 for identical strings (up to upper-/lowercase), 1 for 
   * non-identical strings.
   */
  static int StringCmp( const char * s1, const char * s2 );

  /** Separate next token.
   * This function identifies the next token in the given string, sets a NULL
   * character to mark its end and returns a pointer to the token's first
   * character. The state between calls is saved in the "SplitPosition" field.
   * Calling the function with NULL as a parameter resets the internal state.
   */
  char* StringSplit( char * s1 ) const;

  /// Internal position pointer for "StringSplit()".
  mutable char *SplitPosition;

  /// Return the identifier for the generated archive format (version).
  static const char* GetTypedStreamIdent() 
  {
    return "! TYPEDSTREAM 1.1\n";
  }
  
  /// Debug flag.
  TypedStreamDebugFlag DebugFlag;
  
  /// Output diagnostic message if debug flag is set.
  void DebugOutput( const char* format, ... );
};

//@}

} // namespace cmtk

//@}

#endif // #ifndef __cmtkTypedstream_h_included_
