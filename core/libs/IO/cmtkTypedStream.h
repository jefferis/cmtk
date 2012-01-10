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

#ifndef __cmtkTypedStream_h_included_
#define __cmtkTypedStream_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>

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

/**\name TypedStream.h
 */
//@{

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
  /// This class.
  typedef TypedStream Self;

  /// Access modes for archives.
  typedef enum 
  {
    /// Currently unset.
    MODE_UNSET,
    /// Read-only access.
    MODE_READ,
    /// Write-only access.
    MODE_WRITE,
    /// Write-only access piped through zlib/gzip compression.
    MODE_WRITE_ZLIB,
    /// Open existing archive and append to it.
    MODE_APPEND
  } Mode;
  
  /// Condition upon function return.
  typedef enum 
  {
    /// There was an error.
    CONDITION_ERROR,
    /// No error encountered; operation completed successfully.
    CONDITION_OK
  } Condition;
  
  /// Classes of error conditions
  typedef enum 
  {
    /// No error.
    ERROR_NONE,
    /// Unknown error.
    ERROR_UNKNOWN,
    /** A call to a system function returned an error condition.
     * To find out more details, "errno" may be consulted.
     */
    ERROR_SYSTEM,
    /** Error in the format of the open archive.
     */
    ERROR_FORMAT,
    /** Wrong or invalid arguments given to a function.
     */
    ERROR_ARG,
    /** The requested operation is not available in the current stream mode.
     * This usually means that a write access was requested on a read-only
     * archive or a Seek() operation on a write-only archive.
     */
    ERROR_MODE,
    /** Error in a primitive data object.
     * A value in the archive does not have the correct syntax for the expected
     * type.
     */
    ERROR_TYPE,
    /** An internal limit was exhausted.
     * As we are now using a proper STL stack for keeping track of open levels,
     * this condition should not occur any more.
     */
    ERROR_LIMIT,
    /** Close of a level was requested when none was open.
     */
    ERROR_LEVEL,
    /** The current stream is invalid.
     * This condition is set when an access is tried without opening a file 
     * first.
     */
    ERROR_INVALID,
    ERROR_MAX
  } Status;
  
  /// Identifiers for supported primitive data types. 
  typedef enum 
  {
    /// Interger.
    TYPE_INT,
    /// Boolean (Yes/No).
    TYPE_BOOL,
    /// Binary boolean (0/1).
    TYPE_BINARYBOOL,
    /// Single-precision float.
    TYPE_FLOAT,
    /// Double-precision float
    TYPE_DOUBLE,
    /// String (char*).
    TYPE_STRING
  } Type;

  /// Identifiers for tokens in archives.
  typedef enum 
  {
    /// End-of-file.
    TOKEN_EOF,
    /// Section beginning "{".
    TOKEN_BEGIN,
    /// Section end "}".
    TOKEN_END,
    /// Key (field name).
    TOKEN_KEY,
    /// Field value.
    TOKEN_VALUE,
    /// Comment.
    TOKEN_COMMENT
  } Token;
  
  /// Debug flag values.
  typedef enum 
  {
    /// There was an error.
    DEBUG_OFF,
    /// No error encountered; operation completed successfully.
    DEBUG_ON
  } DebugFlag;
  
  /// Default constructor.
  TypedStream();

  /** Open constructor.
   *\param filename Name of the archive to open.
   *\param mode Access mode, ie. read-only, write-only, etc.
   */
  TypedStream( const char* filename, const Self::Mode mode );

  /** Open constructor for separate path and archive names.
   *\param dir Directory to open archive in.
   *\param archive Name of the archive to open.
   *\param mode Access mode, ie. read-only, write-only, etc.
   */
  TypedStream( const char* dir, const char* archive, const Self::Mode mode );

  /** Destructor.
   * Close() is called to close a possibly open archive.
   */
  virtual ~TypedStream();

  /** Open another archive without constructing a new object.
   */
  void Open( const char* filename, const Self::Mode mode );

  /** Open another archive in explicit directory.
   */
  void Open( const char* dir, const char* archive, const Self::Mode mode );

  /** Close an open archive.
   */
  void Close();

  /** Move to a particular section in the open archive.
   * The named section is found if it is either inside the currently open
   * section or after it on the same level.
   *
   * This function may only be called for read-only archive, ie. for such that
   * were opened in MODE_READONLY mode. For writeable archive, it 
   * will return an error.
   */
  Self::Condition Seek( const char* section /*!< Name of the section whose beginning stream pointer is moved to. */, 
			     const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

  /** Rewind archive.
   * This function resets filepointer of an open archive to the beginning of
   * the current section.
   */
  Self::Condition Rewind();

  /** Return validity of archive.
   *\return 1 if an archive is currently open, 0 if not.
   */
  int IsValid() 
  {
    return (this->File != NULL) || (this->GzFile != NULL); 
  }

  /** Return status of last operation.
   */
  Self::Status GetStatus() const 
  { 
    return this->m_Status; 
  }

  /** Begin a section.
   * In an archive opened for writing, this function will start a new section
   * and increase the indentation level by one. For a read-only archive, this
   * function will generate an error condition.
   *\param section Name of the new section.
   *\return Error condition.
   */
  Self::Condition Begin( const char* section = NULL );

  /** End a section.
   * In the open archive, this function will close the last section and 
   * decrease the nesting level by one.
   *\param flush If this flag is set, the output file buffer will be flushed
   * after closing the section.
   *\return Error condition.
   */
  Self::Condition End( const bool flush = false );

  /** Read boolean value from an open archive.
   * This function recognizes both yes/no and 0/1 entries in the archive.
   * First, "yes" and "no" is tried, if that doesn't work the function reads
   * an integer value from the same key.
   *\return If reading was succesful, the value from the archive is returned.
   * Otherwise the value given as the "defaultValue" parameter is returned.
   */
  bool ReadBool( const char* key /*!< The name of the boolean entry in the archive.*/, 
		 const bool defaultValue = false /*!< Default value returned if no valid entry can be read. This parameter can be omitted and defaults to false.*/,
		 const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

  /** Read array of boole values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  Self::Condition ReadBoolArray( const char* key /*!< The name of the array in the archive.*/, 
				      byte *const array /*!< Pointer to allocated storage for the array to be read into.*/, 
				      const int size /*!< Size of the array.*/, 
				      const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);
  
  /** Read integer value from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  int ReadInt( const char* key /*!< The name of the field in the archive.*/, 
	       const int defaultValue = 0 /*!< Default value returned if no valid entry can be read. This parameter can be omitted and defaults to zero.*/, 
	       const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

  /** Read array of integer values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  Self::Condition ReadIntArray( const char* key /*!< The name of the array in the archive.*/, 
				     int *const array /*!< Pointer to allocated storage for the array to be read into.*/, 
				     const int size /*!< Size of the array.*/, 
				     const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);
  
  /** Read single-precision value from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  float ReadFloat( const char* key /*!< The name of the field in the archive.*/, 
		   const float defaultValue = 0 /*!< Default value returned if no valid entry can be read. This parameter can be omitted and defaults to zero.*/, 
		   const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

  /** Read array of single-precision values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  Self::Condition ReadFloatArray( const char* key /*!< The name of the array in the archive.*/, 
				       float *const array /*!< Pointer to allocated storage for the array to be read into.*/, 
				       const int size /*!< Size of the array.*/, 
				       const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

  /** Read double-precision value from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  double ReadDouble( const char* key /*!< The name of the field in the archive.*/,
		     const double defaultValue = 0 /*!< Default value returned if the field is not found in the archive. */, 
		     const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

  /** Read array of double-precision values from an open archive.
   * For a description of parameters and return value see ReadBool.
   */
  Self::Condition ReadDoubleArray( const char* key /*!< The name of the array in the archive.*/, 
					double *const array /*!< Pointer to allocated storage for the array to be read into.*/, 
					const int size /*!< Size of the array.*/, 
					const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);
  
  /** Read double- or single precision value from an open archive.
   * Whether double- or single-precision data is read depends on the definition
   * of the CMTK_COORDINATES_DOUBLE preprocessor symbol. This function is thus
   * guaranteed to always match the Types::Coordinate type.
   *\see CMTK_COORDINATES_DOUBLE
   *\see Types::Coordinate
   */
  Types::Coordinate ReadCoordinate( const char* key /*!< The name of the field in the archive.*/, 
				    const Types::Coordinate defaultValue = 0 /*!< Default value if the field is not found.*/, 
				    const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */) 
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
   *\see CMTK_DATA_DOUBLE
   *\see Types::DataItem
   */
  Types::DataItem ReadItem( const char* key /*!< The name of the field in the archive.*/, 
			    const Types::DataItem defaultValue = 0 /*!< Default value returned if the field is not found in the archive. */, 
			    const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */) 
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
   *\see CMTK_COORDINATES_DOUBLE
   *\see Types::Coordinate
   */
  Self::Condition ReadCoordinateArray( const char* key /*!< The name of the array in the archive.*/, 
					    Types::Coordinate *const array /*!< Pointer to allocated storage for the array to be read into.*/, 
					    const int size /*!< Size of the array.*/, 
					    const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */)
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
   *\see CMTK_DATA_DOUBLE
   *\see Types::DataItem
   */
  Self::Condition ReadItemArray( const char* key /*!< The name of the array in the archive.*/, 
				      Types::DataItem *const array /*!< Pointer to allocated storage for the array to be read into.*/, 
				      const int size /*!< Size of the array.*/, 
				      const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */) 
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
   *\return A pointer to a newly allocated string is returned if reading was
   * succesful. If no valid entry could be read from the archive, a copy of
   * the string given as "defaultValue" parameter is returned. If that 
   * parameter was NULL, the same value is also returned.
   */
  char* ReadString( const char* key /*!< The name of the field in the archive.*/, 
		    const char* defaultValue = NULL /*!< Default value returned if the field is not found in the archive. */, 
		    const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

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

  /// Write a formated comment to an open archive.
  Self::Condition WriteComment( const char* fmt /*!< printf-style format string for the remaining variable number of parameters */, ... );

  /// Write string array as comment to an open archive.
  Self::Condition WriteComment( const int argc /*!< Number of strings in the array. */, const char* argv[] /*!< Array of string pointers. */ );

  /// Write string array as comment to an open archive.
  Self::Condition WriteComment( int argc /*!< Number of strings in the array. */, char* argv[] /*!< Array of string pointers. */ );

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

  /// Set debugging flag.
  void SetDebugFlag( const Self::DebugFlag debugFlag = Self::DEBUG_ON /*!< Set the debug flag to this value. */ )
  { 
    this->m_DebugFlag = debugFlag;
  }
  
private:
  /// Internal: Length of the read buffer for one archive line.  
  static const int LIMIT_BUFFER = 1024;

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
  Self::Mode m_Mode;

  /** Holds the status of the last operation.
   */
  Self::Status m_Status;

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
  char Buffer[Self::LIMIT_BUFFER];

  /// Pointer to the "key" part of the line currently in Buffer.
  char* BufferKey;

  /// Pointer to the "value" part of the line currently in Buffer.
  char* BufferValue;

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
   */
  Self::Condition GenericReadArray( const char* key /*!< Field key (name)*/, 
				    const int type /*!< Array data type ID */, 
				    void *const array /*!< Target storage space for read data */, 
				    const int arraySize /*!< Number of array elements */, 
				    const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */ );
  
  /// Read the next archive line to the buffer.
  Self::Token ReadLineToken();
  
  /** Compare two strings.
   * Other than the standard library's strcmp() function, this implementation
   * ignores upper and lowercase. Also, strings are terminated by either NULL
   * characters or any white space or newline.
   *\return 0 for identical strings (up to upper-/lowercase), 1 for 
   * non-identical strings.
   */
  static int StringCmp( const char* s1, const char* s2 );

  /** Separate next token.
   * This function identifies the next token in the given string, sets a NULL
   * character to mark its end and returns a pointer to the token's first
   * character. The state between calls is saved in the "SplitPosition" field.
   * Calling the function with NULL as a parameter resets the internal state.
   */
  char* StringSplit( char* s1 /*!< String to split into tokens. */ ) const;

  /// Internal position pointer for "StringSplit()".
  mutable char* SplitPosition;

  /// Return the identifier for the generated archive format (version).
  static const char* GetTypedStreamIdent() 
  {
    return "! TYPEDSTREAM 1.1\n";
  }
  
  /// Debug flag.
  Self::DebugFlag m_DebugFlag;
  
  /// Output diagnostic message if debug flag is set.
  void DebugOutput( const char* format /*!< printf-style format string for the remaining variable number of function arguments.*/, ... );
};

//@}

} // namespace cmtk

//@}

#endif // #ifndef __cmtkTypedstream_h_included_
