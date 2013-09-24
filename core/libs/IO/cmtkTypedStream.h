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

/** base class for reading and writing of "typedstream" archives.
 */
class TypedStream 
{
public:
  /// This class.
  typedef TypedStream Self;

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

  /// Set debugging flag.
  void SetDebugFlag( const Self::DebugFlag debugFlag = Self::DEBUG_ON /*!< Set the debug flag to this value. */ )
  { 
    this->m_DebugFlag = debugFlag;
  }
  
protected:
  /// Internal: Length of the read buffer for one archive line.  
  static const int LIMIT_BUFFER = 1024;

  /// Pointer to the actual file.
  FILE *File;

  /// Pointer to the compressed file in decompression mode.
  gzFile GzFile;

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
