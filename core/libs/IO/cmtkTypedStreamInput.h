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

#ifndef __cmtkTypedStreamInput_h_included_
#define __cmtkTypedStreamInput_h_included_

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

//@{

/** Class for reading "typedstream" archives.
 */
class TypedStreamInput : public TypedStream
{
public:
  /// This class.
  typedef TypedStreamInput Self;

  /// Base class.
  typedef TypedStream Superclass;

  /// Default constructor.
  TypedStreamInput() : TypedStream() {}

  /** Open constructor.
   *\param filename Name of the archive to open.
   *\param mode Access mode, ie. read-only, write-only, etc.
   */
  TypedStreamInput( const std::string& filename );

  /** Open constructor for separate path and archive names.
   *\param dir Directory to open archive in.
   *\param archive Name of the archive to open.
   */
  TypedStreamInput( const std::string& dir, const std::string& archive );

  /** Destructor.
   * Close() is called to close a possibly open archive.
   */
  virtual ~TypedStreamInput();

  /** Open another archive without constructing a new object.
   */
  void Open( const std::string& filename );

  /** Open another archive in explicit directory.
   */
  void Open( const std::string& dir, const std::string& archive );

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

  /** Move to beginning of section.
   * This function will set the file read pointer to the beginning of the current section
   *\return Error condition.
   */
  Self::Condition Begin();

  /** Close current section.
   * In the open archive, this function will close the last section and decrease the nesting level by one.
   *\return Error condition.
   */
  Self::Condition End();

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

  /** Read STL string from an open archive.
   *\return The string that was read. If no valid entry could be read from the archive, a copy of
   * the string given as "defaultValue" parameter is returned.
   */
  std::string ReadStdString( const char* key /*!< The name of the field in the archive.*/, 
			     const std::string& defaultValue = "" /*!< Default value returned if the field is not found in the archive. */, 
			     const bool forward = false /*!< Flag: read forward from current position in stream (if false, reset to current section start) */);

  /// Get open stream major version.
  int GetReleaseMajor() const
  {
    if ( this->m_Status == Self::ERROR_NONE )
      {
      return this->m_ReleaseMajor;
      }
    else
      {
      return -1;
      }
  }

  /// Get open stream minor version.
  int GetReleaseMinor() const
  {
    if ( this->m_Status == Self::ERROR_NONE )
      {
      return this->m_ReleaseMinor;
      }
    else
      {
      return -1;
      }
  }

private:
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

  /// Release major version.
  int m_ReleaseMajor;

  /// Release minor version.
  int m_ReleaseMinor;
};

//@}

} // namespace cmtk

//@}

#endif // #ifndef __cmtkTypedStreamInput_h_included_
