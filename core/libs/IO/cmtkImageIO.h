/*
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkImageIO_h_included_
#define __cmtkImageIO_h_included_

#include <cmtkconfig.h>

#include <cmtkImageInfo.h>
#include <cmtkFileFormat.h>

#include <cmtkStudy.h>
#include <cmtkScalarImage.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Constants for image format capabilities.
typedef enum {
  /// Image files provide image dimension (number of pixels).
  IMAGEFORMAT_DIMS = 1,
  /// Image files provide data structure (bytes per pixel, data type).
  IMAGEFORMAT_STRUCTURE = 2,
  /// Image files provide image calibration (pixel size, slice position).
  IMAGEFORMAT_CALIBRATION = 4
} ImageFormatCapability;

/** Virtual base class for image readers/writers.
 * This class provides a unified access to all supported (2D) image formats.
 * The image data read by an object of this class is internally allocated.
 */
class ImageIO 
{
public:
  /** Create a new image reader/writer.
   *@param format The name of the image file format to instance a reader/writer
   * for. Permitted values so far are "DICOM" and "PGM". For historical
   * reasons, "ACR-NEMA" is also recognized, but the created object is actually
   * a DICOM handler.
   *@return The newly created instace of the image reader/writer handling the
   * given format. If this format is not supported, NULL is returned.
   *@see DICOM
   *@see PGM
   */
  static ImageIO* Create( const char* format );

  /** Create a new image reader/writer.
   *@param format The ID of the image file format to instance a reader/writer
   * for.
   *@return The newly created instace of the image reader/writer handling the
   * given format. If this format is not supported, NULL is returned.
   */
  static ImageIO* Create( const FileFormatID format );

  /** Virtual destructor.
   * The FreeDataPtr() member function is called to free storage allocated by
   * this instance.
   *@see #FreeDataPtr
   */
  virtual ~ImageIO() { this->FreeDataPtr(); }

  /** Read an image.
   * A pointer to the image data is stored in the "DataPtr" field after return
   * from thsi function. It can then be retrieved by a call to 
   * GetReleaseDataPtr().
   */
  virtual ScalarImage* Read( const char*, const Study* = NULL, const int = 0 ) const
  { return NULL; }


  /** Write an image.
   * Derived classes must override this function with an implementation that
   * writes the image data held by this instance to a file of the respective
   * format.
   *@param filename Name (full path) of the file to write to. If only a
   * relative filename is given, it is interpreted relative to the current
   * working directory. There is no guarantee that the working directory
   * remains unchanged between subsequent calls to this function.
   *@param imageInfo Geometry information corresponding to the image data.
   *@param anonymize If this optional parameter is non-zero, no
   * patient-specific information is written to the image file. This parameter
   * defaults to zero, resulting in all information stored to be written.
   */
  virtual void Write ( const char*, const ImageInfo&, const int = 0 ) {};

  /** Return the pointer to the image data.
   * Using this function, the data pointer remains under control of the
   * reader/writer object. It is therefore read-only. For takeover of storage
   * control, use GetReleaseDataPtr() instead.
   */
  const void* GetDataPtr() const { return DataPtr; }

  /** Return and release pointer to image data.
   * This function returns the current pointer to the image data. Also, it
   * releases the pointer into the control of the caller. The reader/writer
   * object will neither free the storage afterwards nor access it in any other
   * way.
   */
  void* GetReleaseDataPtr();

  /** Explicitly set data pointer.
   * This function can be used to transfer control over an existing image data
   * pointer to this object. For example the DICOM library allocates storage
   * for image data itself. By handing the pointer to this area over to the
   * ImageIO class, additional allocation and copying of data are avoided.
   * It must be ensured however, that the creator of the respective memory
   * region releases control over the pointer.
   *@param dataPtr Pointer to the existing image data.
   *@param size Size of the allocated region in data items (pixels). The number
   * of bytes per pixel is defined by the next parameter.
   *@param itemSize This parameter gives the number of bytes per item in the
   * array pointed to by "dataPtr". It defaults to 1, thus if it is omitted
   * the parameter "size" gives the exact size in bytes.
   */
  void SetDataPtr( void *const dataPtr, const size_t size, const size_t itemSize = 1 );

  /** Link data pointer to externally controlled image.
   * This function works in a similar way as 'SetDataPtr()' does. However, 
   * using 'LinkDataPtr()' the control over the memory holding the image data
   * remains with the caller. Thus, it will not be freed by this object.
   */
  void LinkDataPtr( void *const dataPtr, const size_t size, const size_t itemSize = 1 );

  /** Get the size of the allocated image data.
   *@param itemSize This parameter, defaulting to 1, gives the number of bytes
   * stored per pixel of the current image. The allocated size is divided by
   * this value to determine the number of pixels. Omitting this parameters
   * yields the allocated size in bytes.
   *@return The number of data items (bytes by default) stored in the currently
   * allocated image.
   */
  size_t GetSize( const size_t itemSize = 1 ) const;

  /** Return status of last read or write operation.
   * This function may and should be called after calls to Read() or Write()
   * to determine whether the respective operation was completed successfully.
   *@return 1 if the previous read or write operation completed successfully,
   * 0 if an error occurred. In the latter case GetErrorMsg can be called to
   * retrieve the textual description of the error.
   @*see #GetErrorMsg
   */
  int GetStatus() const { return ! Error; }

  /** Get error message from previous operation.
   * In case the last read or write operation failed (Success() returned zero),
   * this function will provide a textual description of what caused the
   * failure.
   *@return A pointer to a string describing the last error, or NULL if the
   * last operation completed succesfully.
   */
  const char* GetErrorMsg() const;

  /// Return flags indicating which information is natively.
  virtual byte GetFormatCapabilities() const { return 0; }

protected:
  /** Default constructor.
   * Initialize all fields with zero values.
   */
  ImageIO();

  /** Set error status.
   * In case of an error, derived classes can call this function to set the
   * inherited error flag.
   *@param error The internal error number. Zero represents successful 
   * completion of the last operation. This parameter defaults to 1. So far,
   * there is no distinction between different non-zero error codes.
   */
  void SetError( const int error = 1 ) { Error = error; }

  /** Set error description.
   * This function is defined with a variable number of arguments so callers
   * can give an arbitrary number of additional parameters. These are handled
   * in the same way the printf() function does.
   *@param format Format string for the remaining parameters in printf() style.
   * If no additional parameters are given, this parameter is the actual error
   * message.
   */
  void SetErrorMsg( const char* format, ... );

  /// The current error status.
  int Error;

  /// Storage for the current error message.
  char ErrorMsg[2048];

  /** Allocate memory for image data.
   *@param size The number of data items to allocate. The number of bytes per
   * item is defined by the next parameter.
   *@param itemSize This parameter gives the number of bytes required to hold
   * each of the "size" items in the newly allocated image data array. It
   * defaults to 1, ie. one byte is exactly one item (pixel).
   */
  void* AllocDataPtr( const size_t size, const size_t itemSize = 1 );

  /** Free image data memory.
   * If the current instance has allocated and still holds image data, it is
   * freed by this function.
   */
  void FreeDataPtr();

  /// Pointer to the actual image data.
  void *DataPtr;

  /// Size of the memory region in bytes (!) pointed to by "DataPtr".
  size_t Size;

  /** Data pointer control flag.
   * If this flag is non-zero, the memory region pointed to by 'DataPtr' is
   * controlled by this object and can be freed if necessary.
   */
  bool FreeDataPtrFlag;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageIO_h_included_
