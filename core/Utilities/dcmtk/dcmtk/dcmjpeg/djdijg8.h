/*
 *
 *  Copyright (C) 1997-2005, OFFIS
 *
 *  This software and supporting documentation were developed by
 *
 *    Kuratorium OFFIS e.V.
 *    Healthcare Information and Communication Systems
 *    Escherweg 2
 *    D-26121 Oldenburg, Germany
 *
 *  THIS SOFTWARE IS MADE AVAILABLE,  AS IS,  AND OFFIS MAKES NO  WARRANTY
 *  REGARDING  THE  SOFTWARE,  ITS  PERFORMANCE,  ITS  MERCHANTABILITY  OR
 *  FITNESS FOR ANY PARTICULAR USE, FREEDOM FROM ANY COMPUTER DISEASES  OR
 *  ITS CONFORMITY TO ANY SPECIFICATION. THE ENTIRE RISK AS TO QUALITY AND
 *  PERFORMANCE OF THE SOFTWARE IS WITH THE USER.
 *
 *  Module:  dcmjpeg
 *
 *  Author:  Norbert Olges, Marco Eichelberg
 *
 *  Purpose: decompression routines of the IJG JPEG library configured for 8 bits/sample.
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 16:59:24 $
 *  Source File:      $Source: /share/dicom/cvs-depot/dcmtk/dcmjpeg/include/dcmtk/dcmjpeg/djdijg8.h,v $
 *  CVS/RCS Revision: $Revision: 1.4 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */

#ifndef DJDIJG8_H
#define DJDIJG8_H

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmjpeg/djdecabs.h" /* for class DJDecoder */

extern "C"
{
  struct jpeg_decompress_struct;
}

class DJCodecParameter;

/** this class encapsulates the decompression routines of the
 *  IJG JPEG library configured for 8 bits/sample.
 */
class DJDecompressIJG8Bit : public DJDecoder
{
public:

  /** constructor
   *  @param cp codec parameters
   *  @param isYBR flag indicating if DICOM photometric interpretation is YCbCr
   */
  DJDecompressIJG8Bit(const DJCodecParameter& cp, OFBool isYBR);

  /// destructor
  virtual ~DJDecompressIJG8Bit();

  /** initializes internal object structures.
   *  Must be called before a new frame is decompressed.
   *  @return EC_Normal if successful, an error code otherwise
   */
  virtual OFCondition init();

  /** suspended decompression routine. Decompresses a JPEG frame
   *  until finished or out of data. Can be called with new data
   *  until a frame is complete.
   *  @param compressedFrameBuffer pointer to compressed input data, must not be NULL
   *  @param compressedFrameBufferSize size of buffer, in bytes
   *  @param uncompressedFrameBuffer pointer to uncompressed output data, must not be NULL.
   *     This buffer must not change between multiple decode() calls for a single frame.
   *  @param uncompressedFrameBufferSize size of buffer, in bytes (!)
   *     Buffer must be large enough to contain a complete frame.
   *  @param isSigned OFTrue, if uncompressed pixel data is signed, OFFalse otherwise
   *  @return EC_Normal if successful, EC_Suspend if more data is needed, an error code otherwise.
   */
  virtual OFCondition decode(
    Uint8 *compressedFrameBuffer,
    Uint32 compressedFrameBufferSize,
    Uint8 *uncompressedFrameBuffer,
    Uint32 uncompressedFrameBufferSize,
    OFBool isSigned);

  /** returns the number of bytes per sample that will be written when decoding.
   */
  virtual Uint16 bytesPerSample() const
  {
    return sizeof(Uint8);
  }

  /** after successful compression,
   *  returns the color model of the decompressed image
   */
  virtual EP_Interpretation getDecompressedColorModel() const
  {
    return decompressedColorModel;
  }

  /** callback function used to report warning messages and the like.
   *  Should not be called by user code directly.
   */
  virtual void outputMessage() const;

private:

  /// private undefined copy constructor
  DJDecompressIJG8Bit(const DJDecompressIJG8Bit&);

  /// private undefined copy assignment operator
  DJDecompressIJG8Bit& operator=(const DJDecompressIJG8Bit&);

  /// cleans up cinfo structure, called from destructor and error handlers
  void cleanup();

  /// codec parameters
  const DJCodecParameter *cparam;

  /// decompression structure
  jpeg_decompress_struct *cinfo;

  /// position of last suspend
  int suspension;

  /// temporary storage for row buffer during suspension
  void *jsampBuffer;

  /// Flag indicating if DICOM photometric interpretation is YCbCr
  OFBool dicomPhotometricInterpretationIsYCbCr;

  /// color model after decompression
  EP_Interpretation decompressedColorModel;

};

#endif

/*
 * CVS/RCS Log
 * $Log: djdijg8.h,v $
 * Revision 1.4  2005/12/08 16:59:24  meichel
 * Changed include path schema for all DCMTK header files
 *
 * Revision 1.3  2005/11/30 14:08:57  onken
 * Added check to decline automatic IJG color space conversion of signed pixel
 * data, because IJG lib only handles unsigned input for conversions.
 *
 * Revision 1.2  2001/11/19 15:13:28  meichel
 * Introduced verbose mode in module dcmjpeg. If enabled, warning
 *   messages from the IJG library are printed on ofConsole, otherwise
 *   the library remains quiet.
 *
 * Revision 1.1  2001/11/13 15:56:23  meichel
 * Initial release of module dcmjpeg
 *
 *
 */
