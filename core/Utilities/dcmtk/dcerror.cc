/*
 *
 *  Copyright (C) 1994-2005, OFFIS
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
 *  Module:  dcmdata
 *
 *  Author:  Andrew Hewett
 *
 *  Purpose: Error handling, codes and strings
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:41:09 $
 *  Source File:      $Source: /share/dicom/cvs-depot/dcmtk/dcmdata/libsrc/dcerror.cc,v $
 *  CVS/RCS Revision: $Revision: 1.16 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */

#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */
#include "dcmtk/dcmdata/dcerror.h"

const OFConditionConst ECC_InvalidTag(                 OFM_dcmdata,  1, OF_error, "Invalid Tag"                            );
const OFConditionConst ECC_TagNotFound(                OFM_dcmdata,  2, OF_error, "Tag Not Found"                          );
const OFConditionConst ECC_InvalidVR(                  OFM_dcmdata,  3, OF_error, "Invalid VR"                             );
const OFConditionConst ECC_InvalidStream(              OFM_dcmdata,  4, OF_error, "Invalid Stream"                         );
const OFConditionConst ECC_EndOfStream(                OFM_dcmdata,  5, OF_error, "End Of Stream"                          );
const OFConditionConst ECC_CorruptedData(              OFM_dcmdata,  6, OF_error, "Corrupted Data"                         );
const OFConditionConst ECC_IllegalCall(                OFM_dcmdata,  7, OF_error, "Illegal Call, perhaps wrong parameters" );
const OFConditionConst ECC_SequEnd(                    OFM_dcmdata,  8, OF_error, "Sequence End"                           );
const OFConditionConst ECC_DoubledTag(                 OFM_dcmdata,  9, OF_error, "Doubled Tag"                            );
const OFConditionConst ECC_StreamNotifyClient(         OFM_dcmdata, 10, OF_error, "I/O suspension or premature end of stream" );
const OFConditionConst ECC_WrongStreamMode(            OFM_dcmdata, 11, OF_error, "Mode (R/W, random/sequence) is wrong"   );
const OFConditionConst ECC_ItemEnd(                    OFM_dcmdata, 12, OF_error, "Item End"                               );
const OFConditionConst ECC_RepresentationNotFound(     OFM_dcmdata, 13, OF_error, "Pixel representation not found"         );
const OFConditionConst ECC_CannotChangeRepresentation( OFM_dcmdata, 14, OF_error, "Pixel representation cannot be changed" );
const OFConditionConst ECC_UnsupportedEncoding(        OFM_dcmdata, 15, OF_error, "Unsupported compression or encryption"  );
// error code 16 is reserved for zlib-related error messages
const OFConditionConst ECC_PutbackFailed(              OFM_dcmdata, 17, OF_error, "Parser failure: Putback operation failed" );
// error code 18 is reserved for file read error messages
// error code 19 is reserved for file write error messages
const OFConditionConst ECC_DoubleCompressionFilters(   OFM_dcmdata, 20, OF_error, "Too many compression filters"           );
const OFConditionConst ECC_ApplicationProfileViolated( OFM_dcmdata, 21, OF_error, "Storage media application profile violated" );
// error code 22 is reserved for dcmodify error messages

const OFCondition EC_InvalidTag(                 ECC_InvalidTag);
const OFCondition EC_TagNotFound(                ECC_TagNotFound);
const OFCondition EC_InvalidVR(                  ECC_InvalidVR);
const OFCondition EC_InvalidStream(              ECC_InvalidStream);
const OFCondition EC_EndOfStream(                ECC_EndOfStream);
const OFCondition EC_CorruptedData(              ECC_CorruptedData);
const OFCondition EC_IllegalCall(                ECC_IllegalCall);
const OFCondition EC_SequEnd(                    ECC_SequEnd);
const OFCondition EC_DoubledTag(                 ECC_DoubledTag);
const OFCondition EC_StreamNotifyClient(         ECC_StreamNotifyClient);
const OFCondition EC_WrongStreamMode(            ECC_WrongStreamMode);
const OFCondition EC_ItemEnd(                    ECC_ItemEnd);
const OFCondition EC_RepresentationNotFound(     ECC_RepresentationNotFound);
const OFCondition EC_CannotChangeRepresentation( ECC_CannotChangeRepresentation);
const OFCondition EC_UnsupportedEncoding(        ECC_UnsupportedEncoding);
const OFCondition EC_PutbackFailed(              ECC_PutbackFailed);
const OFCondition EC_DoubleCompressionFilters(   ECC_DoubleCompressionFilters);
const OFCondition EC_ApplicationProfileViolated( ECC_ApplicationProfileViolated);

const char *dcmErrorConditionToString(OFCondition cond)
{
  return cond.text();
}


/*
** CVS/RCS Log:
** $Log: dcerror.cc,v $
** Revision 1.16  2005/12/08 15:41:09  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.15  2004/11/05 17:20:31  onken
** Added reservation for dcmodify error messages.
**
** Revision 1.14  2002/12/06 12:18:57  joergr
** Added new error status "EC_ApplicationProfileViolated".
**
** Revision 1.13  2002/08/27 16:55:47  meichel
** Initial release of new DICOM I/O stream classes that add support for stream
**   compression (deflated little endian explicit VR transfer syntax)
**
** Revision 1.12  2001/09/25 17:19:50  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.11  2001/06/01 15:49:04  meichel
** Updated copyright header
**
** Revision 1.10  2000/03/08 16:26:35  meichel
** Updated copyright header.
**
** Revision 1.9  2000/02/23 15:11:52  meichel
** Corrected macro for Borland C++ Builder 4 workaround.
**
** Revision 1.8  2000/02/01 10:12:07  meichel
** Avoiding to include <stdlib.h> as extern "C" on Borland C++ Builder 4,
**   workaround for bug in compiler header files.
**
** Revision 1.7  1999/03/31 09:25:27  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.6  1997/10/01 08:44:12  meichel
** Including <unistd.h> if available.
**
** Revision 1.5  1997/07/21 08:17:41  andreas
** - New environment for encapsulated pixel representations. DcmPixelData
**   can contain different representations and uses codecs to convert
**   between them. Codecs are derived from the DcmCodec class. New error
**   codes are introduced for handling of representations. New internal
**   value representation (only for ident()) for PixelData
**
** Revision 1.4  1997/05/22 16:55:05  andreas
** - Added new error code EC_NotImplemented
**
** Revision 1.3  1996/01/29 13:38:26  andreas
** - new put method for every VR to put value as a string
** - better and unique print methods
**
** Revision 1.2  1996/01/05 13:27:36  andreas
** - changed to support new streaming facilities
** - unique read/write methods for file and block transfer
** - more cleanups
**
** Revision 1.1  1995/11/23 17:02:44  hewett
** Updated for loadable data dictionary.  Some cleanup (more to do).
**
*/

