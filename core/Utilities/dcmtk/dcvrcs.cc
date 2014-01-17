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
 *  Author:  Gerd Ehlers, Andreas Barth
 *
 *  Purpose: class DcmCodeString
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:41:48 $
 *  Source File:      $Source: /share/dicom/cvs-depot/dcmtk/dcmdata/libsrc/dcvrcs.cc,v $
 *  CVS/RCS Revision: $Revision: 1.15 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */

#define INCLUDE_CCTYPE
#include "dcmtk/ofstd/ofstdinc.h"

#include "dcmtk/dcmdata/dcvrcs.h"


#define MAX_CS_LENGTH 16


// ********************************


DcmCodeString::DcmCodeString(const DcmTag &tag,
                             const Uint32 len)
  : DcmByteString(tag, len)
{
    maxLength = MAX_CS_LENGTH;
}


DcmCodeString::DcmCodeString(const DcmCodeString &old)
  : DcmByteString(old)
{
}


DcmCodeString::~DcmCodeString()
{
}


DcmCodeString &DcmCodeString::operator=(const DcmCodeString &obj)
{
    DcmByteString::operator=(obj);
    return *this;
}


// ********************************


DcmEVR DcmCodeString::ident() const
{
    return EVR_CS;
}


// ********************************


OFCondition DcmCodeString::getOFString(OFString &stringVal,
                                       const unsigned long pos,
                                       OFBool normalize)
{
    OFCondition l_error = DcmByteString::getOFString(stringVal, pos, normalize);
    if (l_error.good() && normalize)
        normalizeString(stringVal, !MULTIPART, DELETE_LEADING, DELETE_TRAILING);
    return l_error;
}


// ********************************


OFBool DcmCodeString::checkVR(const OFString &value,
                              size_t *pos,
                              const OFBool checkLength)
{
    char c;
    size_t i;
    const size_t length = value.length();
    const size_t maxlen = (length < MAX_CS_LENGTH) || (!checkLength) ? length : MAX_CS_LENGTH;
    /* iterate over all characters (up to the maximum) */
    for (i = 0; i < maxlen; i++)
    {
        c = value.at(i);
        /* check for valid CS character: A-Z, 0-9, _ and ' ' (space) */
        if ((c != ' ') && (c != '_') && !isdigit(c) && !(isalpha(c) && isupper(c)))
            break;
    }
    /* return position of first invalid character (eos if all valid) */
    if (pos != NULL)
        *pos = i;
    /* OFFalse in case of any invalid character */
    return (i == length);
}


/*
** CVS/RCS Log:
** $Log: dcvrcs.cc,v $
** Revision 1.15  2005/12/08 15:41:48  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.14  2003/06/12 15:07:13  joergr
** Added static function checkVR().
**
** Revision 1.13  2002/12/06 13:20:49  joergr
** Enhanced "print()" function by re-working the implementation and replacing
** the boolean "showFullData" parameter by a more general integer flag.
** Made source code formatting more consistent with other modules/files.
**
** Revision 1.12  2002/04/25 10:28:07  joergr
** Removed getOFStringArray() implementation.
**
** Revision 1.11  2001/09/25 17:19:55  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.10  2001/06/01 15:49:15  meichel
** Updated copyright header
**
** Revision 1.9  2000/03/08 16:26:45  meichel
** Updated copyright header.
**
** Revision 1.8  1999/03/31 09:25:48  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.7  1998/11/12 16:48:22  meichel
** Implemented operator= for all classes derived from DcmObject.
**
** Revision 1.6  1997/08/29 13:11:44  andreas
** Corrected Bug in getOFStringArray Implementation
**
** Revision 1.5  1997/08/29 08:32:56  andreas
** - Added methods getOFString and getOFStringArray for all
**   string VRs. These methods are able to normalise the value, i. e.
**   to remove leading and trailing spaces. This will be done only if
**   it is described in the standard that these spaces are not relevant.
**   These methods do not test the strings for conformance, this means
**   especially that they do not delete spaces where they are not allowed!
**   getOFStringArray returns the string with all its parts separated by \
**   and getOFString returns only one value of the string.
**   CAUTION: Currently getString returns a string with trailing
**   spaces removed (if dcmEnableAutomaticInputDataCorrection == OFTrue) and
**   truncates the original string (since it is not copied!). If you rely on this
**   behaviour please change your application now.
**   Future changes will ensure that getString returns the original
**   string from the DICOM object (NULL terminated) inclusive padding.
**   Currently, if you call getOF... before calling getString without
**   normalisation, you can get the original string read from the DICOM object.
**
** Revision 1.4  1997/07/03 15:10:09  andreas
** - removed debugging functions Bdebug() and Edebug() since
**   they write a static array and are not very useful at all.
**   Cdebug and Vdebug are merged since they have the same semantics.
**   The debugging functions in dcmdata changed their interfaces
**   (see dcmdata/include/dcdebug.h)
**
** Revision 1.3  1996/01/05 13:27:46  andreas
** - changed to support new streaming facilities
** - unique read/write methods for file and block transfer
** - more cleanups
**
*/
