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
 *  Purpose: Implementation of class DcmDecimalString
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:41:50 $
 *  CVS/RCS Revision: $Revision: 1.18 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */

#include "dcmtk/dcmdata/dcvrds.h"
#include "dcmtk/ofstd/ofstring.h"
#include "dcmtk/ofstd/ofstd.h"


// ********************************


DcmDecimalString::DcmDecimalString(const DcmTag &tag,
                                   const Uint32 len)
  : DcmByteString(tag, len)
{
    maxLength = 16;
}


DcmDecimalString::DcmDecimalString(const DcmDecimalString &old)
  : DcmByteString(old)
{
}


DcmDecimalString::~DcmDecimalString()
{
}


DcmDecimalString &DcmDecimalString::operator=(const DcmDecimalString &obj)
{
    DcmByteString::operator=(obj);
    return *this;
}


// ********************************


DcmEVR DcmDecimalString::ident() const
{
    return EVR_DS;
}


// ********************************


OFCondition DcmDecimalString::getFloat64(Float64 &doubleVal,
                                         const unsigned long pos)
{
    /* get decimal string value */
    OFString str;
    OFCondition l_error = getOFString(str, pos, OFTrue /*normalize*/);
    if (l_error.good())
    {
        OFBool success = OFFalse;
        /* convert string to float value */
        doubleVal = OFStandard::atof(str.c_str(), &success);
        if (!success)
            l_error = EC_CorruptedData;
    }
    return l_error;
}


// ********************************


OFCondition DcmDecimalString::getOFString(OFString &stringVal,
                                          const unsigned long pos,
                                          OFBool normalize)
{
    /* call inherited method */
    OFCondition l_error = DcmByteString::getOFString(stringVal, pos);
    /* normalize string if required */
    if (l_error.good() && normalize)
        normalizeString(stringVal, !MULTIPART, DELETE_LEADING, DELETE_TRAILING);
    return l_error;
}


/*
** CVS/RCS Log:
** $Log: dcvrds.cc,v $
** Revision 1.18  2005/12/08 15:41:50  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.17  2004/01/16 13:49:00  joergr
** Removed acknowledgements with e-mail addresses from CVS log.
**
** Revision 1.16  2002/12/06 13:20:49  joergr
** Enhanced "print()" function by re-working the implementation and replacing
** the boolean "showFullData" parameter by a more general integer flag.
** Made source code formatting more consistent with other modules/files.
**
** Revision 1.15  2002/06/20 12:06:16  meichel
** Changed toolkit to use OFStandard::atof instead of atof, strtod or
**   sscanf for all string to double conversions that are supposed to
**   be locale independent
**
** Revision 1.14  2002/04/25 10:28:42  joergr
** Removed getOFStringArray() implementation.
**
** Revision 1.13  2001/09/25 17:19:56  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.12  2001/06/01 15:49:15  meichel
** Updated copyright header
**
** Revision 1.11  2000/03/08 16:26:46  meichel
** Updated copyright header.
**
** Revision 1.10  1999/03/31 09:25:49  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.9  1998/11/12 16:48:23  meichel
** Implemented operator= for all classes derived from DcmObject.
**
** Revision 1.8  1997/08/29 13:11:45  andreas
** Corrected Bug in getOFStringArray Implementation
**
** Revision 1.7  1997/08/29 08:32:57  andreas
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
** Revision 1.6  1997/07/03 15:10:10  andreas
** - removed debugging functions Bdebug() and Edebug() since
**   they write a static array and are not very useful at all.
**   Cdebug and Vdebug are merged since they have the same semantics.
**   The debugging functions in dcmdata changed their interfaces
**   (see dcmdata/include/dcdebug.h)
**
** Revision 1.5  1997/07/03 10:22:45  andreas
** - corrected Bug in DcmDecimalString::getFloat64
**
** Revision 1.4  1997/05/12 07:38:26  andreas
** - new get-Methods for DcmDecimalString: getFloat64 and
**   DcmIntegerString: getSint32
**
** Revision 1.3  1996/01/05 13:27:47  andreas
** - changed to support new streaming facilities
** - unique read/write methods for file and block transfer
** - more cleanups
**
*/
