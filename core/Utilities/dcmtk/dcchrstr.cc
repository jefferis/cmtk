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
 *  Author:  Andreas Barth
 *
 *  Purpose: Implementation of class DcmCharString
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:40:57 $
 *  Source File:      $Source: /share/dicom/cvs-depot/dcmtk/dcmdata/libsrc/dcchrstr.cc,v $
 *  CVS/RCS Revision: $Revision: 1.9 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */

//
// This implementation does not support 16 bit character sets. Since 8 bit
// character sets are supported by the class DcmByteString the class
// DcmCharString is derived from DcmByteString without any extensions.
// No special implementation is necessary.
//
// If the extension for 16 bit character sets will be implemented this class
// must be derived directly from DcmElement. This class is designed to support
// the value representations (LO, LT, PN, SH, ST, UT). They are a problem because
// their value width (1, 2, .. Bytes) is specified by the element
// SpecificCharacterSet (0008, 0005) and an implementation must support
// different value widths that cannot be derived from the value representation.
//


#include "dcmtk/dcmdata/dcchrstr.h"


DcmCharString::DcmCharString(const DcmTag &tag, const Uint32 len)
  : DcmByteString(tag, len)
{
}

DcmCharString::DcmCharString(const DcmCharString &old)
  : DcmByteString(old)
{
}

DcmCharString::~DcmCharString(void)
{
}


DcmCharString &DcmCharString::operator=(const DcmCharString &obj)
{
    DcmByteString::operator=(obj);
    return *this;
}


/*
 * CVS/RCS Log:
 * $Log: dcchrstr.cc,v $
 * Revision 1.9  2005/12/08 15:40:57  meichel
 * Changed include path schema for all DCMTK header files
 *
 * Revision 1.8  2002/12/06 13:08:18  joergr
 * Made source code formatting more consistent with other modules/files.
 *
 * Revision 1.7  2001/06/01 15:48:59  meichel
 * Updated copyright header
 *
 * Revision 1.6  2000/03/08 16:26:30  meichel
 * Updated copyright header.
 *
 * Revision 1.5  1999/03/31 09:25:17  meichel
 * Updated copyright header in module dcmdata
 *
 *
 */
