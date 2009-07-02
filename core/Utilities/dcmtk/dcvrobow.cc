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
 *  Purpose: Implementation of class DcmOtherByteOtherWord
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:41:57 $
 *  CVS/RCS Revision: $Revision: 1.45 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */
#include "dcmtk/ofstd/ofstd.h"
#include "dcmtk/ofstd/ofstream.h"
#include "dcmtk/dcmdata/dcvrobow.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmdata/dcswap.h"
#include "dcmtk/dcmdata/dcvm.h"

#define INCLUDE_CSTDIO
#define INCLUDE_CSTDLIB
#define INCLUDE_CSTRING
#include "dcmtk/ofstd/ofstdinc.h"


// ********************************


DcmOtherByteOtherWord::DcmOtherByteOtherWord(const DcmTag &tag,
                                             const Uint32 len)
  : DcmElement(tag, len)
{
}


DcmOtherByteOtherWord::DcmOtherByteOtherWord(const DcmOtherByteOtherWord &old)
  : DcmElement(old)
{
}


DcmOtherByteOtherWord::~DcmOtherByteOtherWord()
{
}


DcmOtherByteOtherWord &DcmOtherByteOtherWord::operator=(const DcmOtherByteOtherWord &obj)
{
    DcmElement::operator=(obj);
    return *this;
}


// ********************************


DcmEVR DcmOtherByteOtherWord::ident() const
{
    return Tag.getEVR();
}


unsigned long DcmOtherByteOtherWord::getVM()
{
    /* value multiplicity for OB/OW is defined as 1 */
    return 1;
}


OFCondition DcmOtherByteOtherWord::setVR(DcmEVR vr)
{
    Tag.setVR(vr);
    return EC_Normal;
}


// ********************************


void DcmOtherByteOtherWord::print(ostream &out,
                                  const size_t flags,
                                  const int level,
                                  const char * /*pixelFileName*/,
                                  size_t * /*pixelCounter*/)
{
    if (valueLoaded())
    {
        const DcmEVR evr = Tag.getEVR();
        Uint16 *wordValues = NULL;
        Uint8 *byteValues = NULL;
        /* get 8 or 16 bit data respectively */
        if (evr == EVR_OW || evr == EVR_lt)
            errorFlag = getUint16Array(wordValues);
        else
            errorFlag = getUint8Array(byteValues);
        /* check data */
        if ((wordValues != NULL) || (byteValues != NULL))
        {
            /* determine number of values to be printed */
            const unsigned int vrSize = (evr == EVR_OW || evr == EVR_lt) ? 4 : 2;
            const unsigned long count = (evr == EVR_OW || evr == EVR_lt) ? (Length / 2) : Length;
            unsigned long expectedLength = count * (vrSize + 1) - 1;
            const unsigned long printCount =
                ((expectedLength > DCM_OptPrintLineLength) && (flags & DCMTypes::PF_shortenLongTagValues)) ?
                (DCM_OptPrintLineLength - 3 /* for "..." */ + 1 /* for last "\" */) / (vrSize + 1) : count;
            unsigned long printedLength = printCount * (vrSize + 1) - 1;
            /* print line start with tag and VR */
            printInfoLineStart(out, flags, level);
            /* print multiple values */
            if (printCount > 0)
            {
                out << hex << setfill('0');
                if (evr == EVR_OW || evr == EVR_lt)
                {
                    /* print word values in hex mode */
                    out << setw(vrSize) << (*(wordValues++));
                    for (unsigned long i = 1; i < printCount; i++)
                        out << "\\" << setw(vrSize) << (*(wordValues++));
                } else {
                    /* print byte values in hex mode */
                    out << setw(vrSize) << OFstatic_cast(int, *(byteValues++));
                    for (unsigned long i = 1; i < printCount; i++)
                        out << "\\" << setw(vrSize) << OFstatic_cast(int, *(byteValues++));
                }
                /* reset i/o manipulators */
                out << dec << setfill(' ');
            }
            /* print trailing "..." if data has been truncated */
            if (printCount < count)
            {
                out << "...";
                printedLength += 3;
            }
            /* print line end with length, VM and tag name */
            printInfoLineEnd(out, flags, printedLength);
        } else
            printInfoLine(out, flags, level, "(no value available)");
    } else
        printInfoLine(out, flags, level, "(not loaded)");
}


void DcmOtherByteOtherWord::printPixel(ostream &out,
                                       const size_t flags,
                                       const int level,
                                       const char *pixelFileName,
                                       size_t *pixelCounter)
{
    if (pixelFileName != NULL)
    {
        /* create filename for pixel data file */
        OFString fname = pixelFileName;
        fname += ".";
        if (pixelCounter != NULL)
        {
            char num[20];
            sprintf(num, "%ld", OFstatic_cast(long, (*pixelCounter)++));
            fname += num;
        }
        fname += ".raw";
        /* create reference to pixel data file in dump output */
        OFString str = "=";
        str += fname;
        printInfoLine(out, flags, level, str.c_str());
        /* check whether pixel data file already exists */
        if (!OFStandard::fileExists(fname))
        {
            /* create binary file for pixel data */
            FILE *file = fopen(fname.c_str(), "wb");
            if (file != NULL)
            {
                if (Tag.getEVR() == EVR_OW || Tag.getEVR() == EVR_lt)
                {
                    /* write 16 bit data in little endian byte-order */
                    Uint16 *data = NULL;
                    getUint16Array(data);
                    if (data != NULL)
                    {
                        swapIfNecessary(EBO_LittleEndian, gLocalByteOrder, data, Length, sizeof(Uint16));
                        fwrite(data, sizeof(Uint16), OFstatic_cast(size_t, Length / sizeof(Uint16)), file);
                        swapIfNecessary(gLocalByteOrder, EBO_LittleEndian, data, Length, sizeof(Uint16));
                    }
                } else {
                    Uint8 *data = NULL;
                    getUint8Array(data);
                    if (data != NULL)
                        fwrite(data, sizeof(Uint8), OFstatic_cast(size_t, Length), file);
                }
                fclose(file);
            } else {
                ofConsole.lockCerr() << "Warning: can't open output file for pixel data: " << fname << endl;
                ofConsole.unlockCerr();
            }
        } else {
            ofConsole.lockCerr() << "Warning: output file for pixel data already exists: " << fname << endl;
            ofConsole.unlockCerr();
        }
    } else
        DcmOtherByteOtherWord::print(out, flags, level, pixelFileName, pixelCounter);
}


// ********************************


OFCondition DcmOtherByteOtherWord::alignValue()
{
    errorFlag = EC_Normal;
    /* add padding byte in case of 8 bit data */
    if ((Tag.getEVR() != EVR_OW && Tag.getEVR() != EVR_lt) && (Length > 0))
    {
        Uint8 *bytes = NULL;
        bytes = OFstatic_cast(Uint8 *, getValue(fByteOrder));
        /* check for odd length */
        if ((bytes != NULL) && ((Length & 1) != 0))
        {
            bytes[Length] = 0;
            Length++;
        }
    }
    return errorFlag;
}


void DcmOtherByteOtherWord::postLoadValue()
{
    if (dcmEnableAutomaticInputDataCorrection.get())
        alignValue();
}


// ********************************


OFCondition DcmOtherByteOtherWord::putUint8Array(const Uint8 *byteValue,
                                                 const unsigned long numBytes)
{
    errorFlag = EC_Normal;
    if (numBytes > 0)
    {
        /* check for valid 8 bit data */
        if ((byteValue != NULL) && (Tag.getEVR() != EVR_OW && Tag.getEVR() != EVR_lt))
        {
            errorFlag = putValue(byteValue, sizeof(Uint8) * OFstatic_cast(Uint32, numBytes));
            alignValue();
        } else
            errorFlag = EC_CorruptedData;
    } else
        putValue(NULL, 0);
    return errorFlag;
}


OFCondition DcmOtherByteOtherWord::putUint16Array(const Uint16 *wordValue,
                                                  const unsigned long numWords)
{
    errorFlag = EC_Normal;
    if (numWords > 0)
    {
        /* check for valid 16 bit data */
        if ((wordValue != NULL) && (Tag.getEVR() == EVR_OW || Tag.getEVR() == EVR_lt))
            errorFlag = putValue(wordValue, sizeof(Uint16) * OFstatic_cast(Uint32, numWords));
        else
            errorFlag = EC_CorruptedData;
    }
    else
        errorFlag = putValue(NULL, 0);
    return errorFlag;
}


// ********************************


OFCondition DcmOtherByteOtherWord::putString(const char *stringVal)
{
    errorFlag = EC_Normal;
    /* check input string */
    if ((stringVal != NULL) && (strlen(stringVal) > 0))
    {
        unsigned long vm = getVMFromString(stringVal);
        if (vm > 0)
        {
            const DcmEVR evr = Tag.getEVR();
            Uint8 *byteField = NULL;
            Uint16 *wordField = NULL;
            /* create new value field */
            if (evr == EVR_OW || evr == EVR_lt)
                wordField = new Uint16[vm];
            else
                byteField = new Uint8[vm];
            const char *s = stringVal;
            Uint16 intVal = 0;
            char *value;
            /* retrieve binary data from hexa-decimal string */
            for (unsigned long i = 0; (i < vm) && errorFlag.good(); i++)
            {
                /* get first value stored in 's', set 's' to beginning of the next value */
                value = getFirstValueFromString(s);
                if (value != NULL)
                {
                    /* integer overflow is currently not checked! */
                    if (sscanf(value, "%hx", &intVal) != 1)
                        errorFlag = EC_CorruptedData;
                    else if (evr == EVR_OW || evr == EVR_lt)
                        wordField[i] = OFstatic_cast(Uint16, intVal);
                    else
                        byteField[i] = OFstatic_cast(Uint8, intVal);
                    delete[] value;
                } else
                    errorFlag = EC_CorruptedData;
            }
            /* set binary data as the element value */
            if (errorFlag.good())
            {
                if (evr == EVR_OW || evr == EVR_lt)
                    errorFlag = putUint16Array(wordField, vm);
                else
                    errorFlag = putUint8Array(byteField, vm);
            }
            /* delete temporary buffers */
            delete[] byteField;
            delete[] wordField;
        } else
            putValue(NULL, 0);
    } else
        putValue(NULL, 0);
    return errorFlag;
}


// ********************************


OFCondition DcmOtherByteOtherWord::getUint8(Uint8 &byteVal,
                                            const unsigned long pos)
{
    /* get 8 bit data */
    Uint8 *uintValues = NULL;
    errorFlag = getUint8Array(uintValues);
    /* check data before returning */
    if (errorFlag.good())
    {
        if (uintValues == NULL)
            errorFlag = EC_IllegalCall;
        else if (pos >= getLength() /*bytes*/)
            errorFlag = EC_IllegalParameter;
        else
            byteVal = uintValues[pos];
    }
    /* clear value in case of error */
    if (errorFlag.bad())
	    byteVal = 0;
    return errorFlag;
}


OFCondition DcmOtherByteOtherWord::getUint8Array(Uint8 *&byteVals)
{
    errorFlag = EC_Normal;
    if (Tag.getEVR() != EVR_OW && Tag.getEVR() != EVR_lt)
        byteVals = OFstatic_cast(Uint8 *, getValue());
    else
        errorFlag = EC_IllegalCall;
    return errorFlag;
}


// ********************************


OFCondition DcmOtherByteOtherWord::getUint16(Uint16 &wordVal,
                                             const unsigned long pos)
{
    Uint16 *uintValues = NULL;
    errorFlag = getUint16Array(uintValues);
    /* check data before returning */
    if (errorFlag.good())
    {
        if (uintValues == NULL)
            errorFlag = EC_IllegalCall;
        else if (pos >= getLength() / sizeof(Uint16) /*words*/)
            errorFlag = EC_IllegalParameter;
        else
            wordVal = uintValues[pos];
    }
    /* clear value in case of error */
    if (errorFlag.bad())
	    wordVal = 0;
    return errorFlag;
}


OFCondition DcmOtherByteOtherWord::getUint16Array(Uint16 *&wordVals)
{
    errorFlag = EC_Normal;
    if (Tag.getEVR() == EVR_OW || Tag.getEVR() == EVR_lt)
        wordVals = OFstatic_cast(Uint16 *, getValue());
    else
        errorFlag = EC_IllegalCall;
    return errorFlag;
}


// ********************************


OFCondition DcmOtherByteOtherWord::getOFString(OFString &stringVal,
                                               const unsigned long pos,
                                               OFBool /*normalize*/)
{
    if (Tag.getEVR() == EVR_OW || Tag.getEVR() == EVR_lt)
    {
        Uint16 uint16Val;
        /* get the specified numeric value (16 bit) */
        errorFlag = getUint16(uint16Val, pos);
        if (errorFlag.good())
        {
            /* ... and convert it to a character string (hex mode) */
            char buffer[32];
            sprintf(buffer, "%4.4hx", uint16Val);
            /* assign result */
            stringVal = buffer;
        }
    } else {
        Uint8 uint8Val;
        /* get the specified numeric value (8 bit) */
        errorFlag = getUint8(uint8Val, pos);
        if (errorFlag.good())
        {
            /* ... and convert it to a character string (hex mode) */
            char buffer[32];
            sprintf(buffer, "%2.2hx", uint8Val);
            /* assign result */
            stringVal = buffer;
        }
    }
    return errorFlag;
}


OFCondition DcmOtherByteOtherWord::getOFStringArray(OFString &stringVal,
                                                    OFBool /*normalize*/)
{
    if (Tag.getEVR() == EVR_OW || Tag.getEVR() == EVR_lt)
    {
        /* get array of 16 bit values */
        Uint16 *uint16Vals = OFstatic_cast(Uint16 *, getValue());
        const size_t count = OFstatic_cast(size_t, getLength() / sizeof(Uint16));
        if ((uint16Vals != NULL) && (count > 0))
        {
            /* reserve number of bytes expected */
            stringVal.reserve(5 * count);
            char *bufPtr = OFconst_cast(char *, stringVal.c_str());
            /* for all array elements ... */
            for (size_t i = 0; i < count; i++)
            {
                /* ... convert numeric value to hexadecimal string representation */
                sprintf(bufPtr, "%4.4hx\\", uint16Vals[i]);
                bufPtr += 5;
            }
            /* remove last '\' */
            *(--bufPtr) = '\0';
            errorFlag = EC_Normal;
        } else
            errorFlag = EC_IllegalCall;
    } else {
        /* get array of 8 bit values */
        Uint8 *uint8Vals = OFstatic_cast(Uint8 *, getValue());
        const size_t count = OFstatic_cast(size_t, getLength());
        if ((uint8Vals != NULL) && (count > 0))
        {
            /* reserve number of bytes expected */
            stringVal.reserve(3 * count);
            char *bufPtr = OFconst_cast(char *, stringVal.c_str());
            /* for all array elements ... */
            for (size_t i = 0; i < count; i++)
            {
                /* ... convert numeric value to hexadecimal string representation */
                sprintf(bufPtr, "%2.2hx\\", uint8Vals[i]);
                bufPtr += 3;
            }
            /* remove last '\' */
            *(--bufPtr) = '\0';
            errorFlag = EC_Normal;
        } else
            errorFlag = EC_IllegalCall;
    }
    return errorFlag;
}


// ********************************


OFCondition DcmOtherByteOtherWord::verify(const OFBool autocorrect)
{
    errorFlag = EC_Normal;
    if (autocorrect)
        errorFlag = alignValue();
    return errorFlag;
}


// ********************************


OFBool DcmOtherByteOtherWord::canWriteXfer(const E_TransferSyntax newXfer,
                                           const E_TransferSyntax /*oldXfer*/)
{
    DcmXfer newXferSyn(newXfer);
    return (Tag != DCM_PixelData) || !newXferSyn.isEncapsulated();
}


// ********************************


OFCondition DcmOtherByteOtherWord::write(DcmOutputStream &outStream,
                                         const E_TransferSyntax oxfer,
                                         const E_EncodingType enctype)
{
    if (fTransferState == ERW_notInitialized)
        errorFlag = EC_IllegalCall;
    else
    {
        if (fTransferState == ERW_init) alignValue();
        /* call inherited method */
        errorFlag = DcmElement::write(outStream, oxfer, enctype);
    }
    return errorFlag;
}


OFCondition DcmOtherByteOtherWord::writeSignatureFormat(DcmOutputStream &outStream,
                                                        const E_TransferSyntax oxfer,
                                                        const E_EncodingType enctype)
{
    if (fTransferState == ERW_notInitialized)
        errorFlag = EC_IllegalCall;
    else
    {
        if (fTransferState == ERW_init) alignValue();
        /* call inherited method */
        errorFlag = DcmElement::writeSignatureFormat(outStream, oxfer, enctype);
    }
    return errorFlag;
}


// ********************************


OFCondition DcmOtherByteOtherWord::writeXML(ostream &out,
                                            const size_t flags)
{
    /* XML start tag: <element tag="gggg,eeee" vr="XX" ...> */
    if (!(flags & DCMTypes::XF_writeBinaryData))
        writeXMLStartTag(out, flags, "binary=\"hidden\"");
    else if (flags & DCMTypes::XF_encodeBase64)
        writeXMLStartTag(out, flags, "binary=\"base64\"");
    else
        writeXMLStartTag(out, flags, "binary=\"yes\"");
    /* write element value (if loaded) */
    if (valueLoaded() && (flags & DCMTypes::XF_writeBinaryData))
    {
        OFString value;
        /* encode binary data as Base64 */
        if (flags & DCMTypes::XF_encodeBase64)
        {
            Uint8 *byteValues = OFstatic_cast(Uint8 *, getValue());
            if (Tag.getEVR() == EVR_OW || Tag.getEVR() == EVR_lt)
            {
                /* Base64 encoder requires big endian input data */
                swapIfNecessary(gLocalByteOrder, EBO_BigEndian, byteValues, Length, sizeof(Uint16));
            }
            out << OFStandard::encodeBase64(byteValues, OFstatic_cast(size_t, Length), value);
        } else {
            /* encode as sequence of hexadecimal numbers */
            if (getOFStringArray(value).good())
                out << value;
        }
    }
    /* XML end tag: </element> */
    writeXMLEndTag(out, flags);
    /* always report success */
    return EC_Normal;
}


/*
** CVS/RCS Log:
** $Log: dcvrobow.cc,v $
** Revision 1.45  2005/12/08 15:41:57  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.44  2005/11/15 16:59:25  meichel
** Added new pseudo VR type EVR_lt that is used for LUT Data when read in
**   implicit VR, which may be US, SS or OW. DCMTK always treats EVR_lt like OW.
**
** Revision 1.43  2004/02/04 16:49:49  joergr
** Adapted type casts to new-style typecast operators defined in ofcast.h.
** Removed acknowledgements with e-mail addresses from CVS log.
**
** Revision 1.42  2003/04/17 15:59:45  joergr
** Use method OFString::c_str() instead of OFString::operator[] to avoid range
** checking (which implies an "expensive" strlen() call).
**
** Revision 1.41  2002/12/06 13:12:36  joergr
** Enhanced "print()" function by re-working the implementation and replacing
** the boolean "showFullData" parameter by a more general integer flag.
** Made source code formatting more consistent with other modules/files.
**
** Revision 1.40  2002/11/27 12:06:57  meichel
** Adapted module dcmdata to use of new header file ofstdinc.h
**
** Revision 1.39  2002/08/27 16:55:59  meichel
** Initial release of new DICOM I/O stream classes that add support for stream
**   compression (deflated little endian explicit VR transfer syntax)
**
** Revision 1.38  2002/07/08 14:44:42  meichel
** Improved dcmdata behaviour when reading odd tag length. Depending on the
**   global boolean flag dcmAcceptOddAttributeLength, the parser now either accepts
**   odd length attributes or implements the old behaviour, i.e. assumes a real
**   length larger by one.
**
** Revision 1.37  2002/05/14 08:22:56  joergr
** Added support for Base64 (MIME) encoded binary data.
**
** Revision 1.36  2002/04/25 10:32:16  joergr
** Added getOFString() implementation.
** Added/modified getOFStringArray() implementation.
** Added support for XML output of DICOM objects.
**
** Revision 1.35  2002/04/16 13:43:25  joergr
** Added configurable support for C++ ANSI standard includes (e.g. streams).
**
** Revision 1.34  2001/10/02 11:48:33  joergr
** Added getUint8/16 routines to class DcmOtherByteOtherWord.
**
** Revision 1.33  2001/09/25 17:19:58  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.32  2001/06/01 15:49:18  meichel
** Updated copyright header
**
** Revision 1.31  2000/11/07 16:56:24  meichel
** Initial release of dcmsign module for DICOM Digital Signatures
**
** Revision 1.30  2000/04/14 15:55:09  meichel
** Dcmdata library code now consistently uses ofConsole for error output.
**
** Revision 1.29  2000/03/08 16:26:49  meichel
** Updated copyright header.
**
** Revision 1.28  2000/03/07 15:41:02  joergr
** Added explicit type casts to make Sun CC 2.0.1 happy.
**
** Revision 1.27  2000/03/06 16:08:05  meichel
** Changed a couple of definitions that implied that Uint32 or size_t are long
**
** Revision 1.26  2000/03/03 14:05:39  meichel
** Implemented library support for redirecting error messages into memory
**   instead of printing them to stdout/stderr for GUI applications.
**
** Revision 1.25  2000/02/23 15:12:07  meichel
** Corrected macro for Borland C++ Builder 4 workaround.
**
** Revision 1.24  2000/02/10 16:04:59  joergr
** Enhanced handling of PixelData/Item element. Externally stored raw data is
** now always imported as little endian and swapped if necessary. This change
** reflects the new 'export' feature of dcmdump.
**
** Revision 1.23  2000/02/10 10:52:24  joergr
** Added new feature to dcmdump (enhanced print method of dcmdata): write
** pixel data/item value fields to raw files.
**
** Revision 1.22  2000/02/03 16:55:20  joergr
** Avoid EVR_pixelItem in comparisons (compare with != EVR_OW instead).
**
** Revision 1.21  2000/02/03 16:35:58  joergr
** Corrected bug that caused wrong calculation of group length for sequence
** of items (e.g. encapsulated pixel data).
**
** Revision 1.20  2000/02/01 10:12:11  meichel
** Avoiding to include <stdlib.h> as extern "C" on Borland C++ Builder 4,
**   workaround for bug in compiler header files.
**
** Revision 1.19  1999/03/31 09:25:55  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.18  1997/07/21 08:11:43  andreas
** - Support for CP 14. PixelData and OverlayData can have VR OW or OB
**   (depending on the transfer syntax). New internal value
**   representation (only for ident()) for OverlayData.
** - New environment for encapsulated pixel representations. DcmPixelData
**   can contain different representations and uses codecs to convert
**   between them. Codecs are derived from the DcmCodec class. New error
**   codes are introduced for handling of representations. New internal
**   value representation (only for ident()) for PixelData
** - Replace all boolean types (BOOLEAN, CTNBOOLEAN, DICOM_BOOL, BOOL)
**   with one unique boolean type OFBool.
**
** Revision 1.17  1997/07/07 07:51:35  andreas
** - Changed type for Tag attribute in DcmObject from prointer to value
** - Enhanced (faster) byte swapping routine. swapIfNecessary moved from
**   a method in DcmObject to a general function.
**
** Revision 1.16  1997/07/03 15:10:15  andreas
** - removed debugging functions Bdebug() and Edebug() since
**   they write a static array and are not very useful at all.
**   Cdebug and Vdebug are merged since they have the same semantics.
**   The debugging functions in dcmdata changed their interfaces
**   (see dcmdata/include/dcdebug.h)
**
** Revision 1.15  1997/06/13 13:07:31  andreas
** - Corrected printing of OW values. The length of the output array was
**   computed incorrectly.
**
** Revision 1.14  1997/05/27 13:49:03  andreas
** - Add method canWriteXfer to class DcmObject and all derived classes.
**   This method checks whether it is possible to convert the original
**   transfer syntax to an new transfer syntax. The check is used in the
**   dcmconv utility to prohibit the change of a compressed transfer
**   syntax to a uncompressed.
**
** Revision 1.13  1997/05/22 16:54:20  andreas
** - Corrected wrong output length in print routine
**
** Revision 1.12  1997/05/16 08:31:28  andreas
** - Revised handling of GroupLength elements and support of
**   DataSetTrailingPadding elements. The enumeratio E_GrpLenEncoding
**   got additional enumeration values (for a description see dctypes.h).
**   addGroupLength and removeGroupLength methods are replaced by
**   computeGroupLengthAndPadding. To support Padding, the parameters of
**   element and sequence write functions changed.
**
** Revision 1.11  1997/04/18 08:17:20  andreas
** - The put/get-methods for all VRs did not conform to the C++-Standard
**   draft. Some Compilers (e.g. SUN-C++ Compiler, Metroworks
**   CodeWarrier, etc.) create many warnings concerning the hiding of
**   overloaded get methods in all derived classes of DcmElement.
**   So the interface of all value representation classes in the
**   library are changed rapidly, e.g.
**   OFCondition get(Uint16 & value, const unsigned long pos);
**   becomes
**   OFCondition getUint16(Uint16 & value, const unsigned long pos);
**   All (retired) "returntype get(...)" methods are deleted.
**   For more information see dcmdata/include/dcelem.h
**
** Revision 1.10  1997/03/26 17:15:59  hewett
** Added very preliminary support for Unknown VR (UN) described in
** Supplement 14.  WARNING: handling of unknown attributes with undefined
** length is not yet supported.
**
** Revision 1.9  1996/08/05 08:46:20  andreas
** new print routine with additional parameters:
**         - print into files
**         - fix output length for elements
** corrected error in search routine with parameter ESM_fromStackTop
**
** Revision 1.8  1996/05/06 10:22:47  meichel
** Copy constructor now handles case when VR=unknown.
**
** Revision 1.7  1996/04/16 16:05:24  andreas
** - better support und bug fixes for NULL element value
**
** Revision 1.6  1996/03/26 09:59:36  meichel
** corrected bug (deletion of const char *) which prevented compilation on NeXT
**
** Revision 1.5  1996/01/29 13:38:33  andreas
** - new put method for every VR to put value as a string
** - better and unique print methods
**
** Revision 1.4  1996/01/09 11:06:50  andreas
** New Support for Visual C++
** Correct problems with inconsistent const declarations
** Correct error in reading Item Delimitation Elements
**
** Revision 1.3  1996/01/05 13:27:51  andreas
** - changed to support new streaming facilities
** - unique read/write methods for file and block transfer
** - more cleanups
**
** Revision 1.2  1995/11/23 17:03:07  hewett
** Updated for loadable data dictionary.  Some cleanup (more to do).
**
*/
