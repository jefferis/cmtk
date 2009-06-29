/*
 *
 *  Copyright (C) 1996-2005, OFFIS
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
 *  Module:  dcmimgle
 *
 *  Author:  Joerg Riesmeier
 *
 *  Purpose: DicomDocument (Source)
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:42:48 $
 *  CVS/RCS Revision: $Revision: 1.16 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/ofstd/ofstring.h"

#include "dcmtk/dcmimgle/didocu.h"
#include "dcmtk/dcmimgle/diutils.h"


/*----------------*
 *  constructors  *
 *----------------*/

DiDocument::DiDocument(const char *filename,
                       const unsigned long flags,
                       const unsigned long fstart,
                       const unsigned long fcount)
  : Object(NULL),
    FileFormat(new DcmFileFormat()),
    Xfer(EXS_Unknown),
    FrameStart(fstart),
    FrameCount(fcount),
    Flags(flags)
{
    if (FileFormat)
    {

        if (FileFormat->loadFile(filename).bad())
        {
            if (DicomImageClass::checkDebugLevel(DicomImageClass::DL_Errors))
            {
                ofConsole.lockCerr() << "ERROR: can't read file '" << filename << "' !" << endl;
                ofConsole.unlockCerr();
            }
            delete FileFormat;
            FileFormat = NULL;
        } else {
            Object = FileFormat->getDataset();
            if (Object != NULL)
            {
                Xfer = OFstatic_cast(DcmDataset *, Object)->getOriginalXfer();
                convertPixelData();
            }
        }
    }
}


DiDocument::DiDocument(DcmObject *object,
                       const E_TransferSyntax xfer,
                       const unsigned long flags,
                       const unsigned long fstart,
                       const unsigned long fcount)
  : Object(NULL),
    FileFormat(NULL),
    Xfer(xfer),
    FrameStart(fstart),
    FrameCount(fcount),
    Flags(flags)
{
    if (object != NULL)
    {
        if (object->ident() == EVR_fileFormat)
        {
            /* store reference to DICOM file format to be deleted on object destruction */
            if (Flags & CIF_TakeOverExternalDataset)
                FileFormat = OFstatic_cast(DcmFileFormat *, object);
            Object = OFstatic_cast(DcmFileFormat *, object)->getDataset();
        } else
            Object = object;
        if (Object != NULL)
        {
            if (Xfer == EXS_Unknown)
                Xfer = OFstatic_cast(DcmDataset *, Object)->getOriginalXfer();
            convertPixelData();
        }
    }
}




void DiDocument::convertPixelData()
{
    DcmStack pstack;
    // convert pixel data to uncompressed format if required
    if (search(DCM_PixelData, pstack))
    {
        DcmPixelData *pixel = OFstatic_cast(DcmPixelData *, pstack.top());
        pstack.clear();
        // push reference to DICOM dataset on the stack (required for decompression process)
        pstack.push(Object);
        pstack.push(pixel);                         // dummy stack entry
        if ((pixel != NULL) && pixel->chooseRepresentation(EXS_LittleEndianExplicit, NULL, pstack).good())
        {
            // set transfer syntax to unencapsulated/uncompressed
            if (DcmXfer(Xfer).isEncapsulated())
                Xfer = EXS_LittleEndianExplicit;
        } else {
            if (DicomImageClass::checkDebugLevel(DicomImageClass::DL_Errors))
            {
                ofConsole.lockCerr() << "ERROR: cannot change to unencapsulated representation for pixel data !" << endl;
                ofConsole.unlockCerr();
            }
        }
    } else {
        if (DicomImageClass::checkDebugLevel(DicomImageClass::DL_Errors))
        {
            ofConsole.lockCerr() << "ERROR: no pixel data found in DICOM dataset !" << endl;
            ofConsole.unlockCerr();
        }
    }
}


/*--------------*
 *  destructor  *
 *--------------*/

DiDocument::~DiDocument()
{
    /* DICOM image loaded from file: delete file format (and data set) */
    if (FileFormat != NULL)
        delete FileFormat;
    /* DICOM image loaded from external data set: only delete if flag is set */
    else if (Flags & CIF_TakeOverExternalDataset)
        delete Object;
}


/********************************************************************/

DcmElement *DiDocument::search(const DcmTagKey &tag,
                               DcmObject *obj) const
{
    DcmStack stack;
    if (obj == NULL)
        obj = Object;
    if ((obj != NULL) && (obj->search(tag, stack, ESM_fromHere, OFFalse /* searchIntoSub */) == EC_Normal) &&
        (stack.top()->getLength(Xfer) > 0))
    {
        return OFstatic_cast(DcmElement *, stack.top());
    }
    return NULL;
}


/********************************************************************/

int DiDocument::search(const DcmTagKey &tag,
                       DcmStack &pstack) const
{
    if (pstack.empty())
        pstack.push(Object);
    DcmObject *obj = pstack.top();
    if ((obj != NULL) && (obj->search(tag, pstack, ESM_fromHere, OFFalse /* searchIntoSub */) == EC_Normal) &&
        (pstack.top()->getLength(Xfer) > 0))
            return 1;
    return 0;
}


/********************************************************************/


unsigned long DiDocument::getVM(const DcmTagKey &tag) const
{
    DcmElement *elem = search(tag);
    if (elem != NULL)
        return elem->getVM();
    return 0;
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   Uint16 &returnVal,
                                   const unsigned long pos,
                                   DcmObject *item) const
{
    return getElemValue(search(tag, item), returnVal, pos);
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   Sint16 &returnVal,
                                   const unsigned long pos) const
{
    DcmElement *elem = search(tag);
    if (elem != NULL)
    {
        elem->getSint16(returnVal, pos);
        return elem->getVM();
    }
    return 0;
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   Uint32 &returnVal,
                                   const unsigned long pos) const
{
    DcmElement *elem = search(tag);
    if (elem != NULL)
    {
        elem->getUint32(returnVal, pos);
        return elem->getVM();
    }
    return 0;
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   Sint32 &returnVal,
                                   const unsigned long pos) const
{
    DcmElement *elem = search(tag);
    if (elem != NULL)
    {
        elem->getSint32(returnVal, pos);
        return elem->getVM();
    }
    return 0;
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   double &returnVal,
                                   const unsigned long pos) const
{
    DcmElement *elem = search(tag);
    if (elem != NULL)
    {
        elem->getFloat64(returnVal, pos);
        return elem->getVM();
    }
    return 0;
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   const Uint16 *&returnVal,
                                   DcmObject *item) const
{
    DcmElement *elem = search(tag, item);
    if (elem != NULL)
    {
        Uint16 *val;
        elem->getUint16Array(val);
        returnVal = val;
        if (elem->getVR() == EVR_OW)
            return elem->getLength(Xfer) / sizeof(Uint16);
        return elem->getVM();
    }
    return 0;
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   const char *&returnVal,
                                   DcmObject *item) const
{
    return getElemValue(search(tag, item), returnVal);
}


unsigned long DiDocument::getValue(const DcmTagKey &tag,
                                   OFString &returnVal,
                                   const unsigned long pos,
                                   DcmObject *item) const
{
    return getElemValue(search(tag, item), returnVal, pos);
}


unsigned long DiDocument::getSequence(const DcmTagKey &tag,
                                      DcmSequenceOfItems *&seq) const
{
    DcmElement *elem = search(tag);
    if ((elem != NULL) && (elem->ident() == EVR_SQ))
        return (seq = OFstatic_cast(DcmSequenceOfItems *, elem))->card();
    return 0;
}


unsigned long DiDocument::getElemValue(const DcmElement *elem,
                                       Uint16 &returnVal,
                                       const unsigned long pos)
{
    if (elem != NULL)
    {
        OFconst_cast(DcmElement *, elem)->getUint16(returnVal, pos);   // remove 'const' to use non-const methods
        return OFconst_cast(DcmElement *, elem)->getVM();
    }
    return 0;
}


unsigned long DiDocument::getElemValue(const DcmElement *elem,
                                       const Uint16 *&returnVal)
{
    if (elem != NULL)
    {
        Uint16 *val;                                            // parameter has no 'const' qualifier
        OFconst_cast(DcmElement *, elem)->getUint16Array(val); // remove 'const' to use non-const methods
        returnVal = val;
        if (OFconst_cast(DcmElement *, elem)->getVR() == EVR_OW)
            return OFconst_cast(DcmElement *, elem)->getLength(/*Xfer*/) / sizeof(Uint16);
        return OFconst_cast(DcmElement *, elem)->getVM();
    }
    return 0;
}


unsigned long DiDocument::getElemValue(const DcmElement *elem,
                                       const char *&returnVal)
{
    if (elem != NULL)
    {
        char *val;                                         // parameter has no 'const' qualifier
        OFconst_cast(DcmElement *, elem)->getString(val); // remove 'const' to use non-const methods
        returnVal = val;
        return OFconst_cast(DcmElement *, elem)->getVM();
    }
    return 0;
}


unsigned long DiDocument::getElemValue(const DcmElement *elem,
                                       OFString &returnVal,
                                       const unsigned long pos)
{
    if (elem != NULL)
    {
        OFconst_cast(DcmElement *, elem)->getOFString(returnVal, pos); // remove 'const' to use non-const methods
        return OFconst_cast(DcmElement *, elem)->getVM();
    }
    return 0;
}


/*
 *
 * CVS/RCS Log:
 * $Log: didocu.cc,v $
 * Revision 1.16  2005/12/08 15:42:48  meichel
 * Changed include path schema for all DCMTK header files
 *
 * Revision 1.15  2003/12/08 15:13:18  joergr
 * Adapted type casts to new-style typecast operators defined in ofcast.h.
 *
 * Revision 1.14  2002/08/21 09:51:47  meichel
 * Removed DicomImage and DiDocument constructors that take a DcmStream
 *   parameter
 *
 * Revision 1.13  2002/06/26 16:10:15  joergr
 * Added new methods to get the explanation string of stored VOI windows and
 * LUTs (not only of the currently selected VOI transformation).
 * Added configuration flag that enables the DicomImage class to take the
 * responsibility of an external DICOM dataset (i.e. delete it on destruction).
 *
 * Revision 1.12  2001/11/29 16:59:52  joergr
 * Fixed bug in dcmimgle that caused incorrect decoding of some JPEG compressed
 * images (certain DICOM attributes, e.g. photometric interpretation, might
 * change during decompression process).
 *
 * Revision 1.11  2001/11/19 12:57:04  joergr
 * Adapted code to support new dcmjpeg module and JPEG compressed images.
 *
 * Revision 1.10  2001/06/01 15:49:54  meichel
 * Updated copyright header
 *
 * Revision 1.9  2000/09/12 10:06:14  joergr
 * Corrected bug: wrong parameter for attribute search routine led to crashes
 * when multiple pixel data attributes were contained in the dataset (e.g.
 * IconImageSequence).
 *
 * Revision 1.8  2000/04/28 12:33:42  joergr
 * DebugLevel - global for the module - now derived from OFGlobal (MF-safe).
 *
 * Revision 1.7  2000/04/27 13:10:26  joergr
 * Dcmimgle library code now consistently uses ofConsole for error output.
 *
 * Revision 1.6  2000/03/08 16:24:27  meichel
 * Updated copyright header.
 *
 * Revision 1.5  2000/03/03 14:09:18  meichel
 * Implemented library support for redirecting error messages into memory
 *   instead of printing them to stdout/stderr for GUI applications.
 *
 * Revision 1.4  1999/04/28 15:01:43  joergr
 * Introduced new scheme for the debug level variable: now each level can be
 * set separately (there is no "include" relationship).
 *
 * Revision 1.3  1998/12/16 16:11:14  joergr
 * Added methods to use getOFString from class DcmElement (incl. multi value
 * fields).
 *
 * Revision 1.2  1998/12/14 17:33:45  joergr
 * Added (simplified) methods to return values of a given DcmElement object.
 *
 * Revision 1.1  1998/11/27 15:53:36  joergr
 * Added copyright message.
 * Added static methods to return the value of a given element.
 * Added support of frame start and count for future use (will be explained
 * later if it is fully implemented).
 *
 * Revision 1.9  1998/07/01 08:39:35  joergr
 * Minor changes to avoid compiler warnings (gcc 2.8.1 with additional
 * options), e.g. add copy constructors.
 *
 * Revision 1.8  1998/06/25 08:52:05  joergr
 * Added compatibility mode to support ACR-NEMA images and wrong
 * palette attribute tags.
 *
 * Revision 1.7  1998/05/11 14:52:28  joergr
 * Added CVS/RCS header to each file.
 *
 *
 */
