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
 *  Author:  Gerd Ehlers
 *
 *  Purpose: class DcmItem
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:41:16 $
 *  CVS/RCS Revision: $Revision: 1.97 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */
#include "dcmtk/dcmdata/dcitem.h"

#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#define INCLUDE_CCTYPE
#include "dcmtk/ofstd/ofstdinc.h"

#include "dcmtk/dcmdata/dcdebug.h"
#include "dcmtk/dcmdata/dcdefine.h"    /* for memzero() */
#include "dcmtk/dcmdata/dcdeftag.h"    /* for name constants */
#include "dcmtk/dcmdata/dcistrma.h"    /* for class DcmInputStream */
#include "dcmtk/dcmdata/dcobject.h"
#include "dcmtk/dcmdata/dcostrma.h"    /* for class DcmOutputStream */
#include "dcmtk/dcmdata/dcovlay.h"
#include "dcmtk/dcmdata/dcpixel.h"
#include "dcmtk/dcmdata/dcsequen.h"
#include "dcmtk/dcmdata/dcswap.h"
#include "dcmtk/dcmdata/dcvr.h"
#include "dcmtk/dcmdata/dcvrae.h"
#include "dcmtk/dcmdata/dcvras.h"
#include "dcmtk/dcmdata/dcvrat.h"
#include "dcmtk/dcmdata/dcvrcs.h"
#include "dcmtk/dcmdata/dcvrda.h"
#include "dcmtk/dcmdata/dcvrds.h"
#include "dcmtk/dcmdata/dcvrdt.h"
#include "dcmtk/dcmdata/dcvrfd.h"
#include "dcmtk/dcmdata/dcvrfl.h"
#include "dcmtk/dcmdata/dcvris.h"
#include "dcmtk/dcmdata/dcvrlo.h"
#include "dcmtk/dcmdata/dcvrlt.h"
#include "dcmtk/dcmdata/dcvrobow.h"
#include "dcmtk/dcmdata/dcvrof.h"
#include "dcmtk/dcmdata/dcvrpn.h"
#include "dcmtk/dcmdata/dcvrsh.h"
#include "dcmtk/dcmdata/dcvrsl.h"
#include "dcmtk/dcmdata/dcvrss.h"
#include "dcmtk/dcmdata/dcvrst.h"
#include "dcmtk/dcmdata/dcvrtm.h"
#include "dcmtk/dcmdata/dcvrui.h"
#include "dcmtk/dcmdata/dcvrul.h"
#include "dcmtk/dcmdata/dcvrulup.h"
#include "dcmtk/dcmdata/dcvrus.h"
#include "dcmtk/dcmdata/dcvrut.h"
#include "dcmtk/dcmdata/dcxfer.h"
#include "dcmtk/ofstd/ofstream.h"
#include "dcmtk/ofstd/ofstring.h"
#include "dcmtk/ofstd/ofcast.h"


// ********************************


DcmItem::DcmItem()
  : DcmObject(ItemTag),
    elementList(NULL),
    lastElementComplete(OFTrue),
    fStartPosition(0),
    privateCreatorCache()
{
    elementList = new DcmList;
}


DcmItem::DcmItem(const DcmTag &tag,
                 const Uint32 len)
  : DcmObject(tag, len),
    elementList(NULL),
    lastElementComplete(OFTrue),
    fStartPosition(0),
    privateCreatorCache()
{
    elementList = new DcmList;
}


DcmItem::DcmItem(const DcmItem &old)
  : DcmObject(old),
    elementList(new DcmList),
    lastElementComplete(old.lastElementComplete),
    fStartPosition(old.fStartPosition),
    privateCreatorCache()
{
    if (!old.elementList->empty())
    {
        elementList->seek(ELP_first);
        old.elementList->seek(ELP_first);
        do 
        {
            elementList->insert(old.elementList->get()->clone(), ELP_next);
        } while (old.elementList->seek(ELP_next));
    }
}


DcmItem::~DcmItem()
{
    DcmObject *dO;
    elementList->seek(ELP_first);
    while (!elementList->empty())
    {
        dO = elementList->remove();
        delete dO;
    }
    delete elementList;
}


// ********************************


OFBool DcmItem::foundVR(char *atposition)
{
    char c1 =  atposition[0];
    char c2 = atposition[1];
    OFBool valid = OFFalse;

    if (isalpha(c1) && isalpha(c2))
    {
        char vrName[3];
        vrName[0] = c1;
        vrName[1] = c2;
        vrName[2] = '\0';

        /* is this VR name a standard VR descriptor */
        DcmVR vr(vrName);
        valid = vr.isStandard();
    } else {
        /* cannot be a valid VR name since non-characters */
        valid = OFFalse;
    }
    return valid;
}


// ********************************


E_TransferSyntax DcmItem::checkTransferSyntax(DcmInputStream & inStream)
{
    E_TransferSyntax transferSyntax;
    char tagAndVR[6];

    /* read 6 bytes from the input stream (try to read tag and VR (data type)) */
    inStream.mark();
    inStream.read(tagAndVR, 6);               // check Tag & VR
    inStream.putback();

    /* create two tag variables (one for little, one for big */
    /* endian) in order to figure out, if there is a valid tag */
    char c1 = tagAndVR[0];
    char c2 = tagAndVR[1];
    char c3 = tagAndVR[2];
    char c4 = tagAndVR[3];
    Uint16 t1 = OFstatic_cast(unsigned short, (c1 & 0xff) + ((c2 & 0xff) << 8));  // explicit little endian
    Uint16 t2 = OFstatic_cast(unsigned short, (c3 & 0xff) + ((c4 & 0xff) << 8));  // conversion
    DcmTag taglittle(t1, t2);
    DcmTag tagbig(swapShort(t1), swapShort(t2));

    /* now we want to determine the transfer syntax which was used to code the information in the stream. */
    /* The decision is based on two questions: a) Did we encounter a valid tag? and b) Do the last 2 bytes */
    /* which were read from the stream represent a valid VR? In certain special cases, where the transfer */
    /* cannot be determined without doubt, we want to guess the most probable transfer syntax. */

    /* if both tag variables show an error, we encountered an invalid tag */
    if ((taglittle.error().bad()) && (tagbig.error().bad()))
    {
        /* in case we encounterd an invalid tag, we want to assume that the used transfer syntax */
        /* is a little endian transfer syntax. Now we have to figure out, if it is an implicit or */
        /* explicit transfer syntax. Hence, check if the last 2 bytes represent a valid VR. */
        if (foundVR(&tagAndVR[4]))
        {
            /* if the last 2 bytes represent a valid VR, we assume that the used */
            /* transfer syntax is the little endian explicit transfer syntax. */
            transferSyntax = EXS_LittleEndianExplicit;
        } else {
            /* if the last 2 bytes did not represent a valid VR, we assume that the */
            /* used transfer syntax is the little endian implicit transfer syntax. */
            transferSyntax = EXS_LittleEndianImplicit;
        }
    }
    /* if at least one tag variable did not show an error, we encountered a valid tag */
    else
    {
        /* in case we encounterd a valid tag, we want to figure out, if it is an implicit or */
        /* explicit transfer syntax. Hence, check if the last 2 bytes represent a valid VR. */
        if (foundVR(&tagAndVR[4]))
        {
            /* having figured out that the last 2 bytes represent a valid */
            /* VR, we need to find out which of the two tags was valid */
            if (taglittle.error().bad())
            {
                /* if the litte endian tag was invalid, the transfer syntax is big endian explicit */
                transferSyntax = EXS_BigEndianExplicit;
            }
            else if (tagbig.error().bad())
            {
                /* if the big endian tag was invalid, the transfer syntax is little endian explicit */
                transferSyntax = EXS_LittleEndianExplicit;
            } else {
                /* if both tags were valid, we take a look at the group numbers. Since */
                /* group 0008 is much more probable than group 0800 for the first tag */
                /* we specify the following: */
                if ((taglittle.getGTag() > 0xff)&&(tagbig.getGTag() <= 0xff)) transferSyntax = EXS_BigEndianExplicit;
                else transferSyntax = EXS_LittleEndianExplicit;
            }
        } else {
            /* having figured out that the last 2 bytes do not represent a */
            /* valid VR, we need to find out which of the two tags was valid */
            if (taglittle.error().bad())
            {
                /* if the litte endian tag was invalid, the transfer syntax is big endian implicit */
                transferSyntax = EXS_BigEndianImplicit;
            }
            else if (tagbig.error().bad())
            {
                /* if the big endian tag was invalid, the transfer syntax is little endian implicit */
                transferSyntax = EXS_LittleEndianImplicit;
            } else {
                /* if both tags were valid, we take a look at the group numbers. Since */
                /* group 0008 is much more probable than group 0800 for the first tag */
                /* we specify the following: */
                if ((taglittle.getGTag() > 0xff)&&(tagbig.getGTag() <= 0xff)) transferSyntax = EXS_BigEndianImplicit;
                else transferSyntax = EXS_LittleEndianImplicit;
            }
        }
    }
    /* dump information on a certain debug level */
    DCM_dcmdataDebug(3, ("found TransferSyntax=(%s)", DcmXfer(transferSyntax).getXferName()));

    /* return determined transfer syntax */
    return transferSyntax;
} // DcmItem::checkTransferSyntax


// ********************************


DcmEVR DcmItem::ident() const
{
    return EVR_item;
}


unsigned long DcmItem::getVM()
{
    return 1;
}


// ********************************


void DcmItem::print(ostream &out,
                    const size_t flags,
                    const int level,
                    const char *pixelFileName,
                    size_t *pixelCounter)
{
    /* print item start line */
    if (flags & DCMTypes::PF_showTreeStructure)
    {
        /* empty text */
        printInfoLine(out, flags, level);
        /* print item content */
        if (!elementList->empty())
        {
            /* reset internal flags */
            const size_t newFlags = flags & ~DCMTypes::PF_lastEntry;
            /* print attributes on this level */
            DcmObject *dO;
            elementList->seek(ELP_first);
            OFBool ok;
            do {
                dO = elementList->get();
                ok = (elementList->seek(ELP_next) != NULL);
                if (!ok)
                    dO->print(out, newFlags | DCMTypes::PF_lastEntry, level + 1, pixelFileName, pixelCounter);
                else
                    dO->print(out, newFlags, level + 1, pixelFileName, pixelCounter);
            } while (ok);
        }
    } else {
        OFOStringStream oss;
        oss << "(Item with ";
        if (Length == DCM_UndefinedLength)
            oss << "undefined";
        else
            oss << "explicit";
        oss << " length #=" << card() << ")" << OFStringStream_ends;
        OFSTRINGSTREAM_GETSTR(oss, tmpString)
        printInfoLine(out, flags, level, tmpString);
        OFSTRINGSTREAM_FREESTR(tmpString)
        /* print item content */
        if (!elementList->empty())
        {
            DcmObject *dO;
            elementList->seek(ELP_first);
            do {
                dO = elementList->get();
                dO->print(out, flags, level + 1, pixelFileName, pixelCounter);
            } while (elementList->seek(ELP_next));
        }
        /* print item end line */
        DcmTag delimItemTag(DCM_ItemDelimitationItem);
        if (Length == DCM_UndefinedLength)
            printInfoLine(out, flags, level, "(ItemDelimitationItem)", &delimItemTag);
        else
            printInfoLine(out, flags, level, "(ItemDelimitationItem for re-encoding)", &delimItemTag);
    }
}


// ********************************


OFCondition DcmItem::writeXML(ostream &out,
                              const size_t flags)
{
    /* XML start tag for "item" */
    out << "<item";
    /* cardinality (number of attributes) = 1..n */
    out << " card=\"" << card() << "\"";
    /* value length in bytes = 0..max (if not undefined) */
    if (Length != DCM_UndefinedLength)
        out << " len=\"" << Length << "\"";
    out << ">" << endl;
    /* write item content */
    if (!elementList->empty())
    {
        /* write content of all children */
        DcmObject *dO;
        elementList->seek(ELP_first);
        do {
            dO = elementList->get();
            dO->writeXML(out, flags);
        } while (elementList->seek(ELP_next));
    }
    /* XML end tag for "item" */
    out << "</item>" << endl;
    /* always report success */
    return EC_Normal;
}


// ********************************


OFBool DcmItem::canWriteXfer(const E_TransferSyntax newXfer,
                             const E_TransferSyntax oldXfer)
{
    OFBool canWrite = OFTrue;
    if (newXfer == EXS_Unknown)
        canWrite = OFFalse;
    else if (!elementList->empty())
    {
        DcmObject *dO;
        elementList->seek(ELP_first);
        do {
            dO = elementList->get();
            canWrite = dO->canWriteXfer(newXfer, oldXfer);
        } while (elementList->seek(ELP_next) && canWrite);
    }
    return canWrite;
}


// ********************************


Uint32 DcmItem::calcElementLength(const E_TransferSyntax xfer,
                                  const E_EncodingType enctype)
{
    Uint32 itemlen = 0;
    DcmXfer xferSyn(xfer);
    itemlen = getLength(xfer, enctype) + xferSyn.sizeofTagHeader(getVR());
    if (enctype == EET_UndefinedLength)
        itemlen += 8;
    return itemlen;
}


// ********************************


Uint32 DcmItem::getLength(const E_TransferSyntax xfer,
                          const E_EncodingType enctype)
{
    Uint32 itemlen = 0;
    if (!elementList->empty())
    {
        DcmObject *dO;
        elementList->seek(ELP_first);
        do {
            dO = elementList->get();
            itemlen += dO->calcElementLength(xfer, enctype);
        } while (elementList->seek(ELP_next));
    }
    return itemlen;
}


// ********************************


OFCondition DcmItem::computeGroupLengthAndPadding(const E_GrpLenEncoding glenc,
                                                  const E_PaddingEncoding padenc,
                                                  const E_TransferSyntax xfer,
                                                  const E_EncodingType enctype,
                                                  const Uint32 padlen,
                                                  const Uint32 subPadlen,
                                                  Uint32 instanceLength)
{
    /* if certain conditions are met, this is considered to be an illegal call. */
    if ((padenc == EPD_withPadding && (padlen % 2 || subPadlen % 2)) ||
        ((glenc == EGL_recalcGL || glenc == EGL_withGL ||
          padenc == EPD_withPadding) && xfer == EXS_Unknown))
        return EC_IllegalCall;

    /* if the caller specified that group length tags and padding */
    /* tags are not supposed to be changed, there is nothing to do. */
    if (glenc == EGL_noChange && padenc == EPD_noChange)
        return EC_Normal;

    /* if we get to this point, we need to do something. First of all, set the error indicator to normal. */
    OFCondition l_error = EC_Normal;

    /* if there are elements in this item... */
    if (!elementList->empty())
    {
        /* initialize some variables */
        DcmObject *dO;
        OFBool beginning = OFTrue;
        Uint16 lastGrp = 0x0000;
        Uint16 actGrp;
        DcmUnsignedLong * actGLElem = NULL;
        DcmUnsignedLong * paddingGL = NULL;
        Uint32 grplen = 0;
        DcmXfer xferSyn(xfer);

        /* determine the current seek mode and set the list pointer to the first element */
        E_ListPos seekmode = ELP_next;
        elementList->seek(ELP_first);

        /* start a loop: we want to go through all elements as long as everything is okay */
        do
        {
            /* set the seek mode to "next" again, in case it has been modified in the last iteration */
            seekmode = ELP_next;

            /* get the current element and assign it to a local variable */
            dO = elementList->get();

            /* if the current element is a sequence, compute group length and padding for the sub sequence */
            if (dO->getVR() == EVR_SQ)
            {
                Uint32 templen = instanceLength + xferSyn.sizeofTagHeader(EVR_SQ);
                l_error =
                    OFstatic_cast(DcmSequenceOfItems *, dO)->computeGroupLengthAndPadding
                    (glenc, padenc, xfer, enctype, subPadlen, subPadlen,
                     templen);
            }

            /* if everything is ok so far */
            if (l_error.good())
            {
                /* in case one of the following two conditions is met */
                /*  (i) the caller specified that we want to add or remove group length elements and the current */
                /*      element's tag shows that it is a group length element (tag's element number equals 0x0000) */
                /*  (ii) the caller specified that we want to add or remove padding elements and the current */
                /*      element's tag shows that it is a padding element (tag is (0xfffc,0xfffc) */
                /* then we want to delete the current (group length or padding) element */
                if (((glenc ==  EGL_withGL || glenc == EGL_withoutGL) && dO->getETag() == 0x0000) ||
                    (padenc != EPD_noChange && dO->getTag() == DCM_DataSetTrailingPadding))
                {
                    delete elementList->remove();
                    seekmode = ELP_atpos;           // remove = 1 forward
                    dO = NULL;
                }
                /* if the above mentioned conditions are not met but the caller specified that we want to add group */
                /* length tags for every group or that we want to recalculate values for existing group length tags */
                else  if (glenc == EGL_withGL || glenc == EGL_recalcGL)
                {
                    /* we need to determine the current element's group number */
                    actGrp = dO->getGTag();

                    /* and if the group number is different from the last remembered group number or */
                    /* if this id the very first element that is treated then we've found a new group */
                    if (actGrp!=lastGrp || beginning) // new Group found
                    {
                        /* set beginning to false in order to specify that the */
                        /* very first element has already been treated */
                        beginning = OFFalse;

                        /* if the current element is a group length element) and it's data type */
                        /* is not UL replace this element with one that has a UL datatype since */
                        /* group length elements are supposed to have this data type */
                        if (dO->getETag() == 0x0000 && dO->ident() != EVR_UL)
                        {
                            delete elementList->remove();
                            DcmTag tagUL(actGrp, 0x0000, EVR_UL);
                            DcmUnsignedLong *dUL = new DcmUnsignedLong(tagUL);
                            elementList->insert(dUL, ELP_prev);
                            dO = dUL;
                            ofConsole.lockCerr() << "DcmItem: Group Length with VR other than UL found, corrected." << endl;
                            ofConsole.unlockCerr();
                        }
                        /* if the above mentioned condition is not met but the caller specified */
                        /* that we want to add group length elements, we need to add such an element */
                        else if (glenc == EGL_withGL)
                        {
                            // Create GroupLength element
                            DcmTag tagUL(actGrp, 0x0000, EVR_UL);
                            DcmUnsignedLong *dUL = new DcmUnsignedLong(tagUL);
                            // insert new GroupLength element
                            elementList->insert(dUL, ELP_prev);
                            dO = dUL;
                        }

                        /* in case we want to add padding elements and the current element is a */
                        /* padding element we want to remember the padding element so that the */
                        /* group length of this element can be stored later */
                        if (padenc == EPD_withPadding && actGrp == 0xfffc)
                            paddingGL = OFstatic_cast(DcmUnsignedLong *, dO);

                        /* if actGLElem conatins a valid pointer it was set in one of the last iterations */
                        /* to the group lenght element of the last group. We need to write the current computed */
                        /* group length value to this element. */
                        if (actGLElem != NULL)
                        {
                            actGLElem->putUint32(grplen);
                            DCM_dcmdataDebug(2, ("DcmItem::computeGroupLengthAndPadding() Length of Group 0x%4.4x len=%lu", actGLElem->getGTag(), grplen));
                        }

                        /* set the group length value to 0 since it is the beginning of the new group */
                        grplen = 0;

                        /* if the current element is a group length element, remember its address for later */
                        /* (we need to assign the group length value to this element in a subsequent iteration) */
                        /* in case the current element (at the beginning of the group) is not a group length */
                        /* element, set the actGLElem pointer to NULL. */
                        if (dO->getETag() == 0x0000)
                            actGLElem = OFstatic_cast(DcmUnsignedLong *, dO);
                        else
                            actGLElem = NULL;
                    }
                    /* if this is not a new group, calculate the element's length and add it */
                    /* to the currently computed group length value */
                    else
                        grplen += dO->calcElementLength(xfer, enctype);

                    /* remember the current element's group number so that it is possible to */
                    /* figure out if a new group is treated in the following iteration */
                    lastGrp = actGrp;
                }
            }
        } while (l_error.good() && elementList->seek(seekmode));

        /* if there was no error and the caller specified that we want to add or recalculate */
        /* group length tags and if actGLElem has a valid value, we need to add the above */
        /* computed group length value to the last group's group length element */
        if (l_error.good() && (glenc == EGL_withGL || glenc == EGL_recalcGL) && actGLElem)
            actGLElem->putUint32(grplen);

        /* if the caller specified that we want to add padding elements and */
        /* if the length up to which shall be padded does not equal 0 we might */
        /* have to add a padding element */
        if (padenc == EPD_withPadding && padlen)
        {
            /* calculate how much space the entire padding element is supposed to occupy */
            Uint32 padding;
            if (ident() == EVR_dataset)
            {
                instanceLength += calcElementLength(xfer, enctype);
                padding = padlen - (instanceLength % padlen);
            } else
                padding = padlen - (getLength(xfer, enctype) % padlen);

            /* if now padding does not equal padlen we need to create a padding element. (if both values are equal */
            /* the element does have the exact required padlen length and does not need a padding element.) */
            if (padding != padlen)
            {
                /* Create new padding element */
                DcmOtherByteOtherWord * paddingEl = new DcmOtherByteOtherWord(DCM_DataSetTrailingPadding);

                /* calculate the length of the new element */
                Uint32 tmplen = paddingEl->calcElementLength(xfer, enctype);

                /* in case padding is smaller than the header of the padding element, we */
                /* need to increase padding (the value which specifies how much space the */
                /* entire padding element is supposed to occupy) until it is no longer smaller */
                while (tmplen > padding)
                    padding += padlen;

                /* determine the amount of bytes that have to be added to the */
                /* padding element so that it has the correct size */
                padding -= tmplen;

                /* create an array of a corresponding size and set the arrayfields */
                Uint8 * padBytes = new Uint8[padding];
                memzero(padBytes, size_t(padding));

                /* set information in the above created padding element (size and actual value) */
                paddingEl->putUint8Array(padBytes, padding);

                /* delete the above created array */
                delete[] padBytes;

                /* insert the padding element into this */
                insert(paddingEl);

                /* finally we need to update the group length for the padding element if it exists */
                if (paddingGL)
                {
                    Uint32 len;
                    paddingGL->getUint32(len);
                    len += paddingEl->calcElementLength(xfer, enctype);
                    paddingGL->putUint32(len);
                }
            }
        }
    }
    return l_error;
}


// ********************************


OFCondition DcmItem::readTagAndLength(DcmInputStream &inStream,
                                      const E_TransferSyntax xfer,
                                      DcmTag &tag,
                                      Uint32 &length,
                                      Uint32 &bytesRead)
{
    OFCondition l_error = EC_Normal;
    Uint32 valueLength = 0;
    DcmEVR nxtobj = EVR_UNKNOWN;
    Uint16 groupTag = 0xffff;
    Uint16 elementTag = 0xffff;

    /* Create a DcmXfer object based on the transfer syntax which was passed */
    DcmXfer xferSyn(xfer);

    /* dump some information if required */
    DCM_dcmdataDebug(4, ("DcmItem::readTagAndLength() read transfer syntax %s", xferSyn.getXferName()));

    /* bail out if at end of stream */
    if (inStream.eos())
        return EC_EndOfStream;

    /* check if either 4 (for implicit transfer syntaxes) or 6 (for explicit transfer */
    /* syntaxes) bytes are available in (i.e. can be read from) inStream. if an error */
    /* occured while performing this check return this error */
    if (inStream.avail() < (xferSyn.isExplicitVR() ? 6u:4u))
        return EC_StreamNotifyClient;

    /* determine the byte ordering of the transfer syntax which was passed; */
    /* if the byte ordering is unknown, this is an illegal call. */
    const E_ByteOrder byteOrder = xferSyn.getByteOrder();
    if (byteOrder == EBO_unknown)
        return EC_IllegalCall;

    /* read tag information (4 bytes) from inStream and create a corresponding DcmTag object */
    inStream.mark();
    inStream.read(&groupTag, 2);
    inStream.read(&elementTag, 2);
    swapIfNecessary(gLocalByteOrder, byteOrder, &groupTag, 2, 2);
    swapIfNecessary(gLocalByteOrder, byteOrder, &elementTag, 2, 2);
    // tag has been read
    bytesRead = 4;
    DcmTag newTag(groupTag, elementTag);

    /* if the transfer syntax which was passed is an explicit VR syntax and if the current */
    /* item is not a delimitation item (note that delimitation items do not have a VR), go */
    /* ahead and read 2 bytes from inStream. These 2 bytes contain this item's VR value. */
    if (xferSyn.isExplicitVR() && newTag.getEVR() != EVR_na)
    {
        char vrstr[3];
        vrstr[2] = '\0';

        /* read 2 bytes */
        inStream.read(vrstr, 2);

        /* create a corresponding DcmVR object */
        DcmVR vr(vrstr);

        /* if the VR which was read is not a standard VR, print a warning */
        if (!vr.isStandard())
        {
            ostream &localCerr = ofConsole.lockCerr();
            localCerr << "DcmItem: Non-standard VR '" << vrstr
                      << "' encountered while parsing attribute " << newTag.getXTag() << ", assuming ";
            if (vr.usesExtendedLengthEncoding())
                localCerr << "4 byte length field" << endl;
            else
                localCerr << "2 byte length field" << endl;
            ofConsole.unlockCerr();
        }

        /* set the VR which was read in the above created tag object. */
        newTag.setVR(vr);

        /* increase counter by 2 */
        bytesRead += 2;
    }

    /* special handling for private elements */
    if ((newTag.getGroup() & 1) && (newTag.getElement() >= 0x1000))
    {
        const char *pc = privateCreatorCache.findPrivateCreator(newTag);
        if (pc)
        {
            // we have a private creator for this element
            newTag.setPrivateCreator(pc);

            if (xferSyn.isImplicitVR())
            {
                // try to update VR from dictionary now that private creator is known
                newTag.lookupVRinDictionary();
            }
        }
    }

    /* determine this item's VR */
    nxtobj = newTag.getEVR();

    /* the next thing we want to do is read the value in the length field from inStream. */
    /* determine if there is a corresponging amount of bytes (for the length field) still */
    /* available in inStream. if not, return an error. */
    if (inStream.avail() < xferSyn.sizeofTagHeader(nxtobj) - bytesRead)
    {
        inStream.putback();    // the UnsetPutbackMark is in readSubElement
        bytesRead = 0;
        l_error = EC_StreamNotifyClient;
        return l_error;
    }

    /* read the value in the length field. In some cases, it is 4 bytes wide, in other */
    /* cases only 2 bytes (see DICOM standard (year 2000) part 5, section 7.1.1) (or the */
    /* corresponding section in a later version of the standard) */
    if (xferSyn.isImplicitVR() || nxtobj == EVR_na)   //note that delimitation items don't have a VR
    {
        inStream.read(&valueLength, 4);            //length field is 4 bytes wide
        swapIfNecessary(gLocalByteOrder, byteOrder, &valueLength, 4, 4);
        bytesRead += 4;
    } else {                                       //the transfer syntax is explicit VR
        DcmVR vr(newTag.getEVR());
        if (vr.usesExtendedLengthEncoding())
        {
            Uint16 reserved;
            inStream.read(&reserved, 2);           // 2 reserved bytes
            inStream.read(&valueLength, 4);        // length field is 4 bytes wide
            swapIfNecessary(gLocalByteOrder, byteOrder, &valueLength, 4, 4);
            bytesRead += 6;
        } else {
            Uint16 tmpValueLength;
            inStream.read(&tmpValueLength, 2);     // length field is 2 bytes wide
            swapIfNecessary(gLocalByteOrder, byteOrder, &tmpValueLength, 2, 2);
            bytesRead += 2;
            valueLength = tmpValueLength;
        }
    }
    /* if the value in length is odd, print an error message */
    if ((valueLength & 1)&&(valueLength != OFstatic_cast(Uint32, -1)))
    {
        ofConsole.lockCerr() << "DcmItem: Length of attribute " << newTag << " is odd" << endl;
        ofConsole.unlockCerr();
    }
    /* assign values to out parameter */
    length = valueLength;
    tag = newTag;

    /* return return value */
    return l_error;
}


// ********************************


OFCondition DcmItem::readSubElement(DcmInputStream &inStream,
                                    DcmTag &newTag,
                                    const Uint32 newLength,
                                    const E_TransferSyntax xfer,
                                    const E_GrpLenEncoding glenc,
                                    const Uint32 maxReadLength)
{
    DcmElement *subElem = NULL;

    /* create a new DcmElement* object with corresponding tag and */
    /* length; the object will be accessible through subElem */
    OFBool readAsUN = OFFalse;
    OFCondition l_error = newDicomElement(subElem, newTag, newLength, &privateCreatorCache, readAsUN);

    /* if no error occured and subElem does not equal NULL, go ahead */
    if (l_error.good() && subElem != NULL)
    {
        // inStream.UnsetPutbackMark(); // not needed anymore with new stream architecture

        /* insert the new element into the (sorted) element list and */
        /* assign information which was read from the instream to it */
        subElem->transferInit();
        /* we need to read the content of the attribute, no matter if */
        /* inserting the attribute succeeds or fails */
        l_error = subElem->read(inStream, (readAsUN ? EXS_LittleEndianImplicit : xfer), glenc, maxReadLength);
        // try to insert element into item. Note that
        // "elementList->insert(subElem, ELP_next)" would be faster,
        // but this is better since this insert-function creates a
        // sorted element list.
        // We insert the element even if subElem->read() reported an error
        // because otherwise I/O suspension would fail.
        OFCondition temp_error = insert(subElem, OFFalse, OFTrue);

        if (temp_error.bad())
        {
            // produce diagnostics
            ofConsole.lockCerr() << "DcmItem: Element " << newTag
               << " found twice in one dataset/item, ignoring second entry" << endl;
            ofConsole.unlockCerr();
            delete subElem;
        }
    }
    /* else if an error occured, try to recover from this error */
    else if (l_error == EC_InvalidTag)
    {
        /* This is the second Putback operation on the putback mark in */
        /* readTagAndLength but it is impossible that both can be executed */
        /* without setting the Mark twice. */
        inStream.putback();
        ofConsole.lockCerr() << "DcmItem: Parse error while parsing attribute " <<  newTag << endl;
        ofConsole.unlockCerr();
    }
    else if (l_error != EC_ItemEnd)
    {
        // inStream.UnsetPutbackMark(); // not needed anymore with new stream architecture
        ofConsole.lockCerr() << "DcmItem: Parse error in sequence item, found " << newTag << " instead of an item delimiter" << endl;
        ofConsole.unlockCerr();
    } else {
        // inStream.UnsetPutbackMark(); // not needed anymore with new stream architecture
    }

    /* return result value */
    return l_error;
}


// ********************************


OFCondition DcmItem::read(DcmInputStream & inStream,
                          const E_TransferSyntax xfer,
                          const E_GrpLenEncoding glenc,
                          const Uint32 maxReadLength)
{
    /* check if this is an illegal call; if so set the error flag and do nothing, else go ahead */
    if (fTransferState == ERW_notInitialized)
        errorFlag = EC_IllegalCall;
    else
    {
        /* figure out if the stream reported an error */
        errorFlag = inStream.status();
        /* if the stream reported an error or if it is the end of the */
        /* stream, set the error flag correspondingly; else go ahead */
        if (errorFlag.good() && inStream.eos())
            errorFlag = EC_EndOfStream;
        else if (errorFlag.good() && fTransferState != ERW_ready)
        {
            /* if the transfer state of this item is ERW_init, get its start */
            /* position in the stream and set the transfer state to ERW_inWork */
            if (fTransferState == ERW_init)
            {
                fStartPosition = inStream.tell();  // start position of this item
                fTransferState = ERW_inWork;
            }
            DcmTag newTag;
            /* start a loop in order to read all elements (attributes) which are contained in the inStream */
            while (inStream.good() && (fTransferredBytes < Length || !lastElementComplete))
            {
                /* initialize variables */
                Uint32 newValueLength = 0;
                Uint32 bytes_tagAndLen = 0;
                /* if the reading of the last element was complete, go ahead and read the next element */
                if (lastElementComplete)
                {
                    /* read this element's tag and length information (and */
                    /* possibly also VR information) from the inStream */
                    errorFlag = readTagAndLength(inStream, xfer, newTag, newValueLength, bytes_tagAndLen);
                    /* increase counter correpsondingly */
                    fTransferredBytes += bytes_tagAndLen;

                    /* if there was an error while we were reading from the stream, terminate the while-loop */
                    /* (note that if the last element had been read from the instream in the last iteration, */
                    /* another iteration will be started, and of course then readTagAndLength(...) above will */
                    /* return that it encountered the end of the stream. It is only then (and here) when the */
                    /* while loop will be terminated.) */
                    if (errorFlag.bad())
                        break;
                    /* If we get to this point, we just started reading the first part */
                    /* of an element; hence, lastElementComplete is not longer true */
                    lastElementComplete = OFFalse;
                    /* read the actual data value which belongs to this element */
                    /* (attribute) and insert this information into the elementList */
                    errorFlag = readSubElement(inStream, newTag, newValueLength, xfer, glenc, maxReadLength);
                    /* if reading was successful, we read the entire data value information */
                    /* for this element; hence lastElementComplete is true again */
                    if (errorFlag.good())
                        lastElementComplete = OFTrue;
                } else {
                    /* if lastElementComplete is false, we have only read the current element's */
                    /* tag and length (and possibly VR) information as well as maybe some data */
                    /* data value information. We need to continue reading the data value */
                    /* information for this particular element. */
                    errorFlag = elementList->get()->read(inStream, xfer, glenc, maxReadLength);
                    /* if reading was successful, we read the entire information */
                    /* for this element; hence lastElementComplete is true */
                    if (errorFlag.good())
                        lastElementComplete = OFTrue;
                }
                /* remember how many bytes were read */
                fTransferredBytes = inStream.tell() - fStartPosition;
                if (errorFlag.good())
                {
                    // If we completed one element, update the private tag cache.
                    if (lastElementComplete)
                        privateCreatorCache.updateCache(elementList->get());
                } else
                    break; // if some error was encountered terminate the while-loop
            } //while

            /* determine an appropriate result value; note that if the above called read function */
            /* encountered the end of the stream before all information for this element could be */
            /* read from the stream, the errorFlag has already been set to EC_StreamNotifyClient. */
            if ((fTransferredBytes < Length || !lastElementComplete) && errorFlag.good())
                errorFlag = EC_StreamNotifyClient;
            if (errorFlag.good() && inStream.eos())
                errorFlag = EC_EndOfStream;
        } // else errorFlag
        /* modify the result value: two kinds of special error codes do not count as an error */
        if (errorFlag == EC_ItemEnd || errorFlag == EC_EndOfStream)
            errorFlag = EC_Normal;
        /* if at this point the error flag indicates success, the item has */
        /* been read completely; hence, set the transfer state to ERW_ready. */
        /* Note that all information for this element could be read from the */
        /* stream, the errorFlag is still set to EC_StreamNotifyClient. */
        if (errorFlag.good())
            fTransferState = ERW_ready;
    }
    /* return result value */
    return errorFlag;
} // DcmItem::read()


// ********************************


OFCondition DcmItem::write(DcmOutputStream &outStream,
                           const E_TransferSyntax oxfer,
                           const E_EncodingType enctype)
{
  if (fTransferState == ERW_notInitialized)
    errorFlag = EC_IllegalCall;
  else
  {
    errorFlag = outStream.status();
    if (errorFlag.good() && fTransferState != ERW_ready)
    {
      if (fTransferState == ERW_init)
      {
        if (outStream.avail() >= 8)
        {
          if (enctype == EET_ExplicitLength)
            Length = getLength(oxfer, enctype);
          else
            Length = DCM_UndefinedLength;
          errorFlag = writeTag(outStream, Tag, oxfer);
          Uint32 valueLength = Length;
          DcmXfer outXfer(oxfer);
          const E_ByteOrder oByteOrder = outXfer.getByteOrder();
          if (oByteOrder == EBO_unknown) return EC_IllegalCall;
          swapIfNecessary(oByteOrder, gLocalByteOrder, &valueLength, 4, 4);
          outStream.write(&valueLength, 4); // 4 bytes length
          elementList->seek(ELP_first);
          fTransferState = ERW_inWork;
        } else
          errorFlag = EC_StreamNotifyClient;
      }
      if (fTransferState == ERW_inWork)
      {
        // elementList->get() can be NULL if buffer was full after
        // writing the last item but before writing the sequence delimitation.
        if (!elementList->empty() && (elementList->get() != NULL))
        {
          DcmObject *dO = NULL;
          do
          {
              dO = elementList->get();
              if (dO->transferState() != ERW_ready)
                errorFlag = dO->write(outStream, oxfer, enctype);
          } while (errorFlag.good() && elementList->seek(ELP_next));
        }
        if (errorFlag.good())
        {
          fTransferState = ERW_ready;
          if (Length == DCM_UndefinedLength)
          {
            if (outStream.avail() >= 8)
            {
                // write Item delimitation
                DcmTag delim(DCM_ItemDelimitationItem);
                errorFlag = writeTag(outStream, delim, oxfer);
                Uint32 delimLen = 0L;
                outStream.write(&delimLen, 4); // 4 bytes length
            }
            else
            {
                // Every subelement of the item is written but it
                // is not possible to write the delimination item into the buffer.
                errorFlag = EC_StreamNotifyClient;
                fTransferState = ERW_inWork;
            }
          }
        }
      }
    }
  }
  return errorFlag;
}

// ********************************

OFCondition DcmItem::writeSignatureFormat(DcmOutputStream &outStream,
                                          const E_TransferSyntax oxfer,
                                          const E_EncodingType enctype)
{
  if (fTransferState == ERW_notInitialized)
    errorFlag = EC_IllegalCall;
  else
  {
    errorFlag = outStream.status();
    if (errorFlag.good() && fTransferState != ERW_ready)
    {
      if (fTransferState == ERW_init)
      {
        if (outStream.avail() >= 4)
        {
          if (enctype == EET_ExplicitLength)
            Length = getLength(oxfer, enctype);
          else
            Length = DCM_UndefinedLength;
          errorFlag = writeTag(outStream, Tag, oxfer);
          /* we don't write the item length */
          elementList->seek(ELP_first);
          fTransferState = ERW_inWork;
        } else
          errorFlag = EC_StreamNotifyClient;
      }
      if (fTransferState == ERW_inWork)
      {
        // elementList->get() can be NULL if buffer was full after
        // writing the last item but before writing the sequence delimitation.
        if (!elementList->empty() && (elementList->get() != NULL))
        {
          DcmObject *dO = NULL;
          do
          {
            dO = elementList->get();
            if (dO->transferState() != ERW_ready)
              errorFlag = dO->writeSignatureFormat(outStream, oxfer, enctype);
          } while (errorFlag.good() && elementList->seek(ELP_next));
        }
        if (errorFlag.good())
        {
          fTransferState = ERW_ready;
          /* we don't write an item delimitation even if the item has undefined length */
        }
      }
    }
  }
  return errorFlag;
}


// ********************************


void DcmItem::transferInit()
{
    DcmObject::transferInit();
    fStartPosition = 0;
    lastElementComplete = OFTrue;
    privateCreatorCache.clear();
    if (!elementList->empty())
    {
        elementList->seek(ELP_first);
        do {
            elementList->get()->transferInit();
        } while (elementList->seek(ELP_next));
    }
}


// ********************************


void DcmItem::transferEnd()
{
    DcmObject::transferEnd();
    privateCreatorCache.clear();
    if (!elementList->empty())
    {
        elementList->seek(ELP_first);
        do {
            elementList->get()->transferEnd();
        } while (elementList->seek(ELP_next));
    }
}


// ********************************


unsigned long DcmItem::card() const
{
    return elementList->card();
}


// ********************************


OFCondition DcmItem::insert(DcmElement *elem,
                            OFBool replaceOld,
                            OFBool checkInsertOrder)
{
    /* intialize error flag with ok */
    errorFlag = EC_Normal;
    /* do something only if the pointer which was passed does not equal NULL */
    if (elem != NULL)
    {
        DcmElement *dE;
        E_ListPos seekmode = ELP_last;
        /* iterate through elementList (from the last element to the first) */
        do {
            /* get current element from elementList */
            dE = OFstatic_cast(DcmElement *, elementList->seek(seekmode));
            /* if there is no element, i.e. elementList is empty */
            if (dE == NULL)
            {
                /* insert new element at the beginning of elementList */
                elementList->insert(elem, ELP_first);
                if (checkInsertOrder)
                {
                  // check if we have inserted at the end of the list
                  if (elem != OFstatic_cast(DcmElement *, elementList->seek(ELP_last)))
                  {
                    // produce diagnostics
                    ofConsole.lockCerr()
                       << "DcmItem: Dataset not in ascending tag order, at element "
                       << elem->getTag() << endl;
                    ofConsole.unlockCerr();
                  }
                }
                /* dump some information if required */
                DCM_dcmdataDebug(3, ("DcmItem::Insert() element (0x%4.4x,0x%4.4x) / VR=\"%s\" at beginning inserted",
                        elem->getGTag(), elem->getETag(), DcmVR(elem->getVR()).getVRName()));
                /* terminate do-while-loop */
                break;
            }
            /* else if the new element's tag is greater than the current element's tag */
            /* (i.e. we have found the position where the new element shall be inserted) */
            else if (elem->getTag() > dE->getTag())
            {
                /* insert the new element after the current element */
                elementList->insert(elem, ELP_next);
                if (checkInsertOrder)
                {
                  // check if we have inserted at the end of the list
                  if (elem != OFstatic_cast(DcmElement *, elementList->seek(ELP_last)))
                  {
                    // produce diagnostics
                    ofConsole.lockCerr()
                       << "DcmItem: Dataset not in ascending tag order, at element "
                       << elem->getTag() << endl;
                    ofConsole.unlockCerr();
                  }
                }
                /* dump some information if required */
                DCM_dcmdataDebug(3, ("DcmItem::Insert() element (0x%4.4x,0x%4.4x) / VR=\"%s\" inserted",
                        elem->getGTag(), elem->getETag(),
                        DcmVR(elem->getVR()).getVRName()));
                /* terminate do-while-loop */
                break;
            }
            /* else if the current element and the new element show the same tag */
            else if (elem->getTag() == dE->getTag())
            {
                /* if new and current element are not identical */
                if (elem != dE)
                {
                    /* if the current (old) element shall be replaced */
                    if (replaceOld)
                    {
                        /* remove current element from list */
                        DcmObject *remObj = elementList->remove();

                        /* now the following holds: remObj == dE and elementList */
                        /* points to the element after the former current element. */

                        /* dump some information if required */
                        DCM_dcmdataDebug(3, ("DcmItem::insert:element (0x%4.4x,0x%4.4x) VR=\"%s\" p=%p removed",
                                remObj->getGTag(), remObj->getETag(),
                                DcmVR(remObj->getVR()).getVRName(), remObj));

                        /* if the pointer to the removed object does not */
                        /* equal NULL (the usual case), delete this object */
                        /* and dump some information if required */
                        if (remObj != NULL)
                        {
                            delete remObj;
                            DCM_dcmdataDebug(3, ("DcmItem::insert:element p=%p deleted", remObj));
                        }
                        /* insert the new element before the current element */
                        elementList->insert(elem, ELP_prev);
                        /* dump some information if required */
                        DCM_dcmdataDebug(3, ("DcmItem::insert() element (0x%4.4x,0x%4.4x) VR=\"%s\" p=%p replaced older one",
                                elem->getGTag(), elem->getETag(),
                                DcmVR(elem->getVR()).getVRName(), elem));

                    }   // if (replaceOld)
                    /* or else, i.e. the current element shall not be replaced by the new element */
                    else {
                        /* set the error flag correspondingly; we do not */
                        /* allow two elements with the same tag in elementList */
                        errorFlag = EC_DoubledTag;
                    }   // if (!replaceOld)
                }   // if (elem != dE)
                /* if the new and the current element are identical, the caller tries to insert */
                /* one element twice. Most probably an application error. */
                else {
                    errorFlag = EC_DoubledTag;
                }
                /* terminate do-while-loop */
                break;
            }
            /* set the seek mode to "get the previous element" */
            seekmode = ELP_prev;
        } while (dE);
    }
    /* if the pointer which was passed equals NULL, this is an illegal call */
    else
        errorFlag = EC_IllegalCall;
    /* return result value */
    return errorFlag;
}


// ********************************


DcmElement *DcmItem::getElement(const unsigned long num)
{
    errorFlag = EC_Normal;
    DcmElement *elem;
    elem = OFstatic_cast(DcmElement *, elementList->seek_to(num));
    // liest Element aus Liste
    if (elem == NULL)
        errorFlag = EC_IllegalCall;
    return elem;
}


// ********************************


DcmObject *DcmItem::nextInContainer(const DcmObject *obj)
{
    if (!obj)
        return elementList->get(ELP_first);
    else
    {
        if (elementList->get() != obj)
        {
            for(DcmObject * search_obj = elementList->seek(ELP_first);
                search_obj && search_obj != obj;
                search_obj = elementList->seek(ELP_next)
               ) {
                /* do nothing, just keep iterating */
            }
        }
        return elementList->seek(ELP_next);
    }
}


// ********************************


OFCondition DcmItem::nextObject(DcmStack &stack,
                                const OFBool intoSub)
{
    OFCondition l_error = EC_Normal;
    DcmObject * container = NULL;
    DcmObject * obj = NULL;
    DcmObject * result = NULL;
    OFBool examSub = intoSub;

    if (stack.empty())
    {
        stack.push(this);
        examSub = OFTrue;
    }

    obj = stack.top();
    if (obj->isLeaf() || !intoSub)
    {
        stack.pop();
        if (stack.card() > 0)
        {
            container = stack.top();
            result = container->nextInContainer(obj);
        }
    } else if (examSub)
        result = obj->nextInContainer(NULL);

    if (result)
        stack.push(result);
    else if (intoSub)
        l_error = nextUp(stack);
    else
        l_error = EC_SequEnd;

    return l_error;
}


// ********************************


DcmElement *DcmItem::remove(const unsigned long num)
{
    errorFlag = EC_Normal;
    DcmElement *elem;
    elem = OFstatic_cast(DcmElement *, elementList->seek_to(num));
    // read element from list
    if (elem != NULL)
        elementList->remove();          // removes element from list but does not delete it
    else
        errorFlag = EC_IllegalCall;
    return elem;
}


// ********************************


DcmElement *DcmItem::remove(DcmObject *elem)
{
    errorFlag = EC_IllegalCall;
    if (!elementList->empty() && elem != NULL)
    {
        DcmObject *dO;
        elementList->seek(ELP_first);
        do {
            dO = elementList->get();
            if (dO == elem)
            {
                elementList->remove();     // removes element from list but does not delete it
                errorFlag = EC_Normal;
                break;
            }
        } while (elementList->seek(ELP_next));
    }
    if (errorFlag == EC_IllegalCall)
        return NULL;
    else
        return OFstatic_cast(DcmElement *, elem);
}


// ********************************


DcmElement *DcmItem::remove(const DcmTagKey &tag)
{
    errorFlag = EC_TagNotFound;
    DcmObject *dO = NULL;
    if (!elementList->empty())
    {
        elementList->seek(ELP_first);
        do {
            dO = elementList->get();
            if (dO->getTag() == tag)
            {
                elementList->remove();     // removes element from list but does not delete it
                errorFlag = EC_Normal;
                break;
            }
        } while (elementList->seek(ELP_next));
    }

    if (errorFlag == EC_TagNotFound)
        return NULL;
    else
        return OFstatic_cast(DcmElement *, dO);
}


// ********************************


OFCondition DcmItem::clear()
{
    errorFlag = EC_Normal;
    DcmObject *dO;
    elementList->seek(ELP_first);
    while (!elementList->empty())
    {
        dO = elementList->remove();
        delete dO;                          // also delete sub elements
    }
    Length = 0;

    return errorFlag;
}


// ********************************


OFCondition DcmItem::verify(const OFBool autocorrect)
{
    DCM_dcmdataDebug(3, ("DcmItem::verify() Tag=(0x%4.4x,0x%4.4x) \"%s\" \"%s\"",
            getGTag(), getETag(), DcmVR(getVR()).getVRName(), Tag.getTagName()));

    errorFlag = EC_Normal;
    if (!elementList->empty())
    {
        DcmObject *dO;
        elementList->seek(ELP_first);
        do {
            dO = elementList->get();
            if (dO->verify(autocorrect).bad())
                errorFlag = EC_CorruptedData;
        } while (elementList->seek(ELP_next));
    }
    if (autocorrect == OFTrue)
        Length = getLength();
    return errorFlag;
}


// ********************************

// Precondition: elementList is non-empty!
// Result:       - return EC_Normal;
//                 push element pointer on resultStack
//               - return EC_TagNotFound;
//                 resultStack unmodified
// Search again: push pointer of sub-element on resultStack and
//               start sub-search

OFCondition DcmItem::searchSubFromHere(const DcmTagKey &tag,
                                       DcmStack &resultStack,
                                       OFBool searchIntoSub)
{
    DcmObject *dO;
    OFCondition l_error = EC_TagNotFound;
    if (!elementList->empty())
    {
        elementList->seek(ELP_first);
        do {
            dO = elementList->get();
            if (searchIntoSub)
            {
                resultStack.push(dO);
                if (dO->getTag() == tag)
                    l_error = EC_Normal;
                else
                    l_error = dO->search(tag, resultStack, ESM_fromStackTop, OFTrue);
                if (l_error.bad())
                    resultStack.pop();
            } else {
                if (dO->getTag() == tag)
                {
                    resultStack.push(dO);
                    l_error = EC_Normal;
                }
            }
        } while (l_error.bad() && elementList->seek(ELP_next));
        DCM_dcmdataCDebug(4, l_error==EC_Normal && dO->getTag()==tag,
               ("DcmItem::searchSubFromHere() Search-Tag=(%4.4x,%4.4x)"
                " found!", tag.getGroup(), tag.getElement()));
    }
    return l_error;
}


// ********************************


OFCondition DcmItem::search(const DcmTagKey &tag,
                            DcmStack &resultStack,
                            E_SearchMode mode,
                            OFBool searchIntoSub)
{
    DcmObject *dO = NULL;
    OFCondition l_error = EC_TagNotFound;
    if (mode == ESM_afterStackTop && resultStack.top() == this)
    {
        l_error = searchSubFromHere(tag, resultStack, searchIntoSub);
    }
    else if (!elementList->empty())
    {
        if (mode == ESM_fromHere || resultStack.empty())
        {
            resultStack.clear();
            l_error = searchSubFromHere(tag, resultStack, searchIntoSub);
        }
        else if (mode == ESM_fromStackTop)
        {
            dO = resultStack.top();
            if (dO == this)
                l_error = searchSubFromHere(tag, resultStack, searchIntoSub);
            else
            {   // gehe direkt zu Sub-Baum und suche dort weiter
                l_error = dO->search(tag, resultStack, mode, searchIntoSub);
// The next two lines destroy the stack->so delete them
//                if (l_error.bad()) // raeumt nur die oberste Stackebene
//                    resultStack.pop();      // ab; der Rest ist unveraendert
            }
        }
        else if (mode == ESM_afterStackTop && searchIntoSub)
        {
            // resultStack enthaelt Zielinformationen:
            // - stelle Zustand der letzen Suche in den einzelnen Suchroutinen
            //   wieder her
            // - finde Position von dO in Baum-Struktur
            //   1. suche eigenen Stack-Eintrag
            //      - bei Fehlschlag Suche beenden
            //   2. nehme naechsthoeheren Eintrag dnO
            //   3. stelle eigene Liste auf Position von dnO
            //   4. starte Suche ab dnO

            unsigned long i = resultStack.card();
            while (i > 0 && (dO = resultStack.elem(i-1)) != this)
            {
                i--;
            }
            if (dO != this && resultStack.card() > 0)
            {                            // oberste Ebene steht nie in resultStack
                i = resultStack.card()+1;// zeige jetzt auf hoechste Ebene+1
                dO = this;               // Treffer der hoechsten Ebene!
            }
            if (dO == this)
            {
                if (i == 1)                   // habe resultStack.top() gefunden
                    l_error = EC_TagNotFound; // markiere als kein Treffer, s.o.
                else                          //   siehe oben
                {
                    E_SearchMode submode = mode;
                    OFBool searchNode = OFTrue;
                    DcmObject *dnO;
                    dnO = resultStack.elem(i - 2); // Knoten der naechsten Ebene
                    elementList->seek(ELP_first);
                    do {
                        dO = elementList->get();
                        searchNode = searchNode ? (dO != dnO) : OFFalse;
                        if (!searchNode)
                        {                             // suche jetzt weiter
                            if (submode == ESM_fromStackTop)
                                resultStack.push(dO); // Stack aktualisieren
                            if (submode == ESM_fromStackTop && dO->getTag() == tag)
                                l_error = EC_Normal;
                            else
                                l_error = dO->search(tag, resultStack, submode, OFTrue);
                            if (l_error.bad())
                                resultStack.pop();
                            else
                                break;
                            submode = ESM_fromStackTop; // ab hier normale Suche
                        }
                    } while (elementList->seek(ELP_next));
                }
            } else
                l_error = EC_IllegalCall;
        } // (mode == ESM_afterStackTop
        else
            l_error = EC_IllegalCall;
    }
    return l_error;
}


// ********************************


OFCondition DcmItem::searchErrors(DcmStack &resultStack)
{
    OFCondition l_error = errorFlag;
    DcmObject *dO = NULL;
    if (errorFlag.bad())
        resultStack.push(this);
    if (!elementList->empty())
    {
        elementList->seek(ELP_first);
        do {
            OFCondition err = EC_Normal;
            dO = elementList->get();
            if ((err = dO->searchErrors(resultStack)).bad())
                l_error = err;
        } while (elementList->seek(ELP_next));
    }
    return l_error;
}


// ********************************


OFCondition DcmItem::loadAllDataIntoMemory()
{
    OFCondition l_error = EC_Normal;
    if (!elementList->empty())
    {
        elementList->seek(ELP_first);
        do {
            OFCondition err = EC_Normal;
            DcmObject *dO = elementList->get();
            if ((err = dO->loadAllDataIntoMemory()).bad())
                l_error = err;
        } while (elementList->seek(ELP_next));
    }
    return l_error;
}




// ********************************

//
// Support functions

OFCondition newDicomElement(DcmElement *&newElement,
                            const DcmTag &tag,
                            const Uint32 length)
{
    DcmTag newTag(tag);
    OFBool readAsUN = OFFalse;
    return newDicomElement(newElement, newTag, length, NULL, readAsUN);
}
                            
DcmElement *newDicomElement(const DcmTag &tag,
                            const Uint32 length)
{
    DcmElement *newElement = NULL;
    newDicomElement(newElement, tag, length);
    return newElement;
}


// ********************************

OFCondition newDicomElement(DcmElement *&newElement,
                            DcmTag &tag,
                            const Uint32 length,
                            DcmPrivateTagCache *privateCreatorCache,
                            OFBool& readAsUN)
{
    /* initialize variables */
    OFCondition l_error = EC_Normal;
    newElement = NULL;
    DcmEVR evr = tag.getEVR();
    readAsUN = OFFalse;
    
    /* revert UN elements with finite length back to known VR if possible */
    if ((evr == EVR_UN) && (length != DCM_UndefinedLength) && dcmEnableUnknownVRConversion.get())
    {
      /* look up VR in data dictionary */
      DcmTag newTag(tag.getGroup(), tag.getElement());

      /* special handling for private elements */
      if (privateCreatorCache && (newTag.getGroup() & 1) && (newTag.getElement() >= 0x1000))
      {
        const char *pc = privateCreatorCache->findPrivateCreator(newTag);
        if (pc)
        {
            // we have a private creator for this element
            newTag.setPrivateCreator(pc);
            newTag.lookupVRinDictionary();
        }
      }

      /* update VR for tag, set "readAsUN" flag that makes sure the attribute value
       * is read in Little Endian Implicit VR (i.e. the UN encoding)
       */
      if (newTag.getEVR() != EVR_UNKNOWN)
      {
        tag.setVR(newTag.getVR());
        evr = tag.getEVR();
        readAsUN = OFTrue;
      }
    }

    /* depending on the VR of the tag which was passed, create the new object */
    switch (evr)
    {
        // byte strings:
        case EVR_AE :
            newElement = new DcmApplicationEntity(tag, length);
            break;
        case EVR_AS :
            newElement = new DcmAgeString(tag, length);
            break;
        case EVR_CS :
            newElement = new DcmCodeString(tag, length);
            break;
        case EVR_DA :
            newElement = new DcmDate(tag, length);
            break;
        case EVR_DS :
            newElement = new DcmDecimalString(tag, length);
            break;
        case EVR_DT :
            newElement = new DcmDateTime(tag, length);
            break;
        case EVR_IS :
            newElement = new DcmIntegerString(tag, length);
            break;
        case EVR_TM :
            newElement = new DcmTime(tag, length);
            break;
        case EVR_UI :
            newElement = new DcmUniqueIdentifier(tag, length);
            break;

        // character strings:
        case EVR_LO :
            newElement = new DcmLongString(tag, length);
            break;
        case EVR_LT :
            newElement = new DcmLongText(tag, length);
            break;
        case EVR_PN :
            newElement = new DcmPersonName(tag, length);
            break;
        case EVR_SH :
            newElement = new DcmShortString(tag, length);
            break;
        case EVR_ST :
            newElement = new DcmShortText(tag, length);
            break;
        case EVR_UT:
            newElement = new DcmUnlimitedText(tag, length);
            break;

        // dependent on byte order:
        case EVR_AT :
            newElement = new DcmAttributeTag(tag, length);
            break;
        case EVR_SS :
            newElement = new DcmSignedShort(tag, length);
            break;
        case EVR_xs : // according to Dicom-Standard V3.0
        case EVR_US :
            newElement = new DcmUnsignedShort(tag, length);
            break;
        case EVR_SL :
            newElement = new DcmSignedLong(tag, length);
            break;
        case EVR_up : // for (0004,eeee) according to Dicom-Standard V3.0
        case EVR_UL :
        {
            // generate Tag with VR from dictionary!
            DcmTag ulupTag(tag.getXTag());
            if (ulupTag.getEVR() == EVR_up)
                newElement = new DcmUnsignedLongOffset(ulupTag, length);
            else
                newElement = new DcmUnsignedLong(tag, length);
        }
        break;
        case EVR_FL:
            newElement = new DcmFloatingPointSingle(tag, length);
            break;
        case EVR_FD :
            newElement = new DcmFloatingPointDouble(tag, length);
            break;
        case EVR_OF:
            newElement = new DcmOtherFloat(tag, length);
            break;

        // sequences and items
        case EVR_SQ :
            newElement = new DcmSequenceOfItems(tag, length);
            break;
        case EVR_na :
            if (tag.getXTag() == DCM_Item)
                l_error = EC_InvalidTag;
            else if (tag.getXTag() == DCM_SequenceDelimitationItem)
                l_error = EC_SequEnd;
            else if (tag.getXTag() == DCM_ItemDelimitationItem)
                l_error = EC_ItemEnd;
            else
                l_error = EC_InvalidTag;
            break;

        // pixel sequences (EVR_pixelSQ) are handled through class DcmPixelData 
        // and should never appear here.

        // unclear 8 or 16 bit:
        case EVR_ox :
            if (tag == DCM_PixelData)
                newElement = new DcmPixelData(tag, length);
            else if (((tag.getGTag() & 0xffe1) == 0x6000)&&(tag.getETag() == 0x3000)) // DCM_OverlayData
                newElement = new DcmOverlayData(tag, length);
            else
                /* we don't know this element's real transfer syntax, so we just
                 * use the defaults of class DcmOtherByteOtherWord and let the
                 * application handle it.
                 */
                newElement = new DcmOtherByteOtherWord(tag, length);
            break;

        case EVR_lt :
            newElement = new DcmOtherByteOtherWord(tag, length);
            break;

        case EVR_OB :
        case EVR_OW :
            if (tag == DCM_PixelData)
                newElement = new DcmPixelData(tag, length);
            else if (((tag.getGTag() & 0xffe1) == 0x6000)&&(tag.getETag() == 0x3000)) // DCM_OverlayData
                newElement = new DcmOverlayData(tag, length);
            else
                if (length == DCM_UndefinedLength) {
                    // The attribute is OB or OW but is encoded with undefined
                    // length.  Assume it is really a sequence so that we can
                    // catch the sequence delimitation item.
                    newElement = new DcmSequenceOfItems(tag, length);
                } else {
                    newElement = new DcmOtherByteOtherWord(tag, length);
                }
            break;

        // read unknown types as byte string:
        case EVR_UNKNOWN :
        case EVR_UNKNOWN2B :
        case EVR_UN :
        default :
            if (length == DCM_UndefinedLength)
            {
              // The attribute VR is UN with undefined length. Assume it is 
              // really a sequence so that we can catch the sequence delimitation item.
              DcmVR sqVR(EVR_SQ); // on writing we will handle this element as SQ, not UN
              DcmTag newTag(tag.getXTag(), sqVR);
              newElement = new DcmSequenceOfItems(newTag, length, dcmEnableCP246Support.get());
            } else {
                // defined length UN element, treat like OB
                newElement = new DcmOtherByteOtherWord(tag, length);
            }
            break;
    }

    /* return result value */
    return l_error;
}


OFCondition nextUp(DcmStack &stack)
{
    DcmObject *oldContainer = stack.pop();
    if (oldContainer->isLeaf())
        return EC_IllegalCall;
    else if (!stack.empty())
    {
        DcmObject *container = stack.top();
        DcmObject *result = container->nextInContainer(oldContainer);
        if (result)
        {
            stack.push(result);
            return EC_Normal;
        }
        else
            return nextUp(stack);
    }
    return EC_TagNotFound;
}


/*
** Simple tests for existance
*/

OFBool DcmItem::tagExists(const DcmTagKey &key,
                          OFBool searchIntoSub)
{
    DcmStack stack;

    OFCondition ec = search(key, stack, ESM_fromHere, searchIntoSub);
    return (ec.good());
}


OFBool DcmItem::tagExistsWithValue(const DcmTagKey &key,
                                   OFBool searchIntoSub)
{
    DcmElement *elem = NULL;
    Uint32 len = 0;
    DcmStack stack;

    OFCondition ec = search(key, stack, ESM_fromHere, searchIntoSub);
    elem = OFstatic_cast(DcmElement *, stack.top());
    if (ec.good() && elem != NULL)
        len = elem->getLength();

    return (ec.good()) && (len > 0);
}


// ********************************

/* --- findAndGet functions: find an element and get it or the value, respectively--- */

OFCondition DcmItem::findAndGetElement(const DcmTagKey &tagKey,
                                       DcmElement *&element,
                                       const OFBool searchIntoSub)
{
    DcmStack stack;
    /* find the element */
    OFCondition status = search(tagKey, stack, ESM_fromHere, searchIntoSub);
    if (status.good())
    {
        element = OFstatic_cast(DcmElement *, stack.top());
        /* should never happen but ... */
        if (element == NULL)
            status = EC_CorruptedData;
    } else {
        /* reset element pointer */
        element = NULL;
    }
    return status;
}


OFCondition DcmItem::findAndGetElements(const DcmTagKey &tagKey,
                                        DcmStack &resultStack)
{
    OFCondition status = EC_TagNotFound;
    DcmStack stack;
    DcmObject *object = NULL;
    /* iterate over all elements */
    while (nextObject(stack, OFTrue).good())
    {
        /* get element */
        object = stack.top();
        if (object->getTag() == tagKey)
        {
            /* add it to the result stack */
            resultStack.push(object);
            status = EC_Normal;
        }
    }
    return status;
}


OFCondition DcmItem::findAndGetString(const DcmTagKey& tagKey,
                                      const char *&value,
                                      const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getString(OFconst_cast(char *&, value));
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetOFString(const DcmTagKey& tagKey,
                                        OFString &value,
                                        const unsigned long pos,
                                        const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getOFString(value, pos);
    }
    /* reset value */
    if (status.bad())
        value.clear();
    return status;
}


OFCondition DcmItem::findAndGetOFStringArray(const DcmTagKey& tagKey,
                                             OFString &value,
                                             const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getOFStringArray(value);
    }
    /* reset value */
    if (status.bad())
        value.clear();
    return status;
}


OFCondition DcmItem::findAndGetUint8(const DcmTagKey& tagKey,
                                     Uint8 &value,
                                     const unsigned long pos,
                                     const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getUint8(value, pos);
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetUint8Array(const DcmTagKey& tagKey,
                                          const Uint8 *&value,
                                          unsigned long *count,
                                          const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        Uint8 *array = NULL;
        status = elem->getUint8Array(array);
        value = array;
    }
    /* set optional count parameter */
    if (count != NULL)
    {
        if (status.good())
            *count = elem->getLength() / sizeof(Uint8);
        else
            *count = 0;
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetUint16(const DcmTagKey& tagKey,
                                      Uint16 &value,
                                      const unsigned long pos,
                                      const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getUint16(value, pos);
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetUint16Array(const DcmTagKey& tagKey,
                                           const Uint16 *&value,
                                           unsigned long *count,
                                           const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        Uint16 *array = NULL;
        status = elem->getUint16Array(array);
        value = array;
    }
    /* set optional count parameter */
    if (count != NULL)
    {
        if (status.good())
            *count = elem->getLength() / sizeof(Uint16);
        else
            *count = 0;
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetSint16(const DcmTagKey& tagKey,
                                      Sint16 &value,
                                      const unsigned long pos,
                                      const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getSint16(value, pos);
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetSint16Array(const DcmTagKey& tagKey,
                                           const Sint16 *&value,
                                           unsigned long *count,
                                           const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        Sint16 *array = NULL;
        status = elem->getSint16Array(array);
        value = array;
    }
    /* set optional count parameter */
    if (count != NULL)
    {
        if (status.good())
            *count = elem->getLength() / sizeof(Sint16);
        else
            *count = 0;
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetUint32(const DcmTagKey& tagKey,
                                      Uint32 &value,
                                      const unsigned long pos,
                                      const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getUint32(value, pos);
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetUint32Array(const DcmTagKey& tagKey,
                                           const Uint32 *&value,
                                           unsigned long *count,
                                           const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        Uint32 *array = NULL;
        status = elem->getUint32Array(array);
        value = array;
    }
    /* set optional count parameter */
    if (count != NULL)
    {
        if (status.good())
            *count = elem->getLength() / sizeof(Uint32);
        else
            *count = 0;
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetSint32(const DcmTagKey& tagKey,
                                      Sint32 &value,
                                      const unsigned long pos,
                                      const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getSint32(value, pos);
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetSint32Array(const DcmTagKey& tagKey,
                                           const Sint32 *&value,
                                           unsigned long *count,
                                           const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        Sint32 *array = NULL;
        status = elem->getSint32Array(array);
        value = array;
    }
    /* set optional count parameter */
    if (count != NULL)
    {
        if (status.good())
            *count = elem->getLength() / sizeof(Sint32);
        else
            *count = 0;
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetLongInt(const DcmTagKey& tagKey,
                                       long int &value,
                                       const unsigned long pos,
                                       const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* distinguish supported VRs */
        switch (elem->ident())
        {
            case EVR_UL:
            case EVR_up:
                Uint32 ul;
                status = elem->getUint32(ul, pos);
                value = OFstatic_cast(long int, ul);
                break;
            case EVR_SL:
            case EVR_IS:
                Sint32 sl;
                status = elem->getSint32(sl, pos);
                value = OFstatic_cast(long int, sl);
                break;
            case EVR_US:
            case EVR_xs:
            case EVR_lt:
                Uint16 us;
                status = elem->getUint16(us, pos);
                value = OFstatic_cast(long int, us);
                break;
            case EVR_SS:
                Sint16 ss;
                status = elem->getSint16(ss, pos);
                value = OFstatic_cast(long int, ss);
                break;
            default:
                status = EC_IllegalCall;
                break;
        }
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetFloat32(const DcmTagKey& tagKey,
                                       Float32 &value,
                                       const unsigned long pos,
                                       const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getFloat32(value, pos);
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetFloat32Array(const DcmTagKey& tagKey,
                                            const Float32 *&value,
                                            unsigned long *count,
                                            const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        Float32 *array = NULL;
        status = elem->getFloat32Array(array);
        value = array;
    }
    /* set optional count parameter */
    if (count != NULL)
    {
        if (status.good())
            *count = elem->getLength() / sizeof(Float32);
        else
            *count = 0;
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetFloat64(const DcmTagKey& tagKey,
                                       Float64 &value,
                                       const unsigned long pos,
                                       const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        status = elem->getFloat64(value, pos);
    }
    /* reset value */
    if (status.bad())
        value = 0;
    return status;
}


OFCondition DcmItem::findAndGetFloat64Array(const DcmTagKey& tagKey,
                                            const Float64 *&value,
                                            unsigned long *count,
                                            const OFBool searchIntoSub)
{
    DcmElement *elem;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* get the value */
        Float64 *array = NULL;
        status = elem->getFloat64Array(array);
        value = array;
    }
    /* set optional count parameter */
    if (count != NULL)
    {
        if (status.good())
            *count = elem->getLength() / sizeof(Float64);
        else
            *count = 0;
    }
    /* reset value */
    if (status.bad())
        value = NULL;
    return status;
}


OFCondition DcmItem::findAndGetSequence(const DcmTagKey &seqTagKey,
                                        DcmSequenceOfItems *&sequence,
                                        const OFBool searchIntoSub)
{
    DcmStack stack;
    /* find the element */
    OFCondition status = search(seqTagKey, stack, ESM_fromHere, searchIntoSub);
    if (status.good())
    {
        DcmElement *delem = OFstatic_cast(DcmElement *, stack.top());
        /* should never happen but ... */
        if (delem == NULL)
            status = EC_CorruptedData;
        /* check for correct VR */
        else if ((delem->ident() == EVR_SQ) || (delem->ident() == EVR_pixelSQ))
            sequence = OFstatic_cast(DcmSequenceOfItems *, delem);
        else
            status = EC_InvalidVR;
    }
    if (status.bad())
    {
        /* reset sequence pointer */
        sequence = NULL;
    }
    return status;
}


OFCondition DcmItem::findAndGetSequenceItem(const DcmTagKey &seqTagKey,
                                            DcmItem *&item,
                                            const signed long itemNum)
{
    DcmStack stack;
    /* find sequence */
    OFCondition status = search(seqTagKey, stack, ESM_fromHere, OFFalse /*searchIntoSub*/);
    if (status.good())
    {
        /* get element */
        DcmElement *delem = OFstatic_cast(DcmElement *, stack.top());
        if (delem != NULL)
        {
            /* check VR */
            if ((delem->ident() == EVR_SQ) || (delem->ident() == EVR_pixelSQ))
            {
                DcmSequenceOfItems *seq = OFstatic_cast(DcmSequenceOfItems *, delem);
                const unsigned long count = seq->card();
                /* empty sequence? */
                if (count > 0)
                {
                    /* get last item */
                    if (itemNum == -1)
                        item = seq->getItem(count - 1);
                    /* get specified item */
                    else if ((itemNum >= 0) && (OFstatic_cast(unsigned long, itemNum) < count))
                        item = seq->getItem(OFstatic_cast(unsigned long, itemNum));
                    /* invalid item number */
                    else
                        status = EC_IllegalParameter;
                } else
                    status = EC_IllegalParameter;
            } else
                status = EC_InvalidVR;
        } else
            status = EC_CorruptedData;
    }
    /* reset item value */
    if (status.bad())
        item = NULL;
    else if (item == NULL)
        status = EC_IllegalCall;
    return status;
}


// ********************************

/* --- findOrCreate functions: find an element or create a new one --- */

OFCondition DcmItem::findOrCreateSequenceItem(const DcmTag& seqTag,
                                              DcmItem *&item,
                                              const signed long itemNum)
{
    DcmStack stack;
    /* find sequence */
    OFCondition status = search(seqTag, stack, ESM_fromHere, OFFalse /*searchIntoSub*/);
    DcmSequenceOfItems *seq = NULL;
    /* sequence found? */
    if (status.good())
    {
        /* get element */
        DcmElement *delem = OFstatic_cast(DcmElement *, stack.top());
        if (delem != NULL)
        {
            /* check VR */
            if ((delem->ident() == EVR_SQ) || (delem->ident() == EVR_pixelSQ))
                seq = OFstatic_cast(DcmSequenceOfItems *, delem);
            else
                status = EC_InvalidVR;
        } else
            status = EC_CorruptedData;
    } else {
        /* create new sequence element */
        seq = new DcmSequenceOfItems(seqTag);
        if (seq != NULL)
        {
            /* insert into item/dataset */
            status = insert(seq, OFTrue /*replaceOld*/);
            if (status.bad())
                delete seq;
        } else
            status = EC_MemoryExhausted;
    }
    if (status.good())
    {
        if (seq != NULL)
        {
            const unsigned long count = seq->card();
            /* existing item? */
            if ((count > 0) && (itemNum >= -1) && (itemNum < OFstatic_cast(signed long, count)))
            {
                if (itemNum == -1)
                {
                    /* get last item */
                    item = seq->getItem(count - 1);
                } else {
                    /* get specified item */
                    item = seq->getItem(OFstatic_cast(unsigned long, itemNum));
                }
            /* create new item(s) */
            } else {
                unsigned long i = 0;
                /* create empty trailing items if required */
                const unsigned long itemCount = (itemNum > OFstatic_cast(signed long, count)) ? (itemNum - count + 1) : 1;
                while ((i < itemCount) && (status.good()))
                {
                    item = new DcmItem();
                    if (item != NULL)
                    {
                        /* append new item to end of sequence */
                        status = seq->append(item);
                        if (status.bad())
                            delete item;
                    } else
                        status = EC_MemoryExhausted;
                    i++;
                }
            }
        } else
            status = EC_IllegalCall;
    }
    /* reset item value */
    if (status.bad())
        item = NULL;
    else if (item == NULL)
        status = EC_IllegalCall;
    return status;
}


// ********************************

/* --- findAndXXX functions: find an element and do something with it --- */

OFCondition DcmItem::findAndDeleteElement(const DcmTagKey &tagKey,
                                          const OFBool allOccurrences,
                                          const OFBool searchIntoSub)
{
    OFCondition status = EC_TagNotFound;
    DcmStack stack;
    DcmObject *object = NULL;
    OFBool intoSub = OFTrue;
    /* iterate over all elements */
    while (nextObject(stack, intoSub).good())
    {
        /* get element */
        object = stack.top();
        if (object->getTag() == tagKey)
        {
            stack.pop();
            /* remove element from dataset and free memory */
            delete OFstatic_cast(DcmItem *, stack.top())->remove(object);
            status = EC_Normal;
            /* delete only the first element? */
            if (!allOccurrences)
                break;
        }
        intoSub = searchIntoSub || allOccurrences;
    }
    return status;
}


OFCondition DcmItem::findAndCopyElement(const DcmTagKey &tagKey,
                                        DcmElement *&newElement,
                                        const OFBool searchIntoSub)
{
    DcmElement *elem = NULL;
    /* find the element */
    OFCondition status = findAndGetElement(tagKey, elem, searchIntoSub);
    if (status.good())
    {
        /* create copy of element */
        newElement = OFstatic_cast(DcmElement *, elem->clone());
        if (newElement == NULL)
            status = EC_MemoryExhausted;
    } else
        newElement = NULL;
    return status;
}


// ********************************

/* --- putAndInsert functions: put value and insert new element --- */

OFCondition DcmItem::putAndInsertString(const DcmTag& tag,
                                        const char *value,
                                        const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    /* create new element */
    DcmElement *elem = NULL;
    switch(tag.getEVR())
    {
        case EVR_AE:
            elem = new DcmApplicationEntity(tag);
            break;
        case EVR_AS:
            elem = new DcmAgeString(tag);
            break;
        case EVR_AT:
            elem = new DcmAttributeTag(tag);
            break;
        case EVR_CS:
            elem = new DcmCodeString(tag);
            break;
        case EVR_DA:
            elem = new DcmDate(tag);
            break;
        case EVR_DS:
            elem = new DcmDecimalString(tag);
            break;
        case EVR_DT:
            elem = new DcmDateTime(tag);
            break;
        case EVR_FL:
            elem = new DcmFloatingPointSingle(tag);
            break;
        case EVR_FD:
            elem = new DcmFloatingPointDouble(tag);
            break;
        case EVR_IS:
            elem = new DcmIntegerString(tag);
            break;
        case EVR_LO:
            elem = new DcmLongString(tag);
            break;
        case EVR_LT:
            elem = new DcmLongText(tag);
            break;
        case EVR_PN:
            elem = new DcmPersonName(tag);
            break;
        case EVR_OB:
        case EVR_OW:
            elem = new DcmOtherByteOtherWord(tag);
            break;
        case EVR_OF:
            elem = new DcmOtherFloat(tag);
            break;
        case EVR_SH:
            elem = new DcmShortString(tag);
            break;
        case EVR_SL:
            elem = new DcmSignedLong(tag);
            break;
        case EVR_SS:
            elem = new DcmSignedShort(tag);
            break;
        case EVR_ST:
            elem = new DcmShortText(tag);
            break;
        case EVR_TM:
            elem = new DcmTime(tag);
            break;
        case EVR_UI:
            elem = new DcmUniqueIdentifier(tag);
            break;
        case EVR_UL:
            elem = new DcmUnsignedLong(tag);
            break;
        case EVR_US:
            elem = new DcmUnsignedShort(tag);
            break;
        case EVR_UT:
            elem = new DcmUnlimitedText(tag);
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putString(value);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else if (status.good())
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertOFStringArray(const DcmTag& tag,
                                               const OFString &value,
                                               const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    /* create new element */
    DcmElement *elem = NULL;
    switch(tag.getEVR())
    {
        case EVR_AE:
            elem = new DcmApplicationEntity(tag);
            break;
        case EVR_AS:
            elem = new DcmAgeString(tag);
            break;
        case EVR_CS:
            elem = new DcmCodeString(tag);
            break;
        case EVR_DA:
            elem = new DcmDate(tag);
            break;
        case EVR_DS:
            elem = new DcmDecimalString(tag);
            break;
        case EVR_DT:
            elem = new DcmDateTime(tag);
            break;
        case EVR_IS:
            elem = new DcmIntegerString(tag);
            break;
        case EVR_LO:
            elem = new DcmLongString(tag);
            break;
        case EVR_LT:
            elem = new DcmLongText(tag);
            break;
        case EVR_PN:
            elem = new DcmPersonName(tag);
            break;
        case EVR_SH:
            elem = new DcmShortString(tag);
            break;
        case EVR_ST:
            elem = new DcmShortText(tag);
            break;
        case EVR_TM:
            elem = new DcmTime(tag);
            break;
        case EVR_UI:
            elem = new DcmUniqueIdentifier(tag);
            break;
        case EVR_UT:
            elem = new DcmUnlimitedText(tag);
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putOFStringArray(value);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else if (status.good())
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertUint8Array(const DcmTag& tag,
                                            const Uint8 *value,
                                            const unsigned long count,
                                            const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    /* create new element */
    DcmElement *elem = NULL;
    switch(tag.getEVR())
    {
        case EVR_OB:
            elem = new DcmOtherByteOtherWord(tag);
            break;
        case EVR_ox:
            /* special handling for Pixel Data */
            if (tag == DCM_PixelData)
            {
                elem = new DcmPixelData(tag);
                if (elem != NULL)
                    elem->setVR(EVR_OB);
            } else
                elem = new DcmPolymorphOBOW(tag);
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putUint8Array(value, count);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else if (status.good())
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertUint16(const DcmTag& tag,
                                        const Uint16 value,
                                        const unsigned long pos,
                                        const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    DcmElement *elem = NULL;
    /* create new element */
    switch(tag.getEVR())
    {
        case EVR_US:
            elem = new DcmUnsignedShort(tag);
            break;
        case EVR_lt:
        case EVR_xs:
            /* special handling */
            elem = new DcmUnsignedShort(DcmTag(tag, EVR_US));
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putUint16(value, pos);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertUint16Array(const DcmTag& tag,
                                             const Uint16 *value,
                                             const unsigned long count,
                                             const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    /* create new element */
    DcmElement *elem = NULL;
    switch(tag.getEVR())
    {
        case EVR_AT:
            elem = new DcmAttributeTag(tag);
            break;
        case EVR_lt:
        case EVR_OW:
            elem = new DcmOtherByteOtherWord(tag);
            break;
        case EVR_US:
            elem = new DcmUnsignedShort(tag);
            break;
        case EVR_ox:
            /* special handling */
            if (tag == DCM_PixelData)
                elem = new DcmPixelData(tag);
            else
                elem = new DcmPolymorphOBOW(tag);
            break;
        case EVR_xs:
            /* special handling */
            elem = new DcmUnsignedShort(DcmTag(tag, EVR_US));
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putUint16Array(value, count);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else if (status.good())
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertSint16(const DcmTag& tag,
                                        const Sint16 value,
                                        const unsigned long pos,
                                        const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    DcmElement *elem = NULL;
    /* create new element */
    switch(tag.getEVR())
    {
        case EVR_SS:
            elem = new DcmSignedShort(tag);
            break;
        case EVR_lt:
        case EVR_xs:
            /* special handling */
            elem = new DcmSignedShort(DcmTag(tag, EVR_SS));
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putSint16(value, pos);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertSint16Array(const DcmTag& tag,
                                             const Sint16 *value,
                                             const unsigned long count,
                                             const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    DcmElement *elem = NULL;
    /* create new element */
    switch(tag.getEVR())
    {
        case EVR_SS:
            elem = new DcmSignedShort(tag);
            break;
        case EVR_lt:
        case EVR_xs:
            /* special handling */
            elem = new DcmSignedShort(DcmTag(tag, EVR_SS));
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putSint16Array(value, count);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertUint32(const DcmTag& tag,
                                        const Uint32 value,
                                        const unsigned long pos,
                                        const OFBool replaceOld)
{
    OFCondition status = EC_IllegalCall;
    /* create new element */
    if (tag.getEVR() == EVR_UL)
    {
        DcmElement *elem = new DcmUnsignedLong(tag);
        if (elem != NULL)
        {
            /* put value */
            status = elem->putUint32(value, pos);
            /* insert into dataset/item */
            if (status.good())
                status = insert(elem, replaceOld);
            /* could not be inserted, therefore, delete it immediately */
            if (status.bad())
                delete elem;
        } else
            status = EC_MemoryExhausted;
    }
    return status;
}


OFCondition DcmItem::putAndInsertSint32(const DcmTag& tag,
                                        const Sint32 value,
                                        const unsigned long pos,
                                        const OFBool replaceOld)
{
    OFCondition status = EC_IllegalCall;
    /* create new element */
    if (tag.getEVR() == EVR_SL)
    {
        DcmElement *elem = new DcmSignedLong(tag);
        if (elem != NULL)
        {
            /* put value */
            status = elem->putSint32(value, pos);
            /* insert into dataset/item */
            if (status.good())
                status = insert(elem, replaceOld);
            /* could not be inserted, therefore, delete it immediately */
            if (status.bad())
                delete elem;
        } else
            status = EC_MemoryExhausted;
    }
    return status;
}


OFCondition DcmItem::putAndInsertFloat32(const DcmTag& tag,
                                         const Float32 value,
                                         const unsigned long pos,
                                         const OFBool replaceOld)
{
    OFCondition status = EC_IllegalCall;
    /* create new element */
    DcmElement *elem = NULL;
    switch(tag.getEVR())
    {
        case EVR_FL:
            elem = new DcmFloatingPointSingle(tag);
            break;
        case EVR_OF:
            elem = new DcmOtherFloat(tag);
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* put value */
        status = elem->putFloat32(value, pos);
        /* insert into dataset/item */
        if (status.good())
            status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else
        status = EC_MemoryExhausted;
    return status;
}


OFCondition DcmItem::putAndInsertFloat64(const DcmTag& tag,
                                         const Float64 value,
                                         const unsigned long pos,
                                         const OFBool replaceOld)
{
    OFCondition status = EC_IllegalCall;
    /* create new element */
    if (tag.getEVR() == EVR_FD)
    {
        DcmElement *elem = new DcmFloatingPointDouble(tag);
        if (elem != NULL)
        {
            /* put value */
            status = elem->putFloat64(value, pos);
            /* insert into dataset/item */
            if (status.good())
                status = insert(elem, replaceOld);
            /* could not be inserted, therefore, delete it immediately */
            if (status.bad())
                delete elem;
        } else
            status = EC_MemoryExhausted;
    }
    return status;
}


OFCondition DcmItem::insertEmptyElement(const DcmTag& tag,
                                        const OFBool replaceOld)
{
    OFCondition status = EC_Normal;
    /* create new element */
    DcmElement *elem = NULL;
    switch(tag.getEVR())
    {
        case EVR_AE:
            elem = new DcmApplicationEntity(tag);
            break;
        case EVR_AS:
            elem = new DcmAgeString(tag);
            break;
        case EVR_AT:
            elem = new DcmAttributeTag(tag);
            break;
        case EVR_CS:
            elem = new DcmCodeString(tag);
            break;
        case EVR_DA:
            elem = new DcmDate(tag);
            break;
        case EVR_DS:
            elem = new DcmDecimalString(tag);
            break;
        case EVR_DT:
            elem = new DcmDateTime(tag);
            break;
        case EVR_FL:
            elem = new DcmFloatingPointSingle(tag);
            break;
        case EVR_FD:
            elem = new DcmFloatingPointDouble(tag);
            break;
        case EVR_OF:
            elem = new DcmOtherFloat(tag);
            break;
        case EVR_IS:
            elem = new DcmIntegerString(tag);
            break;
        case EVR_OB:
        case EVR_OW:
            elem = new DcmOtherByteOtherWord(tag);
            break;
        case EVR_TM:
            elem = new DcmTime(tag);
            break;
        case EVR_UI:
            elem = new DcmUniqueIdentifier(tag);
            break;
        case EVR_LO:
            elem = new DcmLongString(tag);
            break;
        case EVR_LT:
            elem = new DcmLongText(tag);
            break;
        case EVR_UT:
            elem = new DcmUnlimitedText(tag);
            break;
        case EVR_PN:
            elem = new DcmPersonName(tag);
            break;
        case EVR_SH:
            elem = new DcmShortString(tag);
            break;
        case EVR_SQ :
            elem = new DcmSequenceOfItems(tag);
            break;
        case EVR_ST:
            elem = new DcmShortText(tag);
            break;
        default:
            status = EC_IllegalCall;
            break;
    }
    if (elem != NULL)
    {
        /* insert new element into dataset/item */
        status = insert(elem, replaceOld);
        /* could not be inserted, therefore, delete it immediately */
        if (status.bad())
            delete elem;
    } else if (status.good())
        status = EC_MemoryExhausted;
    return status;
}


OFBool DcmItem::containsUnknownVR() const
{
    if (!elementList->empty())
    {
        elementList->seek(ELP_first);
        do {
            if (elementList->get()->containsUnknownVR())
                return OFTrue;
        } while (elementList->seek(ELP_next));
    }
    return OFFalse;
}


/*
** CVS/RCS Log:
** $Log: dcitem.cc,v $
** Revision 1.97  2005/12/08 15:41:16  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.96  2005/11/28 15:53:13  meichel
** Renamed macros in dcdebug.h
**
** Revision 1.95  2005/11/15 18:28:04  meichel
** Added new global flag dcmEnableUnknownVRConversion that enables the automatic
**   re-conversion of defined length UN elements read in an explicit VR transfer
**   syntax, if the real VR is defined in the data dictionary. Default is OFFalse,
**   i.e. to retain the previous behavior.
**
** Revision 1.94  2005/11/15 16:59:25  meichel
** Added new pseudo VR type EVR_lt that is used for LUT Data when read in
**   implicit VR, which may be US, SS or OW. DCMTK always treats EVR_lt like OW.
**
** Revision 1.93  2005/11/07 16:59:26  meichel
** Cleaned up some copy constructors in the DcmObject hierarchy.
**
** Revision 1.92  2005/06/24 10:04:04  joergr
** Added support for internal VR "xs" to putAndInsertXXX() helper methods.
**
** Revision 1.91  2005/05/10 15:27:18  meichel
** Added support for reading UN elements with undefined length according
**   to CP 246. The global flag dcmEnableCP246Support allows to revert to the
**   prior behaviour in which UN elements with undefined length were parsed
**   like a normal explicit VR SQ element.
**
** Revision 1.90  2004/07/01 12:28:27  meichel
** Introduced virtual clone method for DcmObject and derived classes.
**
** Revision 1.89  2004/03/10 10:25:36  joergr
** Translated remaining German comments.
**
** Revision 1.88  2004/02/04 16:02:56  joergr
** Removed pointer declaration from parameter "resultStack" in method
** findAndGetElements(). Removed e-mail addresses from CVS log.
**
** Revision 1.87  2003/10/15 16:55:43  meichel
** Updated error messages for parse errors
**
** Revision 1.86  2003/10/08 10:25:00  joergr
** Added support for AT, OB, OF, OW, SL, SS, UL, US to putAndInsertString().
**
** Revision 1.85  2003/07/16 14:33:43  joergr
** Added new function findAndGetSequence().
** Adapted type casts to new-style typecast operators defined in ofcast.h.
**
** Revision 1.84  2003/06/26 09:17:29  onken
** Added commandline-application dcmodify.
**
** Revision 1.83  2003/06/02 17:25:28  joergr
** Added new helper function DcmItem::findAndCopyElement().
** Fixed bug in findAndDelete() implementation.
** Added explicit support for class DcmPixelData to putAndInsertUintXXArray().
** Changed implementation of findAndGetXXXArray() to avoid problems with MSVC5.
**
** Revision 1.82  2003/05/20 09:17:46  joergr
** Added new helper methods: findAndGetElement(), findAndGetUint32Array(),
** findAndGetSint32Array(), findAndGetFloat64Array(), findAndDeleteElement().
** Enhanced findAndGetSequenceItem() and findOrCreateSequenceItem() by checking
** the return value of ident() - avoids crashes when applied to non-sequence
** elements.
**
** Revision 1.81  2003/03/21 13:08:04  meichel
** Minor code purifications for warnings reported by MSVC in Level 4
**
** Revision 1.80  2002/12/09 09:30:52  wilkens
** Modified/Added doc++ documentation.
**
** Revision 1.79  2002/12/06 12:57:58  joergr
** Enhanced "print()" function by re-working the implementation and replacing
** the boolean "showFullData" parameter by a more general integer flag.
** Made source code formatting more consistent with other modules/files.
** Replaced some German comments by English translations.
**
** Revision 1.78  2002/11/27 12:06:48  meichel
** Adapted module dcmdata to use of new header file ofstdinc.h
**
** Revision 1.77  2002/08/27 16:55:50  meichel
** Initial release of new DICOM I/O stream classes that add support for stream
**   compression (deflated little endian explicit VR transfer syntax)
**
** Revision 1.76  2002/08/02 15:06:33  joergr
** Fixed problems reported by Sun CC 2.0.1.
**
** Revision 1.75  2002/08/02 08:42:33  joergr
** Added optional 'pos' parameter to the putAndInsertXXX() methods.
**
** Revision 1.74  2002/07/23 14:21:33  meichel
** Added support for private tag data dictionaries to dcmdata
**
** Revision 1.73  2002/07/08 16:15:40  meichel
** Unknown undefined length attributes are now converted into SQ instead of UN.
**
** Revision 1.72  2002/07/08 14:44:39  meichel
** Improved dcmdata behaviour when reading odd tag length. Depending on the
**   global boolean flag dcmAcceptOddAttributeLength, the parser now either accepts
**   odd length attributes or implements the old behaviour, i.e. assumes a real
**   length larger by one.
**
** Revision 1.71  2002/07/04 16:35:31  joergr
** Fixed inconsistent formatting of the print() output.
**
** Revision 1.70  2002/06/26 15:49:59  joergr
** Added support for polymorp OB/OW value representation (e.g. pixel data) to
** putAndInsertUint8/16Array() methods.
**
** Revision 1.69  2002/05/29 09:59:37  meichel
** fixed follow-up problem in DcmItem caused by the modifications
**   committed 2002-05-17.
**
** Revision 1.68  2002/05/17 09:58:24  meichel
** fixed bug in DcmItem which caused the parser to fail if the same attribute
**   tag appeared twice within one dataset (which is illegal in DICOM anyway).
**   Added console warning if the attributes read are not in ascending order.
**
** Revision 1.67  2002/04/25 10:15:56  joergr
** Added support for XML output of DICOM objects.
**
** Revision 1.66  2002/04/16 13:43:17  joergr
** Added configurable support for C++ ANSI standard includes (e.g. streams).
**
** Revision 1.65  2002/04/11 12:28:00  joergr
** Enhanced documentation.
**
** Revision 1.64  2001/12/18 11:37:44  joergr
** Added helper method allowing to create and insert empty elements into an
** item/dataset.
**
** Revision 1.63  2001/11/16 15:55:02  meichel
** Adapted digital signature code to final text of supplement 41.
**
** Revision 1.62  2001/11/09 15:53:53  joergr
** Added new helper routines for managing sequences and items.
**
** Revision 1.61  2001/11/01 14:55:39  wilkens
** Added lots of comments.
**
** Revision 1.60  2001/10/10 15:19:51  joergr
** Changed parameter DcmTagKey to DcmTag in DcmItem::putAndInsert... methods
** to support elements which are not in the data dictionary (e.g. private
** extensions).
**
** Revision 1.59  2001/10/02 11:48:01  joergr
** Added functions to get/put 8 bit values/arrays from/to an item/dataset.
**
** Revision 1.58  2001/10/01 15:04:14  joergr
** Introduced new general purpose functions to get/put DICOM element values
** from/to an item/dataset - removed some old and rarely used functions.
** Added "#include <iomanip.h>" to keep gcc 3.0 quiet.
**
** Revision 1.57  2001/09/25 17:19:50  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.56  2001/06/01 15:49:05  meichel
** Updated copyright header
**
** Revision 1.55  2001/05/07 16:08:09  joergr
** Added support for VR=IS to method findIntegerNumber().
**
** Revision 1.54  2001/05/03 08:15:21  meichel
** Fixed bug in dcmdata sequence handling code that could lead to application
**   failure in rare cases during parsing of a correct DICOM dataset.
**
** Revision 1.53  2000/11/07 16:56:20  meichel
** Initial release of dcmsign module for DICOM Digital Signatures
**
** Revision 1.52  2000/04/14 15:55:05  meichel
** Dcmdata library code now consistently uses ofConsole for error output.
**
** Revision 1.51  2000/03/08 16:26:36  meichel
** Updated copyright header.
**
** Revision 1.50  2000/03/03 15:02:09  joergr
** Corrected bug related to padding of file and item size.
**
** Revision 1.49  2000/03/03 14:05:34  meichel
** Implemented library support for redirecting error messages into memory
**   instead of printing them to stdout/stderr for GUI applications.
**
** Revision 1.48  2000/02/29 11:49:29  meichel
** Removed support for VS value representation. This was proposed in CP 101
**   but never became part of the standard.
**
** Revision 1.47  2000/02/23 15:11:55  meichel
** Corrected macro for Borland C++ Builder 4 workaround.
**
** Revision 1.46  2000/02/10 10:52:19  joergr
** Added new feature to dcmdump (enhanced print method of dcmdata): write
** pixel data/item value fields to raw files.
**
** Revision 1.45  2000/02/02 14:32:52  joergr
** Replaced 'delete' statements by 'delete[]' for objects created with 'new[]'.
**
** Revision 1.44  2000/02/01 10:12:07  meichel
** Avoiding to include <stdlib.h> as extern "C" on Borland C++ Builder 4,
**   workaround for bug in compiler header files.
**
** Revision 1.43  1999/03/31 09:25:30  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.42  1999/03/22 15:55:52  meichel
** New handling of unknown (unsupported) VRs encountered when reading explicit
**   VR data. If the VR string consists of uppercase letters, we assume a
**   "future DICOM VR" and decode it expecting an extended length field
**   (4 bytes). Otherwise, we assume an illegal VR string created by some old
**   equipment (i.e.) "??" and decode without extended length (2 bytes).
**
** Revision 1.41  1998/07/15 15:51:59  joergr
** Removed several compiler warnings reported by gcc 2.8.1 with
** additional options, e.g. missing copy constructors and assignment
** operators, initialization of member variables in the body of a
** constructor instead of the member initialization list, hiding of
** methods by use of identical names, uninitialized member variables,
** missing const declaration of char pointers. Replaced tabs by spaces.
**
** Revision 1.40  1998/01/14 15:23:42  hewett
** Added support for the VRs UT (Unlimited Text) and VS (Virtual String).
**
** Revision 1.39  1998/01/14 09:13:53  meichel
** Corrected bug: Overlay Data elements in the groups
**   6002-601f were handled by DcmOtherByteOtherWord
**   instead of the "polymorphous" DcmOverlayData class.
**
** Revision 1.38  1998/01/14 08:42:32  meichel
** Improved algorithm for auto-detection of transfer syntax
**   used when opening a DICOM file without metaheader.
**   Big endian datasets are now detected much more reliably.
**
** Revision 1.37  1997/11/07 08:52:18  meichel
** Corrected bug in the dcmdata read routines which caused incorrect reading
**   of datasets containing attributes with value representation "ox" (= OB or OW)
**   in the dicom dictionary other than PixelData and OverlayData.
**
** Revision 1.36  1997/09/22 14:50:53  hewett
** - Added 2 simple methods to test for the existance of an attribute
**   to DcmItem class (tagExists and tagExistsWithValue).  This code
**   was part of dcmgpdir.cc but is more generally useful.
** - Added 2 methods to find an attribute and retrieve numeric values
**   to DcmItem class (findIntegerNumber and findRealNumber).  The old
**   method findLong is now marked as obsolete and reimplemented using
**   findIntegerNumber.
**
** Revision 1.35  1997/09/12 13:44:53  meichel
** The algorithm introduced on 97.08.28 to detect incorrect odd-length
**   value fields falsely reported undefined length sequences and items
**   to be wrong. Corrected.
**
** Revision 1.34  1997/08/29 08:31:33  andreas
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
** Revision 1.33  1997/08/29 07:52:40  andreas
** - New error messages if length of an element is odd. Previously, no
**   error was reported. But the length is corrected by the method
**   newValueFiel and. so it was impossible for a checking utility to find
**   such an error in DICOM objects.
**
** Revision 1.32  1997/07/24 13:10:52  andreas
** - Removed Warnings from SUN CC 2.0.1
**
** Revision 1.31  1997/07/21 08:11:42  andreas
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
** Revision 1.30  1997/07/07 07:43:59  andreas
** - Changed type for Tag attribute in DcmObject from prointer to value
** - Changed parameter type DcmTag & to DcmTagKey & in all search functions
**   in DcmItem, DcmSequenceOfItems, DcmDirectoryRecord and DcmObject
** - Enhanced (faster) byte swapping routine. swapIfNecessary moved from
**   a method in DcmObject to a general function.
**
** Revision 1.29  1997/07/03 15:09:59  andreas
** - removed debugging functions Bdebug() and Edebug() since
**   they write a static array and are not very useful at all.
**   Cdebug and Vdebug are merged since they have the same semantics.
**   The debugging functions in dcmdata changed their interfaces
**   (see dcmdata/include/dcdebug.h)
**
** Revision 1.28  1997/06/06 09:55:29  andreas
** - corrected error: canWriteXfer returns false if the old transfer syntax
**   was unknown, which causes several applications to prohibit the writing
**   of dataset.
**
** Revision 1.27  1997/05/27 13:49:00  andreas
** - Add method canWriteXfer to class DcmObject and all derived classes.
**   This method checks whether it is possible to convert the original
**   transfer syntax to an new transfer syntax. The check is used in the
**   dcmconv utility to prohibit the change of a compressed transfer
**   syntax to a uncompressed.
**
** Revision 1.26  1997/05/16 08:13:49  andreas
** - Revised handling of GroupLength elements and support of
**   DataSetTrailingPadding elements. The enumeratio E_GrpLenEncoding
**   got additional enumeration values (for a description see dctypes.h).
**   addGroupLength and removeGroupLength methods are replaced by
**   computeGroupLengthAndPadding. To support Padding, the parameters of
**   element and sequence write functions changed.
** - Added a new method calcElementLength to calculate the length of an
**   element, item or sequence. For elements it returns the length of
**   tag, length field, vr field, and value length, for item and
**   sequences it returns the length of the whole item. sequence including
**   the Delimitation tag (if appropriate).  It can never return
**   UndefinedLength.
** - Deleted obsolete method DcmItem::calcHeaderLength because the
**   samce functionality is performed by DcmXfer::sizeofTagHeader
**
** Revision 1.25  1997/05/07 12:27:27  andreas
** Corrected error reading ItemDelimitationItem using explicit transfer syntaxes
**
** Revision 1.24  1997/04/30 16:32:50  andreas
** - Corrected Bug for reading of encapsulated pixel sequences
**
** Revision 1.23  1997/04/24 12:12:18  hewett
** Fixed DICOMDIR generation bug affecting the use of Unknown VR
** attributes (the file offsets were not being computed correctly).
**
** Revision 1.22  1997/04/18 08:10:49  andreas
** - Corrected debugging code
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
** Revision 1.21  1997/03/27 15:52:50  hewett
** Extended preliminary support for Unknown VR (UN) described in
** Supplement 14.  Attributes with undefined length should now
** be handled as a sequence.
**
** Revision 1.20  1997/03/26 17:15:57  hewett
** Added very preliminary support for Unknown VR (UN) described in
** Supplement 14.  WARNING: handling of unknown attributes with undefined
** length is not yet supported.
**
** Revision 1.19  1996/09/24 15:54:14  hewett
** Corrected erroneous setting of an error flag when inserting an
** attribute into an Item (via Item::insert(...)) and the attribute
** was already present.  Now the error flag is only set if replaceOld
** is OFFalse and an attribute already exists.
**
** Revision 1.18  1996/09/13 12:04:12  hewett
** Corrected missing () in function call (stack.card()) used in nextObject(...)
**
** Revision 1.17  1996/08/08 10:15:09  andreas
** Some more testing in nextObject
**
** Revision 1.16  1996/08/08 10:06:23  andreas
** Correct error for intoSub=OFFalse
**
** Revision 1.15  1996/08/05 08:46:12  andreas
** new print routine with additional parameters:
**         - print into files
**         - fix output length for elements
** corrected error in search routine with parameter ESM_fromStackTop
**
** Revision 1.14  1996/07/31 13:14:30  andreas
** - Minor corrections: error code for swapping to or from byteorder unknown
**                      correct read of dataset in fileformat
**
** Revision 1.13  1996/07/17 12:39:38  andreas
** new nextObject for DcmDataSet, DcmFileFormat, DcmItem, ...
**
** Revision 1.12  1996/04/29 15:08:14  hewett
** Replaced DcmItem::findInt(...) with the more general DcmItem::findLong(...).
**
** Revision 1.11  1996/04/27 14:04:55  hewett
** Eliminated compiler warnings when compiling without -DDEBUG.  Very
** minor corrections, mostly unused parameters and uninitialized variables.
**
** Revision 1.10  1996/04/16 16:04:54  andreas
** - const tag Parameter in newDicomElement
**
** Revision 1.9  1996/03/28 18:52:39  hewett
** Added 2 simple find&get methods (findString & findInt).
**
** Revision 1.8  1996/03/12 15:23:27  hewett
** When generating group length tags, the VR of a tag is now explicity
** set to be EVR_UL.  Group length tags not in the dictionary (e.g. for
** private groups) were getting coded incorrectly.
**
** Revision 1.7  1996/03/11 14:16:00  hewett
** Corrected error whereby explicit encoding was being recognised as implicit.
**
** Revision 1.6  1996/03/11 13:03:51  hewett
** Rearranged logic of DcmItem::checkTransferSyntax to make little-endian
** the default if both big and little endian are possible.
**
** Revision 1.5  1996/01/29 13:38:27  andreas
** - new put method for every VR to put value as a string
** - better and unique print methods
**
** Revision 1.4  1996/01/09 11:06:46  andreas
** New Support for Visual C++
** Correct problems with inconsistent const declarations
** Correct error in reading Item Delimitation Elements
**
** Revision 1.3  1996/01/05 13:27:37  andreas
** - changed to support new streaming facilities
** - unique read/write methods for file and block transfer
** - more cleanups
**
*/
