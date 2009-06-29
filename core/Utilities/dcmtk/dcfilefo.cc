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
 *  Purpose: class DcmFileFormat
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:41:10 $
 *  CVS/RCS Revision: $Revision: 1.42 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */

#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CTIME
#define INCLUDE_UNISTD
#include "dcmtk/ofstd/ofstdinc.h"

#include "dcmtk/ofstd/ofstream.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcitem.h"
#include "dcmtk/dcmdata/dcxfer.h"
#include "dcmtk/dcmdata/dcvrobow.h"
#include "dcmtk/dcmdata/dcvrui.h"
#include "dcmtk/dcmdata/dcvrul.h"
#include "dcmtk/dcmdata/dcvrus.h"
#include "dcmtk/dcmdata/dcvrae.h"
#include "dcmtk/dcmdata/dcvrsh.h"
#include "dcmtk/dcmdata/dcmetinf.h"
#include "dcmtk/dcmdata/dcdebug.h"

#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmdata/dcuid.h"
#include "dcmtk/dcmdata/dcostrma.h"    /* for class DcmOutputStream */
#include "dcmtk/dcmdata/dcostrmf.h"    /* for class DcmOutputFileStream */
#include "dcmtk/dcmdata/dcistrma.h"    /* for class DcmInputStream */
#include "dcmtk/dcmdata/dcistrmf.h"    /* for class DcmInputFileStream */


// ********************************


DcmFileFormat::DcmFileFormat()
  : DcmSequenceOfItems(InternalUseTag),
    FileReadMode(ERM_autoDetect)
{
    DcmMetaInfo *MetaInfo = new DcmMetaInfo();
    DcmSequenceOfItems::itemList->insert(MetaInfo);
    DcmDataset *Dataset = new DcmDataset();
    DcmSequenceOfItems::itemList->insert(Dataset);
}


DcmFileFormat::DcmFileFormat(DcmDataset *dataset)
    : DcmSequenceOfItems(InternalUseTag),
      FileReadMode(ERM_autoDetect)
{
    DcmMetaInfo *MetaInfo = new DcmMetaInfo();
    DcmSequenceOfItems::itemList->insert(MetaInfo);
    DcmDataset *newDataset;

    if (dataset == NULL)
        newDataset = new DcmDataset();
    else
        newDataset = new DcmDataset(*dataset);
    DcmSequenceOfItems::itemList->insert(newDataset);
}


DcmFileFormat::DcmFileFormat(const DcmFileFormat &old)
  : DcmSequenceOfItems(old),
    FileReadMode(old.FileReadMode)
{
}


DcmFileFormat::~DcmFileFormat()
{
}


DcmFileFormat &DcmFileFormat::operator=(const DcmFileFormat &obj)
{
    DcmSequenceOfItems::operator=(obj);
    FileReadMode = obj.FileReadMode;
    return *this;
}


// ********************************


DcmEVR DcmFileFormat::ident() const
{
    return EVR_fileFormat;
}


// ********************************


void DcmFileFormat::print(ostream &out,
                          const size_t flags,
                          const int level,
                          const char *pixelFileName,
                          size_t *pixelCounter)
{
    out << endl;
    printNestingLevel(out, flags, level);
    out << "# Dicom-File-Format" << endl;
    if (!itemList->empty())
    {
        DcmObject *dO;
        itemList->seek(ELP_first);
        do {
            dO = itemList->get();
            dO->print(out, flags, level, pixelFileName, pixelCounter);
        } while (itemList->seek(ELP_next));
    } else {
        printNestingLevel(out, flags, level);
        out << "# Dicom-File-Format has been erased";
        out << endl;
    }
}


// ********************************


OFCondition DcmFileFormat::writeXML(ostream &out,
                                    const size_t flags)
{
    OFCondition result = EC_CorruptedData;
    /* XML start tag for "file-format" */
    out << "<file-format";
    if (flags & DCMTypes::XF_useDcmtkNamespace)
        out << " xmlns=\"" << DCMTK_XML_NAMESPACE_URI << "\"";
    out << ">" << endl;
    if (!itemList->empty())
    {
        /* write content of all children */
        DcmObject *dO;
        itemList->seek(ELP_first);
        do {
            dO = itemList->get();
            dO->writeXML(out, flags & ~DCMTypes::XF_useDcmtkNamespace);
        } while (itemList->seek(ELP_next));
        result = EC_Normal;
    }
    /* XML end tag for "file-format" */
    out << "</file-format>" << endl;
    return result;
}


// ********************************


OFCondition DcmFileFormat::checkValue(DcmMetaInfo *metainfo,
                                      DcmDataset *dataset,
                                      const DcmTagKey &atagkey,
                                      DcmObject *obj,
                                      const E_TransferSyntax oxfer)
    /*
     * This function checks if a particular data element of the meta header information is existent
     * in metainfo. If the element is not existent, it will be inserted. Additionally, this function
     * makes sure that the corresponding data element will contain a correct value.
     *
     * Parameters:
     *   metainfo - [in] The meta header information.
     *   dataset  - [in] The data set information.
     *   atagkey  - [in] Tag of the data element which shall be checked.
     *   obj      - [in] Data object from metainfo which represents the data element that shall be checked.
     *                   Equals NULL, if this data element is not existent in the meta header information.
     *   oxfer    - [in] The transfer syntax which shall be used.
     */
{
    /* initialize result value */
    OFCondition l_error = EC_Normal;

    /* if there is meta header information and also data set information, do something */
    if (metainfo != NULL && dataset != NULL)
    {
        /* intitialize variables */
        DcmStack stack;
        DcmTag tag(atagkey);
        if (obj != NULL)
            tag = obj->getTag();

        DcmTagKey xtag = tag.getXTag();
        DcmElement *elem = OFstatic_cast(DcmElement *, obj);

        /* go ahaed and scrutinize one particular data element (depending on xtag) */
        if (xtag == DCM_MetaElementGroupLength)             // (0002,0000)
        {
            if (elem == NULL)
            {
                elem = new DcmUnsignedLong(tag);
                metainfo->insert(elem, OFTrue);
            }
            Uint32 temp = 0;
            if (elem->getLength() == 0 && elem->ident() == EVR_UL)
                OFstatic_cast(DcmUnsignedLong *, elem)->putUint32Array(&temp, 1);
            // the calculation of actual group length value is contained in validateMetaInfo()
        }
        else if (xtag == DCM_FileMetaInformationVersion)    // (0002,0001)
        {
            if (elem == NULL)
            {
                elem = new DcmOtherByteOtherWord(tag);
                metainfo->insert(elem, OFTrue);
            }
            Uint8 version[2] = {0,1};
            if (elem->getLength() == 0 && elem->ident() == EVR_OB)
                OFstatic_cast(DcmOtherByteOtherWord *, elem)->putUint8Array(version, 2);

            // check version of meta header
            Uint8 *currVers;
            l_error = OFstatic_cast(DcmOtherByteOtherWord *, elem)->getUint8Array(currVers);
            if (((currVers[0] & version[0] & 0xff) == version[0]) &&
                ((currVers[1] & version[1] & 0xff) == version[1]))
            {
                DCM_dcmdataDebug(2, ("DcmFileFormat::checkValue() Version of MetaHeader is ok: 0x%2.2x%2.2x",
                        currVers[1], currVers[0]));
            } else {
                currVers[0] = OFstatic_cast(Uint8, currVers[0] | version[0]); // direct manipulation
                currVers[1] = OFstatic_cast(Uint8, currVers[1] | version[1]); // of data
                ostream& localCerr = ofConsole.lockCerr();
                localCerr << "Warning: dcfilefo: unknown Version of MetaHeader detected: 0x";
                localCerr << hex << setfill('0')
                          << setw(2) << OFstatic_cast(int, currVers[1])
                          << setw(2) << OFstatic_cast(int, currVers[0]) << " supported: 0x"
                          << setw(2) << OFstatic_cast(int, version[1])
                          << setw(2) << OFstatic_cast(int, version[0])
                          << dec << setfill(' ') << endl;
                ofConsole.unlockCerr();
            }
        }
        else if (xtag == DCM_MediaStorageSOPClassUID)       // (0002,0002)
        {
            if (elem == NULL)
            {
                elem = new DcmUniqueIdentifier(tag);
                metainfo->insert(elem, OFTrue);
            }
            if (elem->getLength() == 0 && elem->ident() == EVR_UI) {
                if (dataset->search(DCM_SOPClassUID, stack).good())
                {
                    char *uid = NULL;
                    l_error = OFstatic_cast(DcmUniqueIdentifier *, stack.top())->getString(uid);
                    OFstatic_cast(DcmUniqueIdentifier *, elem)->putString(uid);
                    DCM_dcmdataDebug(2, ("DcmFileFormat::checkValue() use SOPClassUID [%s]", uid));

                } else {
                    OFstatic_cast(DcmUniqueIdentifier *, elem)->putString(UID_PrivateGenericFileSOPClass);
                    DCM_dcmdataDebug(2, ("DcmFileFormat::checkValue() No SOP Class UID in Dataset, using PrivateGenericFileSOPClass"));
                }
            }
        }
        else if (xtag == DCM_MediaStorageSOPInstanceUID)    // (0002,0003)
        {
            if (elem == NULL)
            {
                elem = new DcmUniqueIdentifier(tag);
                metainfo->insert(elem, OFTrue);
            }
            if (elem->getLength() == 0 && elem->ident() == EVR_UI)
            {
                if (dataset->search(DCM_SOPInstanceUID, stack).good())
                {
                    char* uid = NULL;
                    l_error =OFstatic_cast(DcmUniqueIdentifier *, stack.top())->getString(uid);
                    OFstatic_cast(DcmUniqueIdentifier *, elem)->putString(uid);
                    DCM_dcmdataDebug(2, ("DcmFileFormat::checkValue() use SOPInstanceUID [%s] from Dataset", uid));
                } else {
                    char uid[128];
                    dcmGenerateUniqueIdentifier(uid);       // from dcuid.h
                    OFstatic_cast(DcmUniqueIdentifier *, elem)->putString(uid);
                    DCM_dcmdataDebug(2, ("DcmFileFormat::checkValue() use new generated SOPInstanceUID [%s]", uid));
                }
            }
        }
        else if (xtag == DCM_TransferSyntaxUID)             // (0002,0010)
        {
            if (elem == NULL)
            {
                elem = new DcmUniqueIdentifier(tag);
                metainfo->insert(elem, OFTrue);
            }
            if (elem->ident() == EVR_UI)
            {
#ifdef DEBUG
char * uidtmp = NULL;
OFstatic_cast(DcmUniqueIdentifier *, elem)->getString(uidtmp);
DCM_dcmdataCDebug(2,  uidtmp != NULL,
       ("DcmFileFormat::checkValue() found old transfer-syntax: [%s]",uidtmp));
#endif
                DcmXfer dcXfer(oxfer);
                const char *uid = dcXfer.getXferID();
                elem->putString(uid);
                DCM_dcmdataDebug(2,("DcmFileFormat::checkValue() use new transfer-syntax [%s] on writing following Dataset",
                        dcXfer.getXferName()));

            }
        }
        else if (xtag == DCM_ImplementationClassUID)        // (0002,0012)
        {
            if (elem == NULL)
            {
                elem = new DcmUniqueIdentifier(tag);
                metainfo->insert(elem, OFTrue);
            }
            if (elem->ident() == EVR_UI)
            {
                const char *uid = OFFIS_IMPLEMENTATION_CLASS_UID;
                OFstatic_cast(DcmUniqueIdentifier *, elem)->putString(uid);
            }

        }
        else if (xtag == DCM_ImplementationVersionName)     // (0002,0013)
        {
            if (elem == NULL)
            {
                elem = new DcmShortString(tag);
                metainfo->insert(elem, OFTrue);
            }
            if (elem->ident() == EVR_SH)
            {
                const char *uid = OFFIS_DTK_IMPLEMENTATION_VERSION_NAME;
                OFstatic_cast(DcmShortString *, elem)->putString(uid);
            }

        }
        else if (xtag == DCM_SourceApplicationEntityTitle)  // (0002,0016)
        {
            if (elem == NULL)
            {
                elem = new DcmApplicationEntity(tag);
                metainfo->insert(elem, OFTrue);
            }
            ofConsole.lockCerr() << "Error: dcfilefo: I don't know how to handle DCM_SourceApplicationEntityTitle!" << endl;
            ofConsole.unlockCerr();
        }
        else if (xtag == DCM_PrivateInformationCreatorUID)  // (0002,0100)
        {
            if (elem == NULL)
            {
                elem = new DcmUniqueIdentifier(tag);
                metainfo->insert(elem, OFTrue);
            }
            ofConsole.lockCerr() << "Error: dcfilefo: I don't know how to handle DCM_PrivateInformationCreatorUID!" << endl;
            ofConsole.unlockCerr();
        }
        else if (xtag == DCM_PrivateInformation)            // (0002,0102)
        {
            if (elem == NULL)
            {
                elem = new DcmOtherByteOtherWord(tag);
                metainfo->insert(elem, OFTrue);
            }
            ofConsole.lockCerr() << "Warning: dcfilefo: I don't know how to handle DCM_PrivateInformation!" << endl;
            ofConsole.unlockCerr();
        } else {
            ofConsole.lockCerr() << "Warning: dcfilefo: I don't know how to handle " << tag.getTagName() << endl;
            ofConsole.unlockCerr();
        }

        /* if at this point elem still equals NULL, something is fishy */
        if (elem == NULL)
            l_error = EC_InvalidVR;
    } else {
    /* (i.e. there is either no meta header information or no data set information */
        l_error = EC_IllegalCall;
    }

    /* return result value */
    return l_error;
}


// ********************************


OFCondition DcmFileFormat::validateMetaInfo(E_TransferSyntax oxfer)
    /*
     * This function makes sure that all data elements of the meta header information are existent
     * in metainfo and contain correct values.
     *
     * Parameters:
     *   oxfer          - [in] The transfer syntax which shall be used.
     */
{
    /* initialize some variables */
    OFCondition l_error = EC_Normal;
    DcmMetaInfo *metinf = getMetaInfo();
    DcmDataset *datset = getDataset();

    /* if there is meta header information and data set information, do something */
    if (metinf != NULL && datset != NULL)
    {
        DcmStack stack;

        /* in the following, we want to make sure all elements of the meta header */
        /* are existent in metinf and contain correct values */

        /* DCM_MetaElementGroupLength */
        metinf->search(DCM_MetaElementGroupLength, stack, ESM_fromHere, OFFalse);
        checkValue(metinf, datset, DCM_MetaElementGroupLength, stack.top(), oxfer);

        /* DCM_FileMetaInformationVersion */
        metinf->search(DCM_FileMetaInformationVersion, stack, ESM_fromHere, OFFalse);
        checkValue(metinf, datset, DCM_FileMetaInformationVersion, stack.top(), oxfer);

        /* DCM_MediaStorageSOPClassUID */
        metinf->search(DCM_MediaStorageSOPClassUID, stack, ESM_fromHere, OFFalse);
        checkValue(metinf, datset, DCM_MediaStorageSOPClassUID, stack.top(), oxfer);

        /* DCM_MediaStorageSOPInstanceUID */
        metinf->search(DCM_MediaStorageSOPInstanceUID, stack, ESM_fromHere, OFFalse);
        checkValue(metinf, datset, DCM_MediaStorageSOPInstanceUID, stack.top(), oxfer);

        /* DCM_TransferSyntaxUID */
        metinf->search(DCM_TransferSyntaxUID, stack, ESM_fromHere, OFFalse);
        checkValue(metinf, datset, DCM_TransferSyntaxUID, stack.top(), oxfer);

        /* DCM_ImplementationClassUID */
        metinf->search(DCM_ImplementationClassUID, stack, ESM_fromHere, OFFalse);
        checkValue(metinf, datset, DCM_ImplementationClassUID, stack.top(), oxfer);

        /* DCM_ImplementationVersionName */
        metinf->search(DCM_ImplementationVersionName, stack, ESM_fromHere, OFFalse);
        checkValue(metinf, datset, DCM_ImplementationVersionName, stack.top(), oxfer);

        /* dump some information if reuqired */
        DCM_dcmdataDebug(2, ("DcmFileFormat: found %ld Elements in DcmMetaInfo metinf.",
                metinf->card()));

        /* calculate new GroupLength for meta header */
        if (metinf->computeGroupLengthAndPadding(EGL_withGL, EPD_noChange,
            META_HEADER_DEFAULT_TRANSFERSYNTAX, EET_UndefinedLength).bad())
        {
            ofConsole.lockCerr() << "Error: DcmFileFormat::validateMetaInfo(): group length of Meta Information Header not adapted." << endl;
            ofConsole.unlockCerr();
        }
    } else {
        /* (i.e. there is either no meta header information or no data set information, or both are missing) */
        l_error = EC_CorruptedData;
    }

    /* return result value */
    return l_error;
}



// ********************************


E_TransferSyntax DcmFileFormat::lookForXfer(DcmMetaInfo *metainfo)
{
    E_TransferSyntax newxfer = EXS_Unknown;
    DcmStack stack;
    if (metainfo && metainfo->search(DCM_TransferSyntaxUID, stack).good())
    {
        DcmUniqueIdentifier *xferUI = OFstatic_cast(DcmUniqueIdentifier *, stack.top());
        if (xferUI->getTag().getXTag() == DCM_TransferSyntaxUID)
        {
            char * xferid = NULL;
            xferUI->getString(xferid);     // auslesen der ID
            DcmXfer localXfer(xferid);      // dekodieren in E_TransferSyntax
            newxfer = localXfer.getXfer();
            DCM_dcmdataDebug(4, ("DcmFileFormat::lookForXfer() detected xfer=%d=[%s] in MetaInfo",
                newxfer, localXfer.getXferName()));
        }
    }
    return newxfer;
}


// ********************************


Uint32 DcmFileFormat::calcElementLength(const E_TransferSyntax xfer,
                                        const E_EncodingType enctype)
{
    return getMetaInfo()->calcElementLength(xfer, enctype) +
        getDataset()->calcElementLength(xfer, enctype);
}


// ********************************


OFBool DcmFileFormat::canWriteXfer(const E_TransferSyntax newXfer,
                                   const E_TransferSyntax oldXfer)
{
    DcmDataset *dataset = getDataset();

    if (dataset)
        return dataset->canWriteXfer(newXfer, oldXfer);
    else
        return OFFalse;
}


// ********************************


OFCondition DcmFileFormat::read(DcmInputStream &inStream,
                                const E_TransferSyntax xfer,
                                const E_GrpLenEncoding glenc,
                                const Uint32 maxReadLength)

{
    if (fTransferState == ERW_notInitialized)
        errorFlag = EC_IllegalCall;
    else
    {
        errorFlag = inStream.status();

        E_TransferSyntax newxfer = xfer;
        DcmDataset *dataset = NULL;

        if (errorFlag.good() && inStream.eos())
            errorFlag = EC_EndOfStream;
        else if (errorFlag.good() && fTransferState != ERW_ready)
        {
            // the new data is added to the end
            itemList->seek(ELP_last);

            DcmMetaInfo *metaInfo = getMetaInfo();
            if (metaInfo == NULL && fTransferState == ERW_init)
            {
                metaInfo = new DcmMetaInfo();
                itemList->insert(metaInfo, ELP_first);
            }
            if (metaInfo && metaInfo->transferState() != ERW_ready)
            {
                errorFlag = metaInfo->read(inStream, xfer, glenc, maxReadLength);
            }

            // read MetaInfo() from Tag(0002,0010) and determine xfer
            newxfer = lookForXfer(metaInfo);
            if (FileReadMode == ERM_fileOnly)
            {
                // reject file if no meta header present
                if (errorFlag.good() && (newxfer == EXS_Unknown))
                    errorFlag = EC_InvalidStream;
            }
            if (errorFlag.good() && (!metaInfo || metaInfo->transferState() == ERW_ready))
            {
                dataset = getDataset();
                if (dataset == NULL && fTransferState == ERW_init)
                {
                    dataset = new DcmDataset();
                    itemList->seek (ELP_first);
                    itemList->insert(dataset, ELP_next);
                }
                if (dataset && dataset->transferState() != ERW_ready)
                {
                    errorFlag = dataset->read(inStream, newxfer, glenc, maxReadLength);
                }
            }
        }
        if (fTransferState == ERW_init)
            fTransferState = ERW_inWork;

        if (dataset && dataset->transferState() == ERW_ready)
            fTransferState = ERW_ready;
    }
    return errorFlag;
}  // DcmFileFormat::read()


// ********************************


OFCondition DcmFileFormat::write(DcmOutputStream &outStream,
                                 const E_TransferSyntax oxfer,
                                 const E_EncodingType enctype)
{
    return write(outStream, oxfer, enctype, EGL_recalcGL, EPD_noChange);
}


OFCondition DcmFileFormat::write(DcmOutputStream &outStream,
                                 const E_TransferSyntax oxfer,
                                 const E_EncodingType enctype,
                                 const E_GrpLenEncoding glenc,
                                 const E_PaddingEncoding padenc,
                                 const Uint32 padlen,
                                 const Uint32 subPadlen,
                                 Uint32 instanceLength)
    /*
     * This function writes data values which are contained in this to the stream which is
     * passed as first argument. With regard to the writing of information, the other parameters
     * which are passed are accounted for. The function will return EC_Normal, if the information
     * from all elements of this data set has been written to the buffer, and it will return some
     * other (error) value if there was an error.
     *
     * Parameters:
     *   outStream      - [inout] The stream that the information will be written to.
     *   oxfer          - [in] The transfer syntax which shall be used.
     *   enctype        - [in] Encoding type for sequences. Specifies how sequences will be handled.
     *   glenc          - [in] Encoding type for group length. Specifies what will be done with group length tags.
     *   padenc         - [in] Encoding type for padding. Specifies what will be done with padding tags.
     *   padlen         - [in] The length up to which the dataset shall be padded, if padding is desired.
     *   subPadlen      - [in] For sequences (i.e. sub elements), the length up to which item shall be padded,
     *                         if padding is desired.
     *   instanceLength - [in] Number of extra bytes added to the item/dataset length used when computing the
     *                         padding. This parameter is for instance used to pass the length of the file meta
     *                         header from the DcmFileFormat to the DcmDataset object.
     */
{
    /* if the transfer state of this is not initialized, this is an illegal call */
    if (fTransferState == ERW_notInitialized)
        errorFlag = EC_IllegalCall;
    else
    {
        /* if this is not an illegal call, do something */

        /* assign data set and the meta information header to local variables */
        DcmDataset *dataset = getDataset();
        DcmMetaInfo *metainfo = getMetaInfo();
        /* Determine the transfer syntax which shall be used. Either we use the one which was passed, */
        /* or (if it's an unknown tranfer syntax) we use the data set's original transfer syntax. */
        E_TransferSyntax outxfer = oxfer;
        if (outxfer == EXS_Unknown && dataset)
            outxfer = dataset->getOriginalXfer();
        /* check if the stream reported an error so far */
        errorFlag = outStream.status();
        /* check if we can actually write data to the stream; in certain cases we cannot. */
        if (outxfer == EXS_Unknown || outxfer == EXS_BigEndianImplicit)
            errorFlag = EC_IllegalCall;
        else if (itemList->empty())
            errorFlag = EC_CorruptedData;
        else if (errorFlag.good() && fTransferState != ERW_ready)
        {
            /* in this case we can write data to the stream */

            /* if this function was called for the first time for the dataset object, the transferState is */
            /* still set to ERW_init. In this case, we need to validate the meta header information, set the */
            /* item list pointer to the fist element and we need to set the transfer state to ERW_inWork. */
            if (fTransferState == ERW_init)
            {
                validateMetaInfo(outxfer);
                itemList->seek(ELP_first);
                fTransferState = ERW_inWork;
            }
            /* if the transfer state is set to ERW_inWork, we need to write the */
            /* information which is included in this to the buffer which was passed. */
            if (fTransferState == ERW_inWork)
            {
                /* write meta header information */
                errorFlag = metainfo->write(outStream, outxfer, enctype);
                /* recalculate the instance length */
                instanceLength += metainfo->calcElementLength(outxfer, enctype);
                /* if everything is ok, write the data set */
                if (errorFlag.good())
                    errorFlag = dataset->write(outStream, outxfer, enctype, glenc, padenc, padlen,
                                               subPadlen, instanceLength);
                /* if everything is ok, set the transfer state to ERW_ready */
                if (errorFlag.good())
                    fTransferState = ERW_ready;
            }
        }
        /* in case the transfer syntax which shall be used is indeed the */
        /* BigEndianImplicit transfer syntax dump some error information */
        if (outxfer == EXS_BigEndianImplicit)
        {
            ofConsole.lockCerr() << "Error: DcmFileFormat::write() illegal TransferSyntax(BI) used" << endl;
            ofConsole.unlockCerr();
        }
    }
    /* return result value */
    return errorFlag;
}


// ********************************


OFCondition DcmFileFormat::loadFile(const char *fileName,
                                    const E_TransferSyntax readXfer,
                                    const E_GrpLenEncoding groupLength,
                                    const Uint32 maxReadLength,
                                    const E_FileReadMode readMode)
{
    if (readMode == ERM_dataset)
        return getDataset()->loadFile(fileName, readXfer, groupLength, maxReadLength);

    OFCondition l_error = EC_IllegalParameter;
    /* check parameters first */
    if ((fileName != NULL) && (strlen(fileName) > 0))
    {
        /* open file for input */
        DcmInputFileStream fileStream(fileName);
        /* check stream status */
        l_error = fileStream.status();
        if (l_error.good())
        {
            /* clear this object */
            l_error = clear();
            if (l_error.good())
            {
                /* save old value */
                const E_FileReadMode oldMode = FileReadMode;
                FileReadMode = readMode;
                /* read data from file */
                transferInit();
                l_error = read(fileStream, readXfer, groupLength, maxReadLength);
                transferEnd();
                /* restore old value */
                FileReadMode = oldMode;
            }
        }
    }
    return l_error;
}


OFCondition DcmFileFormat::saveFile(const char *fileName,
                                    const E_TransferSyntax writeXfer,
                                    const E_EncodingType encodingType,
                                    const E_GrpLenEncoding groupLength,
                                    const E_PaddingEncoding padEncoding,
                                    const Uint32 padLength,
                                    const Uint32 subPadLength,
                                    const OFBool isDataset)
{
    if (isDataset)
    {
        return getDataset()->saveFile(fileName, writeXfer, encodingType, groupLength,
            padEncoding, padLength, subPadLength);
    }
    OFCondition l_error = EC_IllegalParameter;
    /* check parameters first */
    if ((fileName != NULL) && (strlen(fileName) > 0))
    {
        /* open file for output */
        DcmOutputFileStream fileStream(fileName);

        /* check stream status */
        l_error = fileStream.status();
        if (l_error.good())
        {
            /* write data to file */
            transferInit();
            l_error = write(fileStream, writeXfer, encodingType, groupLength, padEncoding, padLength, subPadLength);
            transferEnd();
        }
    }
    return l_error;
}


// ********************************


OFCondition DcmFileFormat::insertItem(DcmItem * /*item*/,
                                      const unsigned long /*where*/)
{
    ofConsole.lockCerr() << "Warning: illegal call of DcmFileFormat::insert(DcmItem*,Uin32)" << endl;
    ofConsole.unlockCerr();
    errorFlag = EC_IllegalCall;
    return errorFlag;
}


// ********************************


DcmItem *DcmFileFormat::remove(const unsigned long /*num*/)
{
    ofConsole.lockCerr() << "Warning: illegal call of DcmFileFormat::remove(Uint32)" << endl;
    ofConsole.unlockCerr();
    errorFlag = EC_IllegalCall;
    return NULL;
}


// ********************************


DcmItem *DcmFileFormat::remove(DcmItem* /*item*/)
{
    ofConsole.lockCerr() << "Warning: illegal call of DcmFileFormat::remove(DcmItem*)" << endl;
    ofConsole.unlockCerr();
    errorFlag = EC_IllegalCall;
    return NULL;
}


// ********************************


OFCondition DcmFileFormat::clear()
{
  getMetaInfo()->clear();
  return getDataset()->clear();
}


// ********************************


DcmMetaInfo *DcmFileFormat::getMetaInfo()
{
    errorFlag = EC_Normal;
    DcmMetaInfo *meta = NULL;
    if (itemList->seek_to(0) != NULL && itemList->get()->ident() == EVR_metainfo)
        meta = OFstatic_cast(DcmMetaInfo *, itemList->get());
    else
        errorFlag = EC_IllegalCall;
    return meta;
}


// ********************************


DcmDataset *DcmFileFormat::getDataset()
{
    errorFlag = EC_Normal;
    DcmDataset *data = NULL;
    if (itemList->seek_to(1) != NULL && itemList->get()->ident() == EVR_dataset)
        data = OFstatic_cast(DcmDataset *, itemList->get());
    else
        errorFlag = EC_IllegalCall;
    return data;
}


// ********************************


DcmDataset *DcmFileFormat::getAndRemoveDataset()
{
    errorFlag = EC_Normal;
    DcmDataset *data = NULL;
    if (itemList->seek_to(1) != NULL && itemList->get()->ident() == EVR_dataset)
    {
        data = OFstatic_cast(DcmDataset *, itemList->remove());
        DcmDataset *Dataset = new DcmDataset();
        DcmSequenceOfItems::itemList->insert(Dataset, ELP_last);
    }
    else
        errorFlag = EC_IllegalCall;
    return data;
}


/*
** CVS/RCS Log:
** $Log: dcfilefo.cc,v $
** Revision 1.42  2005/12/08 15:41:10  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.41  2005/12/02 11:57:44  joergr
** Fixed minor bug.
**
** Revision 1.40  2005/12/02 08:53:36  joergr
** Added new file read mode that makes it possible to distinguish between DICOM
** files, datasets and other non-DICOM files.  For this reason, the last
** parameter of method loadFile() changed from OFBool to E_FileReadMode.
**
** Revision 1.39  2005/11/28 15:53:13  meichel
** Renamed macros in dcdebug.h
**
** Revision 1.38  2005/11/07 17:22:33  meichel
** Implemented DcmFileFormat::clear(). Now the loadFile methods in class
**   DcmFileFormat and class DcmDataset both make sure that clear() is called
**   before new data is read, allowing for a re-use of the object.
**
** Revision 1.37  2004/08/03 11:41:09  meichel
** Headers libc.h and unistd.h are now included via ofstdinc.h
**
** Revision 1.36  2004/02/04 16:32:01  joergr
** Adapted type casts to new-style typecast operators defined in ofcast.h.
** Removed acknowledgements with e-mail addresses from CVS log.
**
** Revision 1.35  2003/04/01 14:57:20  joergr
** Added support for XML namespaces.
**
** Revision 1.34  2003/03/21 13:08:04  meichel
** Minor code purifications for warnings reported by MSVC in Level 4
**
** Revision 1.33  2002/12/06 13:09:26  joergr
** Enhanced "print()" function by re-working the implementation and replacing
** the boolean "showFullData" parameter by a more general integer flag.
** Made source code formatting more consistent with other modules/files.
**
** Revision 1.32  2002/11/27 12:06:47  meichel
** Adapted module dcmdata to use of new header file ofstdinc.h
**
** Revision 1.31  2002/08/27 16:55:47  meichel
** Initial release of new DICOM I/O stream classes that add support for stream
**   compression (deflated little endian explicit VR transfer syntax)
**
** Revision 1.30  2002/08/20 12:18:48  meichel
** Changed parameter list of loadFile and saveFile methods in class
**   DcmFileFormat. Removed loadFile and saveFile from class DcmObject.
**
** Revision 1.29  2002/05/14 08:21:14  joergr
** Renamed some element names.
**
** Revision 1.28  2002/04/25 10:15:35  joergr
** Added support for XML output of DICOM objects.
**
** Revision 1.27  2002/04/16 13:43:17  joergr
** Added configurable support for C++ ANSI standard includes (e.g. streams).
**
** Revision 1.26  2002/04/11 12:27:10  joergr
** Added new methods for loading and saving DICOM files.
**
** Revision 1.25  2002/01/08 10:43:45  joergr
** Corrected spelling of function dcmGenerateUniqueIdentifier().
**
** Revision 1.24  2001/11/01 14:55:37  wilkens
** Added lots of comments.
**
** Revision 1.23  2001/09/28 14:20:26  joergr
** Added "#include <iomanip.h>" to keep gcc 3.0 quiet.
**
** Revision 1.22  2001/09/25 17:19:50  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.21  2001/06/01 15:49:04  meichel
** Updated copyright header
**
** Revision 1.20  2000/04/14 15:55:05  meichel
** Dcmdata library code now consistently uses ofConsole for error output.
**
** Revision 1.19  2000/03/08 16:26:35  meichel
** Updated copyright header.
**
** Revision 1.18  2000/03/03 14:05:33  meichel
** Implemented library support for redirecting error messages into memory
**   instead of printing them to stdout/stderr for GUI applications.
**
** Revision 1.17  2000/02/23 15:11:53  meichel
** Corrected macro for Borland C++ Builder 4 workaround.
**
** Revision 1.16  2000/02/10 10:52:18  joergr
** Added new feature to dcmdump (enhanced print method of dcmdata): write
** pixel data/item value fields to raw files.
**
** Revision 1.15  2000/02/01 10:12:07  meichel
** Avoiding to include <stdlib.h> as extern "C" on Borland C++ Builder 4,
**   workaround for bug in compiler header files.
**
** Revision 1.14  1999/03/31 09:25:28  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.13  1998/07/15 15:51:56  joergr
** Removed several compiler warnings reported by gcc 2.8.1 with
** additional options, e.g. missing copy constructors and assignment
** operators, initialization of member variables in the body of a
** constructor instead of the member initialization list, hiding of
** methods by use of identical names, uninitialized member variables,
** missing const declaration of char pointers. Replaced tabs by spaces.
**
** Revision 1.12  1997/07/21 08:25:28  andreas
** - Replace all boolean types (BOOLEAN, CTNBOOLEAN, DICOM_BOOL, BOOL)
**   with one unique boolean type OFBool.
**
** Revision 1.11  1997/07/03 15:09:57  andreas
** - removed debugging functions Bdebug() and Edebug() since
**   they write a static array and are not very useful at all.
**   Cdebug and Vdebug are merged since they have the same semantics.
**   The debugging functions in dcmdata changed their interfaces
**   (see dcmdata/include/dcdebug.h)
**
** Revision 1.10  1997/05/27 13:48:59  andreas
** - Add method canWriteXfer to class DcmObject and all derived classes.
**   This method checks whether it is possible to convert the original
**   transfer syntax to an new transfer syntax. The check is used in the
**   dcmconv utility to prohibit the change of a compressed transfer
**   syntax to a uncompressed.
**
** Revision 1.9  1997/05/16 08:23:54  andreas
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
**
** Revision 1.8  1997/04/18 08:17:17  andreas
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
** Revision 1.7  1996/08/05 08:46:11  andreas
** new print routine with additional parameters:
**         - print into files
**         - fix output length for elements
** corrected error in search routine with parameter ESM_fromStackTop
**
** Revision 1.6  1996/07/29 15:47:00  andreas
** add method getAndRemoveDataset to remove connection of dataset from fileformat
**
** Revision 1.5  1996/04/12 13:17:22  andreas
** Minor changes to support DEC ALPHA and DEC MIPS
**
** Revision 1.4  1996/03/11 14:48:04  hewett
** When creating a DICOM file, if a SOP Class UID is present in the
** data set then it will be used.  Otherwise, a private UID will now be used
** in the meta header attribute MediaStorageSOPClass UID.  Previously this
** attribute was left empty (an error for type 1 attributes).
**
** Revision 1.3  1996/01/05 13:27:37  andreas
** - changed to support new streaming facilities
** - unique read/write methods for file and block transfer
** - more cleanups
**
**
*/
