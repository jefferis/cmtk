/*
 *
 *  Copyright (C) 1994-2010, OFFIS e.V.
 *  All rights reserved.  See COPYRIGHT file for details.
 *
 *  This software and supporting documentation were developed by
 *
 *    OFFIS e.V.
 *    R&D Division Health
 *    Escherweg 2
 *    D-26121 Oldenburg, Germany
 *
 *
 *  Module:  dcmdata
 *
 *  Author:  Gerd Ehlers, Andreas Barth
 *
 *  Purpose: Interface of class DcmAttributeTag
 *
 *  Last Update:      $Author: joergr $
 *  Update Date:      $Date: 2010-11-05 09:34:11 $
 *  CVS/RCS Revision: $Revision: 1.30 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#ifndef DCVRAT_H
#define DCVRAT_H

#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */

#include "dcmtk/dcmdata/dcelem.h"


/** a class representing the DICOM value representation 'Attribute Tag' (AT)
 */
class DcmAttributeTag
  : public DcmElement
{

  public:

    /** constructor.
     *  Create new element from given tag and length.
     *  @param tag DICOM tag for the new element
     *  @param len value length for the new element
     */
    DcmAttributeTag(const DcmTag &tag,
                    const Uint32 len = 0);

    /** copy constructor
     *  @param old element to be copied
     */
    DcmAttributeTag(const DcmAttributeTag &old);

    /** destructor
     */
    virtual ~DcmAttributeTag();

    /** assignment operator
     *  @param obj element to be assigned/copied
     *  @return reference to this object
     */
    DcmAttributeTag &operator=(const DcmAttributeTag &obj);

    /** clone method
     *  @return deep copy of this object
     */
    virtual DcmObject *clone() const
    {
      return new DcmAttributeTag(*this);
    }

    /** Virtual object copying. This method can be used for DcmObject
     *  and derived classes to get a deep copy of an object. Internally
     *  the assignment operator is called if the given DcmObject parameter
     *  is of the same type as "this" object instance. If not, an error
     *  is returned. This function permits copying an object by value
     *  in a virtual way which therefore is different to just calling the
     *  assignment operator of DcmElement which could result in slicing
     *  the object.
     *  @param rhs - [in] The instance to copy from. Has to be of the same
     *                class type as "this" object
     *  @return EC_Normal if copying was successful, error otherwise
     */
    virtual OFCondition copyFrom(const DcmObject& rhs);

    /** get element type identifier
     *  @return type identifier of this class (EVR_AT)
     */
    virtual DcmEVR ident() const;

    /** check whether stored value conforms to the VR and to the specified VM
     *  @param vm value multiplicity (according to the data dictionary) to be checked for.
     *    (valid values: "1", "1-2", "1-3", "1-8", "1-99", "1-n", "2", "2-n", "2-2n",
     *                   "3", "3-n", "3-3n", "4", "9", "6", "16", "32")
     *  @param oldFormat parameter not used for this VR (only for DA, TM, PN)
     *  @return status of the check, EC_Normal if value is correct, an error code otherwise
     */
    virtual OFCondition checkValue(const OFString &vm = "1-n",
                                   const OFBool oldFormat = OFFalse);

    /** get value multiplicity
     *  @return number of tag value pairs (group,element)
     */
    virtual unsigned long getVM();

    /** print element to a stream.
     *  The output format of the value is a backslash separated sequence of group and
     *  element value pairs, e.g. "(0008,0020)\(0008,0030)"
     *  @param out output stream
     *  @param flags optional flag used to customize the output (see DCMTypes::PF_xxx)
     *  @param level current level of nested items. Used for indentation.
     *  @param pixelFileName not used
     *  @param pixelCounter not used
     */
    virtual void print(STD_NAMESPACE ostream&out,
                       const size_t flags = 0,
                       const int level = 0,
                       const char *pixelFileName = NULL,
                       size_t *pixelCounter = NULL);

    /** get particular tag value
     *  @param tagVal reference to result variable (cleared in case of error)
     *  @param pos index of the value to be retrieved (0..vm-1)
     *  @return status status, EC_Normal if successful, an error code otherwise
     */
    virtual OFCondition getTagVal(DcmTagKey &tagVal,
                                  const unsigned long pos = 0);

    /** get reference to stored integer data.
     *  The array entries with an even-numbered index contain the group numbers
     *  and the odd entries contain the element numbers (see "putUint16Array()").
     *  The number of entries is twice as large as the return value of "getVM()".
     *  @param uintVals reference to result variable
     *  @return status, EC_Normal if successful, an error code otherwise
     */
    virtual OFCondition getUint16Array(Uint16 *&uintVals);

    /** get specified value as a character string.
     *  The output format is "(gggg,eeee)" where "gggg" is the hexa-decimal group
     *  number and "eeee" the hexa-decimal element number of the attribute tag.
     *  @param stringVal variable in which the result value is stored
     *  @param pos index of the value in case of multi-valued elements (0..vm-1)
     *  @param normalize not used
     *  @return status, EC_Normal if successful, an error code otherwise
     */
    virtual OFCondition getOFString(OFString &stringVal,
                                    const unsigned long pos,
                                    OFBool normalize = OFTrue);

    /** set particular tag value
     *  @param tagVal tag value to be set
     *  @param pos index of the value to be set (0 = first position)
     *  @return status, EC_Normal if successful, an error code otherwise
     */
    virtual OFCondition putTagVal(const DcmTagKey &tagVal,
                                  const unsigned long pos = 0);

    /** set element value to given integer array data.
     *  The array entries with an even-numbered index are expected to contain the
     *  group numbers and the odd entries to contain the element numbers, e.g.
     *  {0x0008, 0x0020, 0x0008, 0x0030}. This function uses the same format as
     *  "getUint16Array()".
     *  @param uintVals unsigned integer data to be set
     *  @param numUints number of integer values to be set (should be even)
     *  @return status, EC_Normal if successful, an error code otherwise
     */
    virtual OFCondition putUint16Array(const Uint16 *uintVals,
                                       const unsigned long numUints);

    /** set element value from the given character string.
     *  The input string is expected to be a backslash separated sequence of
     *  attribute tags, e.g. "(0008,0020)\(0008,0030)". This is the same format
     *  as used by "print()".
     *  @param stringVal input character string
     *  @return status, EC_Normal if successful, an error code otherwise
     */
    virtual OFCondition putString(const char *stringVal);

    /** check the currently stored element value
     *  @param autocorrect correct value length if OFTrue
     *  @return status, EC_Normal if value length is correct, an error code otherwise
     */
    virtual OFCondition verify(const OFBool autocorrect = OFFalse);

    /* --- static helper functions --- */

    /** check whether given string value conforms to the VR "AT" (Attribute Tag)
     *  and to the specified VM.
     *  @param value string value to be checked (possibly multi-valued)
     *  @param vm value multiplicity (according to the data dictionary) to be checked for.
     *    (valid values: "1", "1-2", "1-3", "1-8", "1-99", "1-n", "2", "2-n", "2-2n",
     *                   "3", "3-n", "3-3n", "4", "6", "9", "16", "32")
     *  @return status of the check, EC_Normal if value is correct, an error code otherwise
     */
    static OFCondition checkStringValue(const OFString &value,
                                        const OFString &vm = "1-n");
};


#endif // DCVRAT_H


/*
** CVS/RCS Log:
** $Log: dcvrat.h,v $
** Revision 1.30  2010-11-05 09:34:11  joergr
** Added support for checking the value multiplicity "9" (see Supplement 131).
**
** Revision 1.29  2010-10-14 13:15:42  joergr
** Updated copyright header. Added reference to COPYRIGHT file.
**
** Revision 1.28  2010-04-23 15:26:12  joergr
** Specify an appropriate default value for the "vm" parameter of checkValue().
**
** Revision 1.27  2010-04-23 14:25:27  joergr
** Added new method to all VR classes which checks whether the stored value
** conforms to the VR definition and to the specified VM.
**
** Revision 1.26  2009-11-04 09:58:07  uli
** Switched to logging mechanism provided by the "new" oflog module
**
** Revision 1.25  2008-07-17 11:19:49  onken
** Updated copyFrom() documentation.
**
** Revision 1.24  2008-07-17 10:30:23  onken
** Implemented copyFrom() method for complete DcmObject class hierarchy, which
** permits setting an instance's value from an existing object. Implemented
** assignment operator where necessary.
**
** Revision 1.23  2006-08-15 15:49:56  meichel
** Updated all code in module dcmdata to correctly compile when
**   all standard C++ classes remain in namespace std.
**
** Revision 1.22  2005/12/08 16:28:53  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.21  2004/07/01 12:28:25  meichel
** Introduced virtual clone method for DcmObject and derived classes.
**
** Revision 1.20  2003/06/12 13:31:46  joergr
** Fixed inconsistent API documentation reported by Doxygen.
**
** Revision 1.19  2002/12/06 12:49:14  joergr
** Enhanced "print()" function by re-working the implementation and replacing
** the boolean "showFullData" parameter by a more general integer flag.
** Added doc++ documentation.
** Made source code formatting more consistent with other modules/files.
**
** Revision 1.18  2002/04/25 09:50:38  joergr
** Added getOFString() implementation.
**
** Revision 1.17  2001/09/25 17:19:30  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.16  2001/06/01 15:48:48  meichel
** Updated copyright header
**
** Revision 1.15  2000/04/14 15:31:34  meichel
** Removed default value from output stream passed to print() method.
**   Required for use in multi-thread environments.
**
** Revision 1.14  2000/03/08 16:26:22  meichel
** Updated copyright header.
**
** Revision 1.13  2000/03/03 14:05:26  meichel
** Implemented library support for redirecting error messages into memory
**   instead of printing them to stdout/stderr for GUI applications.
**
** Revision 1.12  2000/02/10 10:50:54  joergr
** Added new feature to dcmdump (enhanced print method of dcmdata): write
** pixel data/item value fields to raw files.
**
** Revision 1.11  1999/03/31 09:24:56  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.10  1998/11/12 16:47:45  meichel
** Implemented operator= for all classes derived from DcmObject.
**
** Revision 1.9  1997/07/21 08:25:13  andreas
** - Replace all boolean types (BOOLEAN, CTNBOOLEAN, DICOM_BOOL, BOOL)
**   with one unique boolean type OFBool.
**
** Revision 1.8  1997/04/18 08:13:29  andreas
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
** Revision 1.7  1996/08/05 08:45:29  andreas
** new print routine with additional parameters:
**         - print into files
**         - fix output length for elements
** corrected error in search routine with parameter ESM_fromStackTop
**
** Revision 1.6  1996/04/16 16:01:37  andreas
** - put methods for AttributeTag with DcmTagKey Parameter
** - better support for NULL values
**
** Revision 1.5  1996/01/29 13:38:15  andreas
** - new put method for every VR to put value as a string
** - better and unique print methods
**
** Revision 1.4  1996/01/09 11:06:17  andreas
** New Support for Visual C++
** Correct problems with inconsistent const declarations
**
** Revision 1.3  1996/01/05 13:23:03  andreas
** - changed to support new streaming facilities
** - more cleanups
** - merged read / write methods for block and file transfer
**
*/
