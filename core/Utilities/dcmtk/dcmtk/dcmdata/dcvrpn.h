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
 *  Author:  Gerd Ehlers
 *
 *  Purpose: Interface of class DcmPersonName
 *
 *  Last Update:      $Author: joergr $
 *  Update Date:      $Date: 2010-11-05 09:34:11 $
 *  CVS/RCS Revision: $Revision: 1.26 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */


#ifndef DCVRPN_H
#define DCVRPN_H

#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */

#include "dcmtk/dcmdata/dcchrstr.h"


/** a class representing the DICOM value representation 'Person Name' (PN)
 */
class DcmPersonName
  : public DcmCharString
{

  public:

    /** constructor.
     *  Create new element from given tag and length.
     *  @param tag DICOM tag for the new element
     *  @param len value length for the new element
     */
    DcmPersonName(const DcmTag &tag,
                  const Uint32 len = 0);

    /** copy constructor
     *  @param old element to be copied
     */
    DcmPersonName(const DcmPersonName &old);

    /** destructor
     */
    virtual ~DcmPersonName();

    /** assignment operator
     *  @param obj element to be assigned/copied
     *  @return reference to this object
     */
    DcmPersonName &operator=(const DcmPersonName &obj);

    /** clone method
     *  @return deep copy of this object
     */
    virtual DcmObject *clone() const
    {
      return new DcmPersonName(*this);
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
     *  @return type identifier of this class (EVR_PN)
     */
    virtual DcmEVR ident() const;

    /** check whether stored value conforms to the VR and to the specified VM
     *  @param vm value multiplicity (according to the data dictionary) to be checked for.
     *    (valid values: "1", "1-2", "1-3", "1-8", "1-99", "1-n", "2", "2-n", "2-2n",
     *                   "3", "3-n", "3-3n", "4", "6", "9", "16", "32")
     *  @param oldFormat support old ACR/NEMA format if OFTrue (no '^' separator)
     *  @return status of the check, EC_Normal if value is correct, an error code otherwise
     */
    virtual OFCondition checkValue(const OFString &vm = "1-n",
                                   const OFBool oldFormat = OFFalse);

    /** get a copy of a particular string component
     *  @param stringVal variable in which the result value is stored
     *  @param pos index of the value in case of multi-valued elements (0..vm-1)
     *  @param normalize delete leading and trailing spaces if OFTrue
     *  @return status, EC_Normal if successful, an error code otherwise
     */
    virtual OFCondition getOFString(OFString &stringVal,
                                    const unsigned long pos,
                                    OFBool normalize = OFTrue);

    /** get name components from the element value.
     *  The DICOM PN consists of up to three component groups separated by a "=". The
     *  supported format is "[CG0][=CG1][=CG2]" where the brackets enclose optional
     *  parts and CG0 is a single-byte character representation, CG1 an ideographic
     *  representation, and CG2 a phonetic representation of the name.
     *  Each component group may consist of up to five components separated by a "^".
     *  The format is "[lastName[^firstName[^middleName[^namePrefix[^nameSuffix]]]]";
     *  each component might be empty.
     *  If this function fails the result variables are cleared automatically. If the
     *  format is valid but does not comply with the above described scheme ("=" and "^")
     *  the full person name is returned in the 'lastName' variable.
     *  @param lastName reference to string variable where the "last name" is stored
     *  @param firstName reference to string variable where the "first name" is stored
     *  @param middleName reference to string variable where the "middle name" is stored
     *  @param namePrefix reference to string variable where the "name prefix" is stored
     *  @param nameSuffix reference to string variable where the "name suffix" is stored
     *  @param pos index of the element component in case of value multiplicity (0..vm-1)
     *  @param componentGroup index of the component group (0..2) to be used, see above
     *  @return EC_Normal upon success, an error code otherwise
     */
    OFCondition getNameComponents(OFString &lastName,
                                  OFString &firstName,
                                  OFString &middleName,
                                  OFString &namePrefix,
                                  OFString &nameSuffix,
                                  const unsigned long pos = 0,
                                  const unsigned int componentGroup = 0);

    /** get current element value as a formatted/readable name.
     *  The current element value is expected to be in DICOM PN format as described above.
     *  The output format is "[namePrefix][ firstName][ middleName][ lastName][, nameSuffix]";
     *  the delimiters (" " and ", ") are only inserted if required.
     *  If this function fails the result variable 'formattedName' is cleared automatically.
     *  @param formattedName reference to string variable where the result is stored
     *  @param pos index of the element component in case of value multiplicity (0..vm-1)
     *  @param componentGroup index of the component group (0..2) to be used, see above
     *  @return EC_Normal upon success, an error code otherwise
     */
    OFCondition getFormattedName(OFString &formattedName,
                                 const unsigned long pos = 0,
                                 const unsigned int componentGroup = 0);


    /** put element value from specified name components.
     *  The stored format is "[lastName[^firstName[^middleName[^namePrefix[^nameSuffix]]]]]",
     *  i.e. a DICOM Person Name (PN). Component groups are not (yet) supported.
     *  If this function fails the currently stored value is not modified.
     *  @param lastName reference to string variable where the "last name" is stored
     *  @param firstName reference to string variable where the "first name" is stored
     *  @param middleName reference to string variable where the "middle name" is stored
     *  @param namePrefix reference to string variable where the "name prefix" is stored
     *  @param nameSuffix reference to string variable where the "name suffix" is stored
     *  @return EC_Normal upon success, an error code otherwise
     */
    OFCondition putNameComponents(const OFString &lastName,
                                  const OFString &firstName,
                                  const OFString &middleName,
                                  const OFString &namePrefix,
                                  const OFString &nameSuffix);

    /* --- static helper functions --- */

    /** get name components from specified DICOM person name.
     *  The DICOM PN consists of up to three component groups separated by a "=". The
     *  supported format is "[CG0][=CG1][=CG2]" where the brackets enclose optional
     *  parts and CG0 is a single-byte character representation, CG1 an ideographic
     *  representation, and CG2 a phonetic representation of the name.
     *  Each component group may consist of up to five components separated by a "^".
     *  The format is "[lastName[^firstName[^middleName[^namePrefix[^nameSuffix]]]]";
     *  each component might be empty.
     *  If this function fails the result variables are cleared automatically. If the
     *  format is valid but does not comply with the above described scheme ("=" and "^")
     *  the full person name is returned in the 'lastName' variable.
     *  @param dicomName string value in DICOM PN format to be split into components
     *  @param lastName reference to string variable where the "last name" is stored
     *  @param firstName reference to string variable where the "first name" is stored
     *  @param middleName reference to string variable where the "middle name" is stored
     *  @param namePrefix reference to string variable where the "name prefix" is stored
     *  @param nameSuffix reference to string variable where the "name suffix" is stored
     *  @param componentGroup index of the component group (0..2) to be used, see above
     *  @return EC_Normal upon success, an error code otherwise
     */
    static OFCondition getNameComponentsFromString(const OFString &dicomName,
                                                   OFString &lastName,
                                                   OFString &firstName,
                                                   OFString &middleName,
                                                   OFString &namePrefix,
                                                   OFString &nameSuffix,
                                                   const unsigned int componentGroup = 0);

    /** get specified DICOM person name as a formatted/readable name.
     *  The specified 'dicomName' is expected to be in DICOM PN format as described above.
     *  The output format is "[namePrefix][ firstName][ middleName][ lastName][, nameSuffix]";
     *  the delimiters (" " and ", ") are only inserted if required.
     *  If this function fails the result variable 'formattedName' is cleared automatically.
     *  @param dicomName string value in DICOM PN format to be converted to readable format
     *  @param formattedName reference to string variable where the result is stored
     *  @param componentGroup index of the component group (0..2) to be used, see above
     *  @return EC_Normal upon success, an error code otherwise
     */
    static OFCondition getFormattedNameFromString(const OFString &dicomName,
                                                  OFString &formattedName,
                                                  const unsigned int componentGroup = 0);

    /** get formatted/readable name from specified name components.
     *  The output format is "[namePrefix][ firstName][ middleName][ lastName][, nameSuffix]";
     *  the delimiters (" " and ", ") are only inserted if required.
     *  If this function fails the result variable 'formattedName' is cleared automatically.
     *  @param lastName reference to string variable where the "last name" is stored
     *  @param firstName reference to string variable where the "first name" is stored
     *  @param middleName reference to string variable where the "middle name" is stored
     *  @param namePrefix reference to string variable where the "name prefix" is stored
     *  @param nameSuffix reference to string variable where the "name suffix" is stored
     *  @param formattedName reference to string variable where the result is stored
     *  @return always returns EC_Normal
     */
    static OFCondition getFormattedNameFromComponents(const OFString &lastName,
                                                      const OFString &firstName,
                                                      const OFString &middleName,
                                                      const OFString &namePrefix,
                                                      const OFString &nameSuffix,
                                                      OFString &formattedName);

    /** get DICOM Person Name (PN) from specified name components.
     *  The output format is "[lastName[^firstName[^middleName[^namePrefix[^nameSuffix]]]]]".
     *  Component groups are not (yet) supported.
     *  If this function fails the result variable 'dicomName' is cleared automatically.
     *  @param lastName reference to string variable where the "last name" is stored
     *  @param firstName reference to string variable where the "first name" is stored
     *  @param middleName reference to string variable where the "middle name" is stored
     *  @param namePrefix reference to string variable where the "name prefix" is stored
     *  @param nameSuffix reference to string variable where the "name suffix" is stored
     *  @param dicomName reference to string variable where the result is stored
     *  @return always returns EC_Normal
     */
    static OFCondition getStringFromNameComponents(const OFString &lastName,
                                                   const OFString &firstName,
                                                   const OFString &middleName,
                                                   const OFString &namePrefix,
                                                   const OFString &nameSuffix,
                                                   OFString &dicomName);

    /** check whether given string value conforms to the VR "PN" (Person Name)
     *  and to the specified VM.
     *  @param value string value to be checked (possibly multi-valued)
     *  @param vm value multiplicity (according to the data dictionary) to be checked for.
     *    (valid values: "1", "1-2", "1-3", "1-8", "1-99", "1-n", "2", "2-n", "2-2n",
     *                   "3", "3-n", "3-3n", "4", "6", "9", "16", "32")
     *  @param oldFormat support old ACR/NEMA name format if OFTrue (i.e. without "^" delimiters)
     *  @return status of the check, EC_Normal if value is correct, an error code otherwise
     */
    static OFCondition checkStringValue(const OFString &value,
                                        const OFString &vm = "1-n",
                                        const OFBool oldFormat = OFTrue);
};


#endif // DCVRPN_H


/*
** CVS/RCS Log:
** $Log: dcvrpn.h,v $
** Revision 1.26  2010-11-05 09:34:11  joergr
** Added support for checking the value multiplicity "9" (see Supplement 131).
**
** Revision 1.25  2010-10-14 13:15:43  joergr
** Updated copyright header. Added reference to COPYRIGHT file.
**
** Revision 1.24  2010-04-23 15:26:13  joergr
** Specify an appropriate default value for the "vm" parameter of checkValue().
**
** Revision 1.23  2010-04-23 14:25:27  joergr
** Added new method to all VR classes which checks whether the stored value
** conforms to the VR definition and to the specified VM.
**
** Revision 1.22  2010-04-22 09:31:30  joergr
** Revised misleading parameter documentation for the checkValue() method.
**
** Revision 1.21  2010-04-22 08:59:10  joergr
** Added support for further VM values ("1-8", "1-99", "16", "32") to be checked.
**
** Revision 1.20  2009-08-03 09:05:30  joergr
** Added methods that check whether a given string value conforms to the VR and
** VM definitions of the DICOM standards.
**
** Revision 1.19  2008-07-17 11:19:49  onken
** Updated copyFrom() documentation.
**
** Revision 1.18  2008-07-17 10:30:23  onken
** Implemented copyFrom() method for complete DcmObject class hierarchy, which
** permits setting an instance's value from an existing object. Implemented
** assignment operator where necessary.
**
** Revision 1.17  2005-12-08 16:29:05  meichel
** Changed include path schema for all DCMTK header files
**
** Revision 1.16  2004/07/01 12:28:25  meichel
** Introduced virtual clone method for DcmObject and derived classes.
**
** Revision 1.15  2003/05/20 08:56:20  joergr
** Added methods and static functions to compose a DICOM Person Name from five
** name components.
**
** Revision 1.14  2002/12/06 12:49:17  joergr
** Enhanced "print()" function by re-working the implementation and replacing
** the boolean "showFullData" parameter by a more general integer flag.
** Added doc++ documentation.
** Made source code formatting more consistent with other modules/files.
**
** Revision 1.13  2002/04/25 09:56:19  joergr
** Removed getOFStringArray() implementation.
**
** Revision 1.12  2001/10/10 15:17:38  joergr
** Updated comments.
**
** Revision 1.11  2001/10/01 15:01:39  joergr
** Introduced new general purpose functions to get/set person names, date, time
** and date/time.
**
** Revision 1.10  2001/09/25 17:19:33  meichel
** Adapted dcmdata to class OFCondition
**
** Revision 1.9  2001/06/01 15:48:51  meichel
** Updated copyright header
**
** Revision 1.8  2000/03/08 16:26:25  meichel
** Updated copyright header.
**
** Revision 1.7  1999/03/31 09:25:04  meichel
** Updated copyright header in module dcmdata
**
** Revision 1.6  1998/11/12 16:47:52  meichel
** Implemented operator= for all classes derived from DcmObject.
**
** Revision 1.5  1997/09/11 15:13:16  hewett
** Modified getOFString method arguments by removing a default value
** for the pos argument.  By requiring the pos argument to be provided
** ensures that callers realise getOFString only gets one component of
** a multi-valued string.
**
** Revision 1.4  1997/08/29 08:32:43  andreas
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
** Revision 1.3  1996/01/05 13:23:08  andreas
** - changed to support new streaming facilities
** - more cleanups
** - merged read / write methods for block and file transfer
**
*/
