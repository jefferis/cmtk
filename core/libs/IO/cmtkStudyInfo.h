/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkStudyInfo_h_included_
#define __cmtkStudyInfo_h_included_

#include <stdlib.h>
#include <string.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

enum
{
  INFO_PATNAME = 0,
  INFO_PATID = 1,
  INFO_PATBIRTHD = 2,
  INFO_PATSEX = 3,
  INFO_PATMAIDNAME = 4,
  INFO_PATAGE = 5,
  INFO_PATSIZE = 6,
  INFO_PATWEIGHT = 7,
  INFO_PATADDRESS = 8,
  INFO_PATPLANID = 9,
  INFO_PATMOTHMAIDNAME = 10,

  INFO_ACQSTUDYDATE = 11,
  INFO_ACQSERIESDATE = 12,
  INFO_ACQSTUDYTIME = 13,
  INFO_ACQMODALITY = 14,
  INFO_ACQMANUFACT = 15,
  INFO_ACQINSTID = 16,
  INFO_ACQREFPHYS = 17,
  INFO_ACQSTATID = 18,
  INFO_ACQINSTDEP = 19,
  INFO_ACQPROCDESCR = 20,
  INFO_ACQATTPHYS = 21,
  INFO_ACQRADIOL = 22,
  INFO_ACQOPERID = 23,
  INFO_ACQMANUMODEL = 24,
  INFO_ACQSCANSEQ = 25,
  INFO_ACQSLICETHICK = 26,
  INFO_ACQSLICESPACING = 27,
  INFO_ACQGANTILT = 28,
  INFO_ACQPATPOSIT = 29,

  INFO_ACQNUMBER = 30,
  INFO_RELPATORIENT = 31,
  INFO_RELIMGPOSITION = 32,
  INFO_RELIMGPOSITPAT = 33,
  INFO_RELIMGORIENT = 34,
  INFO_RELIMGORIENTPAT = 35,
  INFO_RELIMGLOCATION = 36,
  INFO_RELIMGLATERAL = 37,

  INFO_NETUID = 38
};

#define INFO_NUMANONS 11
#define INFO_NUMFIELDS 38

#define StudyInfo_NoAnon 0
#define StudyInfo_Anon 1

#define StudyInfo_NoDupString 0
#define StudyInfo_DupString 1

/** DICOM image series information.
 * This class handles the various information elements extracted from DICOM
 * files when they are read.
 */
class StudyInfo 
{
protected:
  /// Array of tag values.
  char* Elements[INFO_NUMFIELDS];

  /** Table of tag assignments.
   * This table defines the standard partitioning of the tags stored into
   * geometry- and content-dependent data.
   */
  static const char DefaultTagAssignmentsTable[];

public:
  /// Offset of the original slice images in mm.
  double SliceOffset;

  /** Slice location direction flag.
   * This flag is -1 if the slice locations decrease for increasing slice
   * indices, and +1 if slice locations increase.
   */
  int SliceDirection;

  /** Clear object fields.
   * All tag value pointers and internal variables are reset to NULL (or 0).
   * No memory is freed, so this function should not be called after the
   * object has been used to store string tags.
   */
  void Clear () {
    memset( Elements, 0, sizeof(Elements) );
    SliceOffset = 0;
    SliceDirection = 0;
  };

  /// Default Constructor.
  StudyInfo () {
    Clear();
  };

  /** Merge Constructor.
   * This constructor merges geometry and content dependent information from
   * to other objects. An optional assignment table defines which tags are to
   * be taken from which source object.
   *@param geometry The object containing the tags of the geometry defining
   * study (ie. reference study).
   *@param content The content dependent tags are taken from this object
   * (ie. reslice study).
   *@param assignments This parameter is a pointer to a table containing one
   * char entry per information tag. Its entries define which tag is to be
   * taken from which study. If this parameter is NULL or omitted, a standard
   * assignment table (StudyInfo::DefaultTagAssignmentTable) is used.
   */
  StudyInfo ( const StudyInfo&, const StudyInfo&, 
		 const char* = NULL );

  /** Destructor.
   * All memory allocated is freed.
   */
  virtual ~StudyInfo();

  /** Get tag value.
   *@param idx The index of the desired tag. Use INFO_xxx constants to specify
   * this value.
   *@param anon If this flag is non-zero, returned data is anonymized. This
   * means that for all data allowing identification of the patient, the value
   * "anonymized" is returned.
   */
  const char* GetField ( const int, const int = 0 ) const;

  /** Return the pointer to an information element.
   * This function does not duplicate the requested information. The returned
   * pointer can therefore only be used for reading.
   */
  const char* ReturnFieldPtr ( const int ) const;

  /** Set tag value.
   *@param idx The index of the tag the value of which to set. Use INFO_xxx
   * constants to specify this value.
   *@param str The string to store as the tag value.
   *@param dup If this flag is non-zero, the tag value string is duplicated
   * before it is stored. If this flag is given as zero, the string pointer
   * must not be used by the calling function anymore but will be freed by a
   * call to free() on destruction of this object.
   */
  void SetField ( const int, char *const, const bool = true );

  /** Merge Operator.
   * This function merges geometry and content dependent information from
   * to other objects. An optional assignment table defines which tags are to
   * be taken from which source object.
   *@param geometry The object containing the tags of the geometry defining
   * study (ie. reference study).
   *@param content The content dependent tags are taken from this object
   * (ie. reslice study).
   *@param assignments This parameter is a pointer to a table containing one
   * char entry per information tag. Its entries define which tag is to be
   * taken from which study. If this parameter is NULL or omitted, a standard
   * assignment table (StudyInfo::DefaultTagAssignmentTable) is used.
   *@return A reference to this object.
   */
  StudyInfo& Merge ( const StudyInfo&, const StudyInfo&, const char* = NULL );

  /// Create a duplicate of this object.
  StudyInfo* Clone() const;

  /// Return default text for tags with unknown values.
  static char* StrNotAvailable () 
  {
    return NULL;
  }
  
  /// Return text for anonymized tags.
  static const char* StrAnonymized () 
  {
    return "anonymized";
  }
  
  /// Return modification device ID.
  static const char* ModDevID () 
  {
    return "IGSL";
  }
  
  /** Assignment operator.
   * This operator copies all information in the given object. All strings are
   * duplicated, so both objects are completely independent afterwards.
   */
  StudyInfo& operator= ( const StudyInfo& other );
};

/** Temporary study info storage.
 * This variable is used for example when a StudyInfo object is needed for a
 * function call.
 */
extern StudyInfo LastImageStudyInfo;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStudyInfo_h_included_
