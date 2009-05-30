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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkDcmTags_h_included_
#define __cmtkDcmTags_h_included_

#  include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Structure to hold the internal DICOM data dictionary.
 * As the OFFIS library also holds a (complete) data dictionary, this data
 * structure and in fact this whole file may become obsolete one day.
 */
typedef struct
{
  /** Internal reference number as used by StudyInfo.
   *@see StudyInfo
   */
  int id;
  
  /** DICOM group identifier.
   */
  int group;

  /** DICOM element identifier within the given group.
   */
  int elem;

  /** Texutal description of this data element.
   */
  const char *name;
} NamedDcmTag;

/** Table of all supported data elements.
 * So this is our acutal data dictionary. Its end is defined by an entry
 * with "id" -1 and all remaining fields set to zero (NULL).
 */
extern const NamedDcmTag NamedDcmTagTable[];

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDcmTags_h_included_
