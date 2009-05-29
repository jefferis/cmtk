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

#include <cmtkDcmTags.h>
#include <cmtkStudyInfo.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

const NamedDcmTag NamedDcmTagTable[] = 
{ { INFO_PATNAME,         0x0010, 0x0010,
    "Patient Name" },
  { INFO_PATID,           0x0010, 0x0020,
    "Patient ID" },
  { INFO_PATBIRTHD,       0x0010, 0x0030,
    "Patient Birtdate" },
  { INFO_PATSEX,          0x0010, 0x0040,
    "Patient Sex"},
  { INFO_PATMAIDNAME,     0x0010, 0x1005,
    "Patient Maidenname" },
  { INFO_PATAGE,          0x0010, 0x1010,
    "Patient Age" },
  { INFO_PATSIZE,         0x0010, 0x1020,
    "Patient Size" },
  { INFO_PATWEIGHT,       0x0010, 0x1030,
    "Patient Weight" },
  { INFO_PATADDRESS,      0x0010, 0x1040,
    "Patient Address" },
  { INFO_PATPLANID,       0x0010, 0x1050,
    "Insurance Plan ID" },
  { INFO_PATMOTHMAIDNAME, 0x0010, 0x1060,
    "Patient's Mother Maidenname" },
  { INFO_ACQSTUDYDATE,    0x0008, 0x0020,
    "Study Date" },
  { INFO_ACQSERIESDATE,   0x0008, 0x0021,
    "Series Date" },
  { INFO_ACQSTUDYTIME,    0x0008, 0x0030,
    "Study Time" },
  { INFO_ACQMODALITY,     0x0008, 0x0060,
    "Modality" },
  { INFO_ACQMANUFACT,     0x0008, 0x0070,
    "Manufacturer" },
  { INFO_ACQINSTID,       0x0008, 0x0080,
    "Institute ID" },
  { INFO_ACQREFPHYS,      0x0008, 0x0090,
    "Referring Physician" },
  { INFO_ACQSTATID,       0x0008, 0x1010,
    "Station ID" },
  { INFO_ACQINSTDEP,      0x0008, 0x1040,
    "Institute Department" },
  { INFO_ACQPROCDESCR,    0x0008, 0x1030,
    "Procedure Description" },
  { INFO_ACQATTPHYS,      0x0008, 0x1050,
    "Attending Physician" },
  { INFO_ACQRADIOL,       0x0008, 0x1060,
    "Physicians Reading Study" },
  { INFO_ACQOPERID,       0x0008, 0x1070,
    "Operator's ID" },
  { INFO_ACQMANUMODEL,    0x0008, 0x1090,
    "Manufacturer Model" },
  { INFO_ACQSCANSEQ,      0x0018, 0x0020,
    "Scanning Sequence" }, 
  { INFO_ACQSLICETHICK,   0x0018, 0x0050,
    "Slice Thickness" }, 
  { INFO_ACQSLICESPACING, 0x0018, 0x0088,
    "Spacing Between Slices" }, 
  { INFO_ACQGANTILT,      0x0018, 0x1120,
    "Gantry Detector Tilt" },
  { INFO_ACQPATPOSIT,     0x0018, 0x5100,
    "Patient Position" },
  { INFO_ACQNUMBER,       0x0020, 0x0012,
    "Acquisition Number" },
  { INFO_RELPATORIENT,    0x0020, 0x0020,
    "Patient Orientation" },
  { INFO_RELIMGPOSITION,  0x0020, 0x0030,
    "Image Position" },
  { INFO_RELIMGPOSITPAT,  0x0020, 0x0032,
    "Image Position Patient" },
  { INFO_RELIMGORIENT,    0x0020, 0x0035,
    "Image Orientation" },
  { INFO_RELIMGORIENTPAT, 0x0020, 0x0037,
    "Image Orientation Patient" },
  { INFO_RELIMGLOCATION,  0x0020, 0x0050,
    "Image Location" },
  { INFO_RELIMGLATERAL,   0x0020, 0x0060,
    "Image Laterality" },
  { -1,                   0x0000, 0x0000, 
    NULL }
};

} // namespace cmtk
