/*
//
//  Copyright 2004-2012 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkVolumeFromStudy.h"

#include <System/cmtkCompressedStream.h>
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkMountPoints.h>
#include <System/cmtkProgress.h>

#include <IO/cmtkDICOM.h>
#include <IO/cmtkVolumeFromFile.h>
#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkLandmarkList.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkVolume.h>

#include <limits.h>
#include <string.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#include <memory>

namespace cmtk {

/** \addtogroup IO */
//@{

const UniformVolume::SmartPtr VolumeFromStudy::Read(
    const Study *study, const Types::Coordinate tolerance) {
  if (!study) return UniformVolume::SmartPtr(NULL);

  const StudyImageSet *studyImageSet =
      dynamic_cast<const StudyImageSet *>(study);
  if (studyImageSet) {
    UniformVolume::SmartPtr volume =
        VolumeFromStudy(tolerance).AssembleVolume(studyImageSet);
    if (!volume)
      StdErr << "ERROR: volume assembly failed in directory "
             << studyImageSet->GetImageDirectory() << "\n";
    return volume;
  } else
    return VolumeIO::Read(study->GetFileSystemPath());
}

const UniformVolume::SmartPtr VolumeFromStudy::AssembleVolume(
    const StudyImageSet *study) {
  UniformVolume::SmartPtr Result(NULL);
  const std::string imageDir =
      MountPoints::Translate(study->GetImageDirectory());

  try {
    DebugOutput(2) << "Reading images from path " << imageDir << "\n";

    Progress::Begin(0, study->size(), 1, "Volume image assembly");

    unsigned int nextPlane = 0;
    StudyImageSet::const_iterator it = study->begin();
    while (it != study->end()) {
      DebugOutput(2) << "\r" << *it;

      char fullpath[PATH_MAX];
      snprintf(fullpath, sizeof(fullpath), "%s%c%s", imageDir.c_str(),
               (int)CMTK_PATH_SEPARATOR, it->c_str());

      ScalarImage::SmartPtr image =
          ScalarImage::SmartPtr(DICOM::Read(fullpath));

      // TODO: when returning NULL here, we also should tell
      // VolumeFromSlices that we give up, so it can free its
      // temporary storage.
      if (!image) return UniformVolume::SmartPtr(NULL);

      if (!nextPlane) {
        // special treatment for first image in sequence.
        if (study->GetMultiFile())
          InitSequence(image, study->size());
        else
          InitSequence(image, study->m_Dims[AXIS_Z]);
      }

      const char *error = FillPlane(nextPlane, image);

      Progress::SetProgress(nextPlane);

      if (error) {
        StdErr.printf("ERROR: %s: %s\n", fullpath, error);
        return UniformVolume::SmartPtr(NULL);
      }

      ++it;
    }
    Progress::Done();

    Result = this->FinishVolume();

    // If seomthing went wrong constructing the volume, return the NULL pointer
    // that we just got
    if (!Result) return Result;

    TypedArray::SmartPtr data = Result->GetData();
    if (data) {
      if (study->GetPadding() && !data->GetPaddingFlag()) {
        data->SetPaddingValue(study->GetPaddingValue());
      }
    }
  } catch (...) {
  }

  return Result;
}

}  // namespace cmtk
