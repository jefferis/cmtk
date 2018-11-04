/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
//
//  Copyright 2015 Google, Inc.
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

#include "cmtkGroupwiseRegistrationOutput.h"

#include <System/cmtkDebugOutput.h>

#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkClassStreamStudyList.h>
#include <IO/cmtkGroupwiseRegistrationFunctionalIO.h>
#include <IO/cmtkStudyList.h>
#include <IO/cmtkVolumeIO.h>

#include <limits.h>

#include <algorithm>
#include <sstream>
#include <vector>

namespace cmtk {

/** \addtogroup Registration */
//@{

bool GroupwiseRegistrationOutput::WriteGroupwiseArchive(
    const char *path) const {
  // create class stream archive.
  if (path) {
    ClassStreamOutput stream;

    if (this->m_OutputRootDirectory) {
      char completePath[PATH_MAX];
      snprintf(completePath, sizeof(completePath), "%s%c%s",
               this->m_OutputRootDirectory, (int)CMTK_PATH_SEPARATOR, path);
      stream.Open(completePath, ClassStreamOutput::MODE_WRITE_ZLIB);
    } else
      stream.Open(path, ClassStreamOutput::MODE_WRITE_ZLIB);

    if (!stream.IsValid()) return false;
    stream << *this->m_Functional;
    stream.Close();
  }

  return true;
}

bool GroupwiseRegistrationOutput::WriteXformsSeparateArchives(
    const std::string &path, const std::string &templatePath) {
  if (!path.empty()) {
    for (size_t img = 0; img < this->m_Functional->GetNumberOfTargetImages();
         ++img) {
      StudyList slist;
      Study::SmartPtr refstudy;
      if (this->m_OutputRootDirectory && !this->m_ExistingTemplatePath) {
        refstudy = slist.AddStudy(std::string(this->m_OutputRootDirectory) +
                                  CMTK_PATH_SEPARATOR + templatePath);
      } else {
        refstudy = slist.AddStudy(templatePath);
      }

      const UniformVolume *image =
          this->m_Functional->GetOriginalTargetImage(img);
      Study::SmartPtr imgstudy =
          slist.AddStudy(image->GetMetaInfo(META_FS_PATH).c_str());

      WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom(
          this->m_Functional->GetGenericXformByIndex(img));
      if (warpXform) {
        AffineXform::SmartPtr affineXform(warpXform->GetInitialAffineXform());
        slist.AddXform(refstudy, imgstudy, affineXform, warpXform);
      } else {
        AffineXform::SmartPtr affineXform =
            AffineXform::SmartPtr::DynamicCastFrom(
                this->m_Functional->GetGenericXformByIndex(img));
        slist.AddXform(refstudy, imgstudy, affineXform);
      }

      std::ostringstream fullPath;
      if (this->m_OutputRootDirectory) {
        fullPath << this->m_OutputRootDirectory << CMTK_PATH_SEPARATOR;
      }

      fullPath << path << CMTK_PATH_SEPARATOR << "target-";
      fullPath.fill('0');
      fullPath.width(3);
      fullPath << img << ".list";

      ClassStreamStudyList::Write(fullPath.str(), &slist);
    }
  }

  return true;
}

bool GroupwiseRegistrationOutput::WriteAverageImage(
    const char *path, const cmtk::Interpolators::InterpolationEnum interp,
    const cmtk::ScalarDataType pixelType, const bool useTemplateData) {
  // reformat output and generate average images
  if (path) {
    UniformVolume::SmartPtr templateGrid =
        this->m_Functional->GetTemplateGrid();
    const size_t numberOfPixels = templateGrid->GetNumberOfPixels();

    TypedArray::SmartPtr average(TypedArray::Create(pixelType, numberOfPixels));
    float *averagePtr = static_cast<float *>(average->GetDataPtr());

    std::vector<byte> count;

    if (useTemplateData) {
      if (!templateGrid->GetData()) {
        UniformVolume::SmartPtr readImage(
            VolumeIO::ReadOriented(templateGrid->GetMetaInfo(META_FS_PATH)));
        templateGrid->SetData(readImage->GetData());
      }

      for (size_t px = 0; px < numberOfPixels; ++px) {
        averagePtr[px] = static_cast<float>(templateGrid->GetDataAt(px));
      }
      count.resize(numberOfPixels, 1);
    } else {
      average->Fill(0);
      count.resize(numberOfPixels, 0);
    }

    DebugOutput(1) << "Reformating output images\n";

    const size_t idxFrom = 0;
    const size_t idxSkip = 1;
    for (size_t idx = idxFrom;
         idx < this->m_Functional->GetNumberOfTargetImages(); idx += idxSkip) {
      UniformVolume::SmartPtr floatingVolume =
          this->m_Functional->GetOriginalTargetImage(idx);
      if (!floatingVolume->GetData())
        floatingVolume = UniformVolume::SmartPtr(
            VolumeIO::ReadOriented(floatingVolume->GetMetaInfo(META_FS_PATH)));

      cmtk::ReformatVolume reformat;
      reformat.SetReferenceVolume(templateGrid);
      reformat.SetFloatingVolume(floatingVolume);
      reformat.SetInterpolation(interp);

      AffineXform::SmartPtr affineXform =
          AffineXform::SmartPtr::DynamicCastFrom(
              this->m_Functional->GetGenericXformByIndex(idx));
      if (affineXform) {
        reformat.SetAffineXform(affineXform);
      }

      WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom(
          this->m_Functional->GetGenericXformByIndex(idx));
      if (warpXform) reformat.SetWarpXform(warpXform);

      UniformVolume::SmartPtr ref(reformat.PlainReformat());
      const TypedArray *data = ref->GetData();
#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(numberOfPixels); ++i) {
        Types::DataItem v;
        if (data->Get(v, i)) {
          averagePtr[i] += static_cast<float>(v);
          ++count[i];
        }
      }
    }

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(numberOfPixels); ++i) {
      if (count[i])
        averagePtr[i] /= count[i];
      else
        average->SetPaddingAt(i);
    }
    templateGrid->SetData(average);

    if (this->m_OutputRootDirectory) {
      char fullPath[PATH_MAX];
      snprintf(fullPath, sizeof(fullPath), "%s%c%s",
               this->m_OutputRootDirectory, CMTK_PATH_SEPARATOR, path);
      VolumeIO::Write(*templateGrid, fullPath);
    } else {
      VolumeIO::Write(*templateGrid, path);
    }
  }

  return 0;
}

}  // namespace cmtk
