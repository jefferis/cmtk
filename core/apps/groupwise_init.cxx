/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>
#include <System/cmtkTimers.h>

#include <Base/cmtkFilterVolume.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkUniformVolume.h>

#include <IO/cmtkVolumeIO.h>

#include <Registration/cmtkGroupwiseRegistrationFunctionalAffineInitializer.h>
#include <Registration/cmtkGroupwiseRegistrationFunctionalBase.h>
#include <Registration/cmtkGroupwiseRegistrationOutput.h>

#include <vector>

const char *PreDefinedTemplatePath = NULL;
cmtk::UniformVolume::SmartPtr PreDefinedTemplate;

const char *OutputRootDirectory = NULL;
const char *OutputArchive = "groupwise.xforms";
const char *OutputStudyListGroup = "groupwise.list";
const char *OutputStudyListIndividual = "pairs";
const char *AverageImagePath = "average.nii";
cmtk::Interpolators::InterpolationEnum AverageImageInterpolation =
    cmtk::Interpolators::LINEAR;

bool AlignCentersOfMass = false;
bool InitScales = false;
bool CenterTemplate = false;

// this vector holds all target image filenames
std::vector<const char *> fileNameList;

// this vector holds the original (not downsampled) images.
std::vector<cmtk::UniformVolume::SmartPtr> imageListOriginal;

cmtk::ScalarDataType DataType = cmtk::TYPE_FLOAT;

int doMain(int argc, const char *argv[]) {
  try {
    cmtk::CommandLine cl;
    cl.SetProgramInfo(cmtk::CommandLine::PRG_TITLE,
                      "Affine initialization for groupwise registration");
    cl.SetProgramInfo(
        cmtk::CommandLine::PRG_DESCR,
        "Compute initial affine alignment for a group of input images, which "
        "can be used as an input for groupwise registration");
    cl.SetProgramInfo(cmtk::CommandLine::PRG_SYNTX,
                      "groupwise_init [options] image0 [image1 ...]");
    cl.SetProgramInfo(cmtk::CommandLine::PRG_CATEG, "CMTK.Image Registration");

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption(Key('t', "template"), &PreDefinedTemplatePath,
                 "Input filename for pre-defined template image.");

    cl.AddOption(Key('O', "output-root"), &OutputRootDirectory,
                 "Root directory for all output files.");
    cl.AddOption(Key('o', "output"), &OutputArchive,
                 "Output filename for groupwise registration archive.");
    cl.AddOption(Key("output-average"), &AverageImagePath,
                 "Output filename for registered average image.");
    cl.AddSwitch(
        Key("average-cubic"), &AverageImageInterpolation,
        cmtk::Interpolators::CUBIC,
        "Use cubic (rather than linear) interpolation for average image.");
    cl.AddSwitch(Key("no-output-average"), &AverageImagePath,
                 (const char *)NULL, "Do not write average image.");

    cl.AddSwitch(Key("align-centers-of-mass"), &AlignCentersOfMass, true,
                 "Initially align centers of mass rather than centers of "
                 "bounding boxes.");
    cl.AddSwitch(Key("init-scales"), &InitScales, true,
                 "Initialize scale factors using first-order moments");
    cl.AddSwitch(Key("center-template"), &CenterTemplate, true,
                 "Center aligned images in template grid field of view.");

    cmtk::CommandLine::EnumGroup<cmtk::ScalarDataType>::SmartPtr typeGroup =
        cl.AddEnum("outputtype", &DataType,
                   "Scalar data type for the output average image.");
    typeGroup->AddSwitch(Key("char"), cmtk::TYPE_CHAR, "8 bits, signed");
    typeGroup->AddSwitch(Key("byte"), cmtk::TYPE_BYTE, "8 bits, unsigned");
    typeGroup->AddSwitch(Key("short"), cmtk::TYPE_SHORT, "16 bits, signed");
    typeGroup->AddSwitch(Key("ushort"), cmtk::TYPE_USHORT, "16 bits, unsigned");
    typeGroup->AddSwitch(Key("int"), cmtk::TYPE_INT, "32 bits signed");
    typeGroup->AddSwitch(Key("uint"), cmtk::TYPE_UINT, "32 bits unsigned");
    typeGroup->AddSwitch(Key("float"), cmtk::TYPE_FLOAT,
                         "32 bits floating point");
    typeGroup->AddSwitch(Key("double"), cmtk::TYPE_DOUBLE,
                         "64 bits floating point\n");

    cl.Parse(argc, const_cast<const char **>(argv));

    const char *next = cl.GetNext();
    while (next) {
      fileNameList.push_back(next);
      next = cl.GetNextOptional();
    }
  } catch (const cmtk::CommandLine::Exception &e) {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException(1);
  }

  // Make sure we don't exceed maximum number of supported images. This is due
  // to using, for example, "byte" for pixelwise image counts in averaging.
  if (fileNameList.size() > 255) {
    cmtk::StdErr << "ERROR: no more than 255 images are supported.\n";
    throw cmtk::ExitException(1);
  }

  cmtk::GroupwiseRegistrationFunctionalBase::SmartPtr initializer(
      new cmtk::GroupwiseRegistrationFunctionalBase);

  int idx = 0;
  for (std::vector<const char *>::const_iterator fnIt = fileNameList.begin();
       fnIt != fileNameList.end(); ++fnIt, ++idx) {
    cmtk::UniformVolume::SmartPtr nextImage;
    cmtk::UniformVolume::SmartPtr image(cmtk::VolumeIO::ReadOriented(*fnIt));
    if (!image || !image->GetData()) {
      cmtk::StdErr << "ERROR: Could not read image " << *fnIt << "\n";
      throw cmtk::ExitException(1);
    }
    nextImage = image;
    imageListOriginal.push_back(nextImage);
  }

  initializer->SetTargetImages(imageListOriginal);
  if (PreDefinedTemplatePath) {
    PreDefinedTemplate = cmtk::UniformVolume::SmartPtr(
        cmtk::VolumeIO::ReadGridOriented(PreDefinedTemplatePath));
  }

  if (PreDefinedTemplate)
    initializer->SetTemplateGrid(PreDefinedTemplate);
  else
    initializer->CreateTemplateGridFromTargets(imageListOriginal);

  cmtk::UniformVolume::SmartPtr templateGrid = initializer->GetTemplateGrid();

  cmtk::DebugOutput(1).GetStream().printf(
      "Template grid is %d x %d x %d pixels of size %f x %f x %f\n",
      templateGrid->m_Dims[0], templateGrid->m_Dims[1], templateGrid->m_Dims[2],
      templateGrid->m_Delta[0], templateGrid->m_Delta[1],
      templateGrid->m_Delta[2]);

  cmtk::GroupwiseRegistrationFunctionalAffineInitializer::InitializeXforms(
      *initializer, true /*alignCenters*/, AlignCentersOfMass, InitScales,
      CenterTemplate);

  cmtk::GroupwiseRegistrationOutput output;
  if (PreDefinedTemplatePath) {
    output.SetExistingTemplatePath(true);
  } else {
    PreDefinedTemplatePath = AverageImagePath;
  }

  output.SetFunctional(initializer);
  output.SetOutputRootDirectory(OutputRootDirectory);
  output.WriteGroupwiseArchive(OutputArchive);
  output.WriteXformsSeparateArchives(OutputStudyListIndividual,
                                     PreDefinedTemplatePath);
  output.WriteAverageImage(AverageImagePath, AverageImageInterpolation);

  return 0;
}

#include "cmtkSafeMain"
