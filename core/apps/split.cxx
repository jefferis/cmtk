/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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
#include <System/cmtkExitException.h>

#include <IO/cmtkClassStreamAffineXform.h>
#include <IO/cmtkClassStreamOutput.h>

#include <Base/cmtkUniformVolume.h>
#include <IO/cmtkVolumeIO.h>

std::string InputFilePath;
std::string OutputFilePath;
std::string OutputXformPath;

int Axis = cmtk::AXIS_Z;
int Factor = 2;
bool Padded = false;

int doMain(const int argc, const char *argv[]) {
  try {
    cmtk::CommandLine cl;
    cl.SetProgramInfo(cmtk::CommandLine::PRG_TITLE, "Split images");
    cl.SetProgramInfo(cmtk::CommandLine::PRG_DESCR,
                      "Split volume image into sub-images, i.e., to separate "
                      "interleaved images into passes");

    typedef cmtk::CommandLine::Key Key;

    cl.BeginGroup("Orientation", "Orientation Options");
    cmtk::CommandLine::EnumGroup<int>::SmartPtr interleaveGroup = cl.AddEnum(
        "slice-orientation", &Axis, "Define slice axis for splitting.");
    interleaveGroup->AddSwitch(Key('a', "axial"), (int)cmtk::AXIS_Z,
                               "Interleaved axial images");
    interleaveGroup->AddSwitch(Key('s', "sagittal"), (int)cmtk::AXIS_X,
                               "Interleaved sagittal images");
    interleaveGroup->AddSwitch(Key('c', "coronal"), (int)cmtk::AXIS_Y,
                               "Interleaved coronal images");
    interleaveGroup->AddSwitch(Key('x', "interleave-x"), (int)cmtk::AXIS_X,
                               "Interleaved along x axis");
    interleaveGroup->AddSwitch(Key('y', "interleave-y"), (int)cmtk::AXIS_Y,
                               "Interleaved along y axis");
    interleaveGroup->AddSwitch(Key('z', "interleave-z"), (int)cmtk::AXIS_Z,
                               "Interleaved along z axis");
    cl.EndGroup();

    cl.BeginGroup("Splitting", "Splitting Options");
    cl.AddOption(Key('f', "factor"), &Factor,
                 "Interleave factor. This is the number of subimages "
                 "generated. If this is set to zero, a separate 2D output "
                 "image is generated for each slice in the input "
                 "(in the given slice orientation).");
    cl.EndGroup();

    cl.BeginGroup("Output", "Output Options");
    cl.AddSwitch(Key('p', "padded"), &Padded, true,
                 "Padded output, i.e., fill in removed slices");
    cl.AddOption(Key("output-xform-path"), &OutputXformPath,
                 "Optional path template (fprintf-style) for output affine "
                 "transformation that maps input image coordinates to each "
                 "output image.");
    cl.EndGroup();

    cl.AddParameter(&InputFilePath, "InputImage", "Input image path")
        ->SetProperties(cmtk::CommandLine::PROPS_IMAGE);
    cl.AddParameter(
          &OutputFilePath, "OutputImagePattern",
          "Output image path pattern. Use '%d' to substitute subimage index.")
        ->SetProperties(cmtk::CommandLine::PROPS_IMAGE |
                        cmtk::CommandLine::PROPS_OUTPUT);

    cl.Parse(argc, argv);
  } catch (const cmtk::CommandLine::Exception &e) {
    cmtk::StdErr << e << "\n";
    return 1;
  }

  cmtk::UniformVolume::SmartPtr volume(
      cmtk::VolumeIO::ReadOriented(InputFilePath));

  // if Factor is zero, set to number of slices and generate 2D output images.
  if (Factor == 0) Factor = volume->m_Dims[Axis];

  for (int i = 0; i < Factor; ++i) {
    cmtk::UniformVolume::SmartPtr subvolume(
        Padded ? volume->GetInterleavedPaddedSubVolume(Axis, Factor, i)
               : volume->GetInterleavedSubVolume(Axis, Factor, i));

    char path[PATH_MAX];
    if (snprintf(path, PATH_MAX, OutputFilePath.c_str(), i) > PATH_MAX) {
      cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
    } else {
      cmtk::VolumeIO::Write(*subvolume, path);
    }

    if (!OutputXformPath.empty()) {
      if (snprintf(path, PATH_MAX, OutputXformPath.c_str(), i) > PATH_MAX) {
        cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
      }
      cmtk::AffineXform xform;
      cmtk::Types::Coordinate xlate[3] = {0, 0, 0};
      xlate[Axis] = -i * volume->m_Delta[Axis];

      try {
        xform.SetXlate(xlate);
      } catch (
          const cmtk::AffineXform::MatrixType::SingularMatrixException &ex) {
        cmtk::StdErr
            << "ERROR: singular matrix in cmtk::AffineXform::SetXlate()\n";
        throw cmtk::ExitException(1);
      }

      cmtk::ClassStreamOutput stream(path, cmtk::ClassStreamOutput::MODE_WRITE);
      if (stream.IsValid()) {
        stream << xform;
        stream.Close();
      }
    }
  }

  return 0;
}

#include "cmtkSafeMain"
