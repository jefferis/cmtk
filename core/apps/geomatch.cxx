/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
//  $Revision: 4818 $
//
//  $LastChangedDate: 2013-09-10 11:28:54 -0700 (Tue, 10 Sep 2013) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkFixedSquareMatrix.h>
#include <Base/cmtkUniformVolume.h>

#include <IO/cmtkVolumeIO.h>

#include <stdio.h>
#include <string>
#include <vector>

cmtk::UniformVolume::SmartPtr readVolume(const std::string &path,
                                         const char *readOrientation) {
  cmtk::UniformVolume::SmartPtr volume;
  try {
    if (readOrientation)
      volume = cmtk::VolumeIO::ReadOriented(path, readOrientation);
    else
      volume = cmtk::VolumeIO::Read(path);
  } catch (...) {
    if (readOrientation)
      volume = cmtk::VolumeIO::ReadGridOriented(path, readOrientation);
    else
      volume = cmtk::VolumeIO::ReadGrid(path);
  }

  if (!volume) {
    cmtk::StdErr << "ERROR: cannot read image from " << path << "\n";
    throw cmtk::ExitException(1);
  }

  return volume;
}

int doMain(int argc, const char *argv[]) {
  const char *readOrientation = NULL;
  std::vector<std::string> imagePaths;

  bool noCheckXforms = false;
  bool noCheckPixels = false;

  double tolerance = 1e-6;
  double toleranceXlate = 1e-3;

  try {
    cmtk::CommandLine cl;
    cl.SetProgramInfo(
        cmtk::CommandLine::PRG_TITLE,
        "Check whether the geometries (e.g., grid dimensions, pixel sizes, "
        "spatial coordinates) or two or more images match.");
    cl.SetProgramInfo(
        cmtk::CommandLine::PRG_DESCR,
        "This tool reads two or more images and tests whether their grid "
        "dimensions, pixel sizes, and image-to-space transformations match. "
        "Optionally, all images are reoriented into standard orientation "
        "before performing the test. If all images match, the tool returns "
        "with exit code 0, otherwise it returns with exit code 2. "
        "In case of an error (e.g., one of the images can not be read), the "
        "exit code is 1.");

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup("Input", "Input Options");
    cl.AddSwitch(Key("read-ras"), &readOrientation, "RAS",
                 "Read all images in RAS orientation");
    cl.EndGroup();

    cl.BeginGroup("Comparison", "Image Comparison Options");
    cl.AddSwitch(Key("no-check-xforms"), &noCheckXforms, true,
                 "Do not check transformation matrices.");
    cl.AddSwitch(Key("no-check-pixels"), &noCheckPixels, true,
                 "Do not check pixelsizes.");
    cl.AddOption(Key("tolerance"), &tolerance,
                 "Numerical tolerance for floating point comparisons of "
                 "transformation matrices.");
    cl.AddOption(Key("tolerance-xlate"), &toleranceXlate,
                 "Numerical tolerance for floating point comparisons of the "
                 "translational components of the transformation matrices.");
    cl.EndGroup();

    cl.AddParameterVector(&imagePaths, "ImagePaths", "List of image files.");

    cl.Parse(argc, const_cast<const char **>(argv));
  } catch (const cmtk::CommandLine::Exception &ex) {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException(1);
  }

  // Check if we have enough parameters
  if (imagePaths.size() < 2) {
    cmtk::StdErr << "ERROR: need at least two image paths.\n";
    throw cmtk::ExitException(1);
  }

  cmtk::UniformVolume::SmartConstPtr firstVolume =
      readVolume(imagePaths[0], readOrientation);
  const cmtk::AffineXform::MatrixType firstVolumeMatrix =
      firstVolume->GetImageToPhysicalMatrix();

  for (size_t i = 1; i < imagePaths.size(); ++i) {
    cmtk::UniformVolume::SmartConstPtr nextVolume =
        readVolume(imagePaths[i], readOrientation);

    // First and always, use the DataGrid member function to check grid
    // dimensions
    if (!firstVolume->DataGrid::GridMatches(*nextVolume)) {
      cmtk::DebugOutput(1) << "MISMATCH: grid dimensions\n";
      return 2;
    }

    // Check pixels - use default (UniformVolume) member function
    if (!noCheckPixels) {
      if (!firstVolume->GridMatches(*nextVolume)) {
        cmtk::DebugOutput(1) << "MISMATCH: pixel size\n";
        return 2;
      }
    }

    if (!noCheckXforms) {
      // Check rotational part of image matrices
      const cmtk::AffineXform::MatrixType nextVolumeMatrix =
          nextVolume->GetImageToPhysicalMatrix();
      if ((firstVolumeMatrix.GetTopLeft3x3() - nextVolumeMatrix.GetTopLeft3x3())
              .FrobeniusNorm() > tolerance) {
        cmtk::DebugOutput(1)
            << "MISMATCH: image-to-space matrix (rotational part)\n";
        return 2;
      }

      // Check translational part of image matrices
      if ((firstVolumeMatrix.GetRowVector(3) - nextVolumeMatrix.GetRowVector(3))
              .MaxAbsValue() > toleranceXlate) {
        cmtk::DebugOutput(1)
            << "MISMATCH: image-to-space matrix (translational part)\n";
        return 2;
      }
    }
  }

  return 0;
}

#include "cmtkSafeMain"
