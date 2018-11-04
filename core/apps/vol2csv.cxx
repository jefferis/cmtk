/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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
#include <System/cmtkStrUtility.h>

#include <IO/cmtkVolumeIO.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

int doMain(const int argc, const char *argv[]) {
  std::string regionsImagePath;
  std::string labelsFilePath;
  double pscaleFactor = 1;
  std::string pscaleImagePath;
  std::string outputFilePath;
  std::string densityLabels;

  std::vector<std::string> densityImagePaths;

  cmtk::Types::DataItem normalizeDensities = 1.0;

  try {
    cmtk::CommandLine cl;
    cl.SetProgramInfo(cmtk::CommandLine::PRG_TITLE,
                      "Compute regional volumes and write to CSV file.");
    cl.SetProgramInfo(
        cmtk::CommandLine::PRG_DESCR,
        "This tool computes the volumes of regions in a label image. "
        "It optionally accepts density maps (e.g., for different tissues) and "
        "computes and prints the per-region content for each. "
        "Also, the tool can accept an optional 'pixel volume' map to account "
        "for local pixel volume variations, e.g., due to spatial distortion.");

    cl.AddParameter(&regionsImagePath, "RegionsImage",
                    "Image of labeled regions.")
        ->SetProperties(cmtk::CommandLine::PROPS_IMAGE);
    cl.AddParameterVector(
          &densityImagePaths, "DensityImages",
          "List of density images. For each image given here, the total "
          "density per region is computed for each label in the regions image.")
        ->SetProperties(cmtk::CommandLine::PROPS_IMAGE |
                        cmtk::CommandLine::PROPS_OPTIONAL);

    typedef cmtk::CommandLine::Key Key;

    cl.BeginGroup("input", "Input Options");
    cl.AddOption(
        Key("normalize-densities"), &normalizeDensities,
        "Optional normalization factor for density images. Typically, the "
        "values in the density images should be in the range 0..1, but often "
        "such images are scaled to "
        "different ranges to accomodate storage as integers. If, for example, "
        "densities are stored as values 0..255, set this paramater to 255.");
    cl.AddOption(Key("labels-file"), &labelsFilePath,
                 "If provided, this text file contains names for all labels in "
                 "the regions image. These names are then used to label the "
                 "rows of the CSV output.");
    cl.EndGroup();

    cl.BeginGroup("correct", "Correction Options");
    cl.AddOption(
        Key("pixel-scale-factor"), &pscaleFactor,
        "If provided, this global scale factor is applied to all pixel volumes "
        "to compensate for deviations between world and image scale.");
    cl.AddOption(Key("pixel-scale-image"), &pscaleImagePath,
                 "If provided, this volume contains scale factors for the "
                 "volume of each pixel. This is typically the Jacobian "
                 "determinant map of a spatial unwarping deformation.")
        ->SetProperties(cmtk::CommandLine::PROPS_IMAGE);
    cl.EndGroup();

    cl.BeginGroup("output", "Output Options");
    cl.AddOption(Key("density-labels"), &densityLabels,
                 "This option can be used to provide labels for the density "
                 "maps, which are used as column labels in the output. "
                 "Labels must be separated by commas and must not contain any "
                 "unescaped spaces.");
    cl.AddOption(Key('o', "output"), &outputFilePath,
                 "If provided, program output is written to this file. If not "
                 "provided, output is written to the STDOUT stream.");
    cl.EndGroup();

    cl.Parse(argc, argv);
  } catch (const cmtk::CommandLine::Exception &e) {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException(1);
  }

  // read ROI image
  cmtk::UniformVolume::SmartConstPtr regionsImage(
      cmtk::VolumeIO::ReadOriented(regionsImagePath));
  if (!regionsImage) {
    cmtk::StdErr << "ERROR: could not read regions image " << regionsImagePath
                 << "\n";
    throw cmtk::ExitException(1);
  }

  // read optional label name text file
  std::map<size_t, std::string> labelToNameMap;
  if (!labelsFilePath.empty()) {
    std::ifstream labelsFile(labelsFilePath.c_str());
    if (!labelsFile.good()) {
      cmtk::StdErr << "ERROR: could not read label file " << labelsFilePath
                   << "\n";
      throw cmtk::ExitException(1);
    }

    size_t idx;
    std::string name;
    std::string restOfLine;
    while (!labelsFile.eof()) {
      labelsFile >> idx >> name;
      labelToNameMap[idx] = name;
      std::getline(labelsFile, restOfLine);
    }
  }

  // read optional pixel volume scale image
  cmtk::UniformVolume::SmartPtr pscaleImage;
  if (!pscaleImagePath.empty()) {
    pscaleImage = cmtk::VolumeIO::ReadOriented(pscaleImagePath);

    if (!pscaleImage) {
      cmtk::StdErr << "ERROR: could not read pixel volume image "
                   << pscaleImagePath << "\n";
      throw cmtk::ExitException(1);
    }

    if (!regionsImage->GridMatches(*pscaleImage)) {
      cmtk::StdErr << "ERROR: grid of pixel volume image " << pscaleImagePath
                   << " does not match that of the regions image.\n";
      throw cmtk::ExitException(1);
    }
  }

  // read optional density images
  std::vector<cmtk::UniformVolume::SmartConstPtr> densityImages;
  for (size_t idx = 0; idx < densityImagePaths.size(); ++idx) {
    cmtk::UniformVolume::SmartPtr nextImage(
        cmtk::VolumeIO::ReadOriented(densityImagePaths[idx]));
    if (!nextImage) {
      cmtk::StdErr << "ERROR: could not read density image "
                   << densityImagePaths[idx] << "\n";
      throw cmtk::ExitException(1);
    }

    if (!regionsImage->GridMatches(*nextImage)) {
      cmtk::StdErr << "ERROR: grid of density image " << densityImagePaths[idx]
                   << " does not match that of the regions image.\n";
      throw cmtk::ExitException(1);
    }

    if (normalizeDensities != 1.0)
      nextImage->GetData()->Rescale(1.0 / normalizeDensities);

    densityImages.push_back(nextImage);
  }

  // parse optional density labels for output columns
  std::vector<std::string> densityLabelsVector;
  if (densityLabels.empty()) {
    for (size_t midx = 0; midx < densityImages.size(); ++midx) {
      std::ostringstream strm;
      strm << "density" << midx;
      densityLabelsVector.push_back(strm.str());
    }
  } else {
    densityLabelsVector = cmtk::StrSplit(densityLabels, ",");
  }

  // make sure we have exactly one label per column
  if (densityLabelsVector.size() != densityImages.size()) {
    cmtk::StdErr << "ERROR: must provide exactly one density label per density "
                    "image (identified "
                 << densityLabelsVector.size() << " labels for "
                 << densityImages.size() << " images)\n";
    throw cmtk::ExitException();
  }

  // compute pixel volume
  const cmtk::Types::Coordinate pixelVolumeRegionsImage =
      regionsImage->m_Delta.Product();

  // compute number of labels in the ROI image
  const size_t maxLabel = std::max(
      1, static_cast<int>(regionsImage->GetData()->GetRange().m_UpperBound));
  std::vector<cmtk::Types::Coordinate> regionVolumes(1 + maxLabel, 0.0);

  // prepare vector for volume per label and density map
  std::vector<std::vector<cmtk::Types::Coordinate>> regionDensities(
      densityImages.size());
  for (size_t midx = 0; midx < densityImages.size(); ++midx) {
    regionDensities[midx].resize(1 + maxLabel, 0.0);
  }

  // go over all pixels and count/compound volumes
  for (size_t px = 0; px < regionsImage->GetNumberOfPixels(); ++px) {
    const size_t label = std::min<size_t>(
        maxLabel, std::max(0, static_cast<int>(regionsImage->GetDataAt(px))));

    // get size for this pixel and scale, depending on whether we have a
    // per-pixel map or not.
    cmtk::Types::Coordinate pixelVolume =
        pixelVolumeRegionsImage * pscaleFactor;

    if (pscaleImage) pixelVolume *= pscaleImage->GetDataAt(px);

    regionVolumes[label] += pixelVolume;

    for (size_t midx = 0; midx < densityImages.size(); ++midx) {
      regionDensities[midx][label] +=
          pixelVolume * densityImages[midx]->GetDataAt(px);
    }
  }

  // select either output file or standard output
  std::ofstream outputFile;
  std::ostream &output = !outputFilePath.empty()
      ? outputFile.open(outputFilePath.c_str(), std::ios::out),
               outputFile : std::cout;

  // write column labels
  output << "label,volume";
  for (size_t midx = 0; midx < densityImages.size(); ++midx) {
    output << "," << densityLabelsVector[midx];
  }
  output << "\n";

  // write rows with label volumes
  for (size_t label = 0; label <= maxLabel; ++label) {
    if (!labelsFilePath.empty()) {
      std::map<size_t, std::string>::const_iterator it =
          labelToNameMap.find(label);
      if (it == labelToNameMap.end()) continue;

      output << "\"" << it->second << "\"";
    } else {
      output << label;
    }

    output << "," << regionVolumes[label];
    for (size_t midx = 0; midx < densityImages.size(); ++midx) {
      output << "," << regionDensities[midx][label];
    }
    output << "\n";
  }

  return 0;
}

#include "cmtkSafeMain"
