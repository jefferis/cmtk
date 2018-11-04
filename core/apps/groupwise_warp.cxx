/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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
//  $Revision: 2697 $
//
//  $LastChangedDate: 2011-01-10 16:40:53 -0800 (Mon, 10 Jan 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>
#include <System/cmtkTimers.h>

#include <Base/cmtkFilterVolume.h>
#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkWarpXform.h>

#include <Registration/cmtkBestDirectionOptimizer.h>
#include <Registration/cmtkGroupwiseRegistrationOutput.h>
#include <Registration/cmtkSplineWarpCongealingFunctional.h>
#include <Registration/cmtkSplineWarpGroupwiseRegistrationRMIFunctional.h>

#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkVolumeIO.h>

#include <vector>

bool OptimizeRMI = false;

int DownsampleFrom = 4;
int DownsampleTo = 1;

const char *PreDefinedTemplatePath = NULL;
bool UseTemplateData = false;

const char *DisableControlPointsMaskPath = NULL;

int RefineTransformationGrid = 0;
float JacobianConstraintWeight = 0.0;
float BendingEnergyWeight = 0.0;

cmtk::Types::Coordinate GridSpacing = 40.0;
bool GridSpacingExact = true;
bool ForceZeroSum = false;
bool ForceZeroSumNoAffine = false;
unsigned int ForceZeroSumFirstN = 0;
unsigned int NormalGroupFirstN = 0;

float SamplingDensity = -1.0;

bool DeactivateUninformative = true;
float PartialGradientThreshold = 0.0;

bool UseNumberOfHistogramBins = false;
unsigned int NumberOfHistogramBins = 0;

bool CropImageHistograms = false;
bool HistogramMatching = false;
bool RepeatHistogramMatching = false;

bool UseSmoothSigmaFactorPixel = false;
cmtk::Types::Coordinate SmoothSigmaFactorPixel = 0.0;

bool UseSmoothSigmaFactorControlPointSpacing = false;
cmtk::Types::Coordinate SmoothSigmaFactorControlPointSpacing = 0.0;

cmtk::Types::Coordinate Accuracy = 0.01;
cmtk::Types::Coordinate Exploration = 0.25;
cmtk::Types::Coordinate OptimizerStepFactor = 0.5;
cmtk::Optimizer::ReturnType OptimizerDeltaFThreshold = 0;
bool OptimizerAggressive = true;
int OptimizerRepeatLevel = 2;

bool DisableOptimization = false;

const char *AffineGroupRegistration = NULL;

const char *OutputRootDirectory = NULL;
const char *OutputArchive = "groupwise.xforms";
const char *OutputStudyListGroup = "groupwise.list";
const char *OutputStudyListIndividual = "pairs";
const char *AverageImagePath = "average.nii";
cmtk::Interpolators::InterpolationEnum AverageImageInterpolation =
    cmtk::Interpolators::LINEAR;

byte UserBackgroundValue = 0;
bool UserBackgroundFlag = false;

int doMain(int argc, const char *argv[]) {
  try {
    cmtk::CommandLine cl;
    cl.SetProgramInfo(cmtk::CommandLine::PRG_TITLE,
                      "Nonrigid population registration");
    cl.SetProgramInfo(
        cmtk::CommandLine::PRG_DESCR,
        "This tool nonrigidly registers a population of input images "
        "simultaneously, without a template, using either the 'congealing' "
        "algorithm or a groupwise similarity measure based on "
        "a continuous approximation of mutual information ('RMI').");
    cl.SetProgramInfo(cmtk::CommandLine::PRG_SYNTX,
                      "groupwise_warp [options] affineGroupRegistration");
    cl.SetProgramInfo(cmtk::CommandLine::PRG_CATEG, "CMTK.Image Registration");

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup("Metric", "Registration Metric Options");
    cl.AddSwitch(Key("rmi"), &OptimizeRMI, true,
                 "Use the RMI (a.k.a. regional mutual information) metric to "
                 "drive the registration).");
    cl.AddSwitch(Key("congeal"), &OptimizeRMI, false,
                 "Use the congealing algorithm using pixelwise stack entropies "
                 "to drive the registration.");
    cl.EndGroup();

    cl.BeginGroup("Template", "Template Image Options");
    cl.AddOption(Key('t', "template"), &PreDefinedTemplatePath,
                 "Override template image with given file.");
    cl.AddOption(
        Key('T', "template-with-data"), &PreDefinedTemplatePath,
        "Use user-supplied template images's pixel data in registration",
        &UseTemplateData);
    cl.EndGroup();

    cl.BeginGroup("Multuresolution", "Multuresolution Parameters");
    cl.AddOption(Key('d', "downsample-from"), &DownsampleFrom,
                 "Initial downsampling factor [4].");
    cl.AddOption(Key('D', "downsample-to"), &DownsampleTo,
                 "Final downsampling factor [1].");
    cl.AddOption(Key('s', "sampling-density"), &SamplingDensity,
                 "Probabilistic sampling density [default: off].");
    cl.EndGroup();

    cl.BeginGroup("Output", "Output Options");
    cl.AddOption(Key('O', "output-root"), &OutputRootDirectory,
                 "Root directory for all output files.");
    cl.AddOption(Key('o', "output"), &OutputArchive,
                 "Output filename for groupwise registration archive.");
    cl.AddOption(Key("output-average"), &AverageImagePath,
                 "Output filename for registered average image.");
    cl.AddSwitch(Key("no-output-average"), &AverageImagePath,
                 (const char *)NULL, "Do not write average image.");
    cl.AddSwitch(Key("average-cubic"), &AverageImageInterpolation,
                 cmtk::Interpolators::CUBIC,
                 "Use cubic interpolation for average image (default: linear)");
    cl.EndGroup();

    cl.BeginGroup("Image", "Image Options and Operations");
    cl.AddOption(Key('B', "force-background"), &UserBackgroundValue,
                 "Force background pixels (outside FOV) to given (bin) value.",
                 &UserBackgroundFlag);
    cl.AddOption(Key('H', "histogram-bins"), &NumberOfHistogramBins,
                 "Set number of histogram bins for entropy evaluation.",
                 &UseNumberOfHistogramBins);
    cl.AddSwitch(Key("crop-histograms"), &CropImageHistograms, true,
                 "Crop image histograms to make better use of histogram bins.");
    cl.AddSwitch(Key("match-histograms"), &HistogramMatching, true,
                 "Match all image histograms to template data (or first image, "
                 "if no template image is given)");
    cl.AddSwitch(Key("repeat-match-histograms"), &RepeatHistogramMatching, true,
                 "Frequetly repeat histogram-based intensity matching to "
                 "account for changing volume proportions.");
    cl.AddOption(Key("smooth-pixels"), &SmoothSigmaFactorPixel,
                 "Sigma of Gaussian smoothing kernel in multiples of template "
                 "image pixel size",
                 &UseSmoothSigmaFactorPixel);
    cl.AddOption(Key("smooth-cps"), &SmoothSigmaFactorControlPointSpacing,
                 "Sigma of Gaussian smoothing kernel in multiples of control "
                 "point delta",
                 &UseSmoothSigmaFactorControlPointSpacing);
    cl.EndGroup();

    cl.BeginGroup("Transformation", "Transformation Parameters");
    cl.AddOption(Key("grid-spacing"), &GridSpacing,
                 "Control point grid spacing.");
    cl.AddSwitch(Key("grid-spacing-fit"), &GridSpacingExact, false,
                 "Use grid spacing that fits volume FOV");
    cl.AddOption(Key('r', "refine-grid"), &RefineTransformationGrid,
                 "Number of times to refine transformation grid [default: 0].");
    cl.EndGroup();

    cl.BeginGroup("Constraints", "Transformation Constraint Options");
    cl.AddSwitch(Key('z', "zero-sum"), &ForceZeroSum, true,
                 "Enforce zero-sum computation.");
    cl.AddSwitch(Key("zero-sum-no-affine"), &ForceZeroSumNoAffine, true,
                 "Enforce zero-sum computation EXCLUDING affine components.");
    cl.AddOption(Key('N', "normal-group-first-n"), &NormalGroupFirstN,
                 "First N images are from the normal group and should be "
                 "registered unbiased.");
    cl.AddOption(Key('Z', "zero-sum-first-n"), &ForceZeroSumFirstN,
                 "Enforce zero-sum computation for first N images.",
                 &ForceZeroSum);
    cl.AddOption(
        Key('J', "jacobian-weight"), &JacobianConstraintWeight,
        "Weight for Jacobian volume preservation constraint [default: off]");
    cl.AddOption(Key('E', "bending-energy-weight"), &BendingEnergyWeight,
                 "Weight for grid bending energy regularization constraint "
                 "[default: off]");
    cl.EndGroup();

    cl.BeginGroup("Optimization", "Optimization Parameters");
    cl.AddOption(Key('e', "exploration"), &Exploration,
                 "Exploration of optimization in pixels");
    cl.AddOption(Key('a', "accuracy"), &Accuracy,
                 "Accuracy of optimization in pixels");
    cl.AddOption(Key('S', "step-factor"), &OptimizerStepFactor,
                 "Step factor for successive optimization passes");
    cl.AddOption(Key("delta-f-threshold"), &OptimizerDeltaFThreshold,
                 "Optional threshold to terminate optimization (level) if "
                 "relative change of target function drops below this value.");
    cl.AddOption(Key('p', "partial-gradient-thresh"), &PartialGradientThreshold,
                 "Threshold factor for partial gradient zeroing [<0 turn off]");
    cl.AddSwitch(Key("activate-uninformative"), &DeactivateUninformative, false,
                 "Activate potentially uninformative control points");
    cl.AddOption(
        Key("disable-cp-mask"), &DisableControlPointsMaskPath,
        "Path to mask image (matching template grid) defining areas in which "
        "control points should be disabled. "
        "This guarantees that mask foreground areas remain undeformed.");
    cl.AddSwitch(Key("disable-optimization"), &DisableOptimization, true,
                 "Disable optimization and output initial configuration.");
    cl.EndGroup();

    cl.Parse(argc, const_cast<const char **>(argv));

    AffineGroupRegistration = cl.GetNext();
  } catch (const cmtk::CommandLine::Exception &e) {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException(1);
  }

  cmtk::GroupwiseRegistrationFunctionalXformTemplate<
      cmtk::SplineWarpXform>::SmartPtr functional;
  if (OptimizeRMI)
    functional = cmtk::SplineWarpGroupwiseRegistrationRMIFunctional::SmartPtr(
        new cmtk::SplineWarpGroupwiseRegistrationRMIFunctional);
  else
    functional = cmtk::SplineWarpCongealingFunctional::SmartPtr(
        new cmtk::SplineWarpCongealingFunctional);

  functional->SetFreeAndRereadImages(
      !HistogramMatching);  // we cannot unload the original images if we still
                            // need them for histogram matching
  functional->SetForceZeroSumFirstN(ForceZeroSumFirstN);
  functional->SetForceZeroSum(ForceZeroSum || ForceZeroSumNoAffine);
  functional->SetForceZeroSumNoAffine(ForceZeroSumNoAffine);
  functional->SetCropImageHistograms(CropImageHistograms);
  if (UserBackgroundFlag)
    functional->SetUserBackgroundValue(UserBackgroundValue);

  if (UseNumberOfHistogramBins) {
    // must be done IMMEDIATELY after creating the functional!
    // otherwise, scaling and conversion of input images goes
    // WRONG!
    functional->SetNumberOfHistogramBins(NumberOfHistogramBins);
  }

  // disable parameters with less than 1% of maximum contribution
  functional->SetPartialGradientMode((PartialGradientThreshold > 0),
                                     PartialGradientThreshold);
  functional->SetDeactivateUninformativeMode(DeactivateUninformative);

  cmtk::ClassStreamInput stream(AffineGroupRegistration);
  if (!stream.IsValid()) {
    cmtk::StdErr << "Input archive " << AffineGroupRegistration
                 << " could not be opened for reading.\n";
    throw cmtk::ExitException(1);
  }
  stream >> *functional;
  stream.Close();

  // Make sure we don't exceed maximum number of supported images. This is due
  // to using, for example, "byte" for pixelwise image counts in averaging.
  if (functional->GetNumberOfTargetImages() > 255) {
    cmtk::StdErr << "ERROR: no more than 255 images are supported.\n";
    throw cmtk::ExitException(1);
  }

  if (PreDefinedTemplatePath) {
    cmtk::UniformVolume::SmartPtr templateImage;
    if (UseTemplateData) {
      templateImage = cmtk::UniformVolume::SmartPtr(
          cmtk::VolumeIO::ReadOriented(PreDefinedTemplatePath));
    } else {
      templateImage = cmtk::UniformVolume::SmartPtr(
          cmtk::VolumeIO::ReadGridOriented(PreDefinedTemplatePath));
    }

    if (!templateImage) {
      cmtk::StdErr << "ERROR: could not read template grid/image "
                   << PreDefinedTemplatePath << "\n";
      throw cmtk::ExitException(2);
    }

    functional->SetTemplateGrid(templateImage, 0 /*downsample*/,
                                false /*useTemplateData: set this later*/);
  }

  // this vector holds the original (not downsampled) images.
  std::vector<cmtk::UniformVolume::SmartPtr> imageListOriginal;
  functional->GetOriginalTargetImages(imageListOriginal);

  cmtk::UniformVolume::SmartPtr originalTemplateGrid =
      functional->GetTemplateGrid();

  if (HistogramMatching) {
    const cmtk::TypedArray *referenceDataForHistogramMatching = NULL;

    bool useTemplateForHistogramMatching = true;
    if (originalTemplateGrid && UseTemplateData) {
      referenceDataForHistogramMatching = originalTemplateGrid->GetData();
    }
    if (!referenceDataForHistogramMatching) {
      useTemplateForHistogramMatching = false;
      referenceDataForHistogramMatching = imageListOriginal[0]->GetData();
    }

    if (referenceDataForHistogramMatching) {
      for (size_t idx = useTemplateForHistogramMatching ? 0 : 1;
           idx < imageListOriginal.size(); ++idx) {
        imageListOriginal[idx]->GetData()->ApplyFunctionObject(
            cmtk::TypedArrayFunctionHistogramMatching(
                *(imageListOriginal[idx]->GetData()),
                *referenceDataForHistogramMatching));
      }
    }
  }

  functional->InitializeXforms(
      GridSpacing,
      GridSpacingExact);  // must do this before downsampling template grid
  functional->SetRepeatIntensityHistogramMatching(RepeatHistogramMatching);
  const cmtk::Types::Coordinate FinestGridSpacing =
      GridSpacing / (1 << RefineTransformationGrid);

  if (DisableControlPointsMaskPath) {
    functional->SetDisableControlPointsMask(cmtk::UniformVolume::SmartConstPtr(
        cmtk::VolumeIO::Read(DisableControlPointsMaskPath)));
  }

  const double timeBaselineProcess = cmtk::Timers::GetTimeProcess();

  if (!DisableOptimization) {
    cmtk::CoordinateVector v;

    const int downsampleFrom = std::max(DownsampleFrom, DownsampleTo);
    const int downsampleTo = std::min(DownsampleFrom, DownsampleTo);

    for (int downsample = downsampleFrom;
         (downsample >= downsampleTo) || RefineTransformationGrid;
         downsample = downsample ? downsample / 2 : -1) {
      if ((RefineTransformationGrid > 0) && (downsample != downsampleFrom)) {
        functional->RefineTransformationGrids();
        --RefineTransformationGrid;
      }
      functional->GetParamVector(v);

      const int actualDownsample = std::max(downsampleTo, downsample);

      functional->SetTemplateGrid(originalTemplateGrid, actualDownsample,
                                  UseTemplateData);
      cmtk::UniformVolume::SmartPtr templateGrid =
          functional->GetTemplateGrid();

      if (UseSmoothSigmaFactorPixel) {
        functional->SetGaussianSmoothImagesSigma(SmoothSigmaFactorPixel *
                                                 templateGrid->GetMinDelta());
        functional->SetTargetImages(imageListOriginal);
      } else {
        if (UseSmoothSigmaFactorControlPointSpacing) {
          functional->SetGaussianSmoothImagesSigma(
              SmoothSigmaFactorControlPointSpacing * FinestGridSpacing *
              (1 << RefineTransformationGrid));
          functional->SetTargetImages(imageListOriginal);
        }
      }

      cmtk::DebugOutput(1).GetStream().printf(
          "Template grid is %d x %d x %d pixels of size %f x %f x %f\n",
          templateGrid->m_Dims[0], templateGrid->m_Dims[1],
          templateGrid->m_Dims[2], templateGrid->m_Delta[0],
          templateGrid->m_Delta[1], templateGrid->m_Delta[2]);

      if (SamplingDensity > 0) {
        functional->SetProbabilisticSampleDensity(SamplingDensity);
        functional->SetProbabilisticSampleUpdatesAfter(10);
      }

      functional->AllocateStorage();

      cmtk::BestDirectionOptimizer optimizer(OptimizerStepFactor);
      optimizer.SetAggressiveMode(OptimizerAggressive);
      optimizer.SetRepeatLevelCount(OptimizerRepeatLevel);
      optimizer.SetDeltaFThreshold(OptimizerDeltaFThreshold);
      optimizer.SetFunctional(functional);

      cmtk::Types::Coordinate exploration =
          Exploration * templateGrid->GetMinDelta();
      cmtk::Types::Coordinate accuracy = Accuracy * templateGrid->GetMinDelta();
      if ((downsample > downsampleTo) || RefineTransformationGrid)
        accuracy =
            std::max<cmtk::Types::Coordinate>(accuracy, .25 * exploration);

      try {
        // do we have a normal subgroup?
        if (NormalGroupFirstN) {
          // yes: first run normal group by itself
          cmtk::StdErr << "Running normal subgroup...\n";
          functional->SetForceZeroSum(ForceZeroSum);
          functional->SetActiveImagesFromTo(0, NormalGroupFirstN);
          functional->SetActiveXformsFromTo(0, NormalGroupFirstN);
          optimizer.Optimize(v, exploration, accuracy);

          // second: run abnormal group, but keep using normal group's data for
          // reference
          cmtk::StdErr << "Running diseased subgroup...\n";
          functional->SetForceZeroSum(false);  // no point here
          functional->SetActiveImagesFromTo(0, imageListOriginal.size());
          functional->SetActiveXformsFromTo(NormalGroupFirstN,
                                            imageListOriginal.size());
          optimizer.Optimize(v, exploration, accuracy);
        } else {
          optimizer.Optimize(v, exploration, accuracy);
        }
      } catch (cmtk::GroupwiseRegistrationFunctionalBase::BadXform) {
        cmtk::StdErr << "FAILED: at least one image has too few pixels in the "
                        "template area.\n";
        return 1;
      }
    }
  }

  // determine and print CPU time
  const double timeElapsedProcess =
      cmtk::Timers::GetTimeProcess() - timeBaselineProcess;
  cmtk::StdErr.printf("Process CPU time [s]: %f\n", timeElapsedProcess);

  functional->SetTargetImages(imageListOriginal);
  functional->SetTemplateGrid(originalTemplateGrid);

  cmtk::GroupwiseRegistrationOutput output;
  output.SetFunctional(functional);
  output.SetOutputRootDirectory(OutputRootDirectory);

  output.WriteGroupwiseArchive(OutputArchive);
  output.WriteXformsSeparateArchives(OutputStudyListIndividual,
                                     AverageImagePath);
  output.WriteAverageImage(AverageImagePath, AverageImageInterpolation,
                           cmtk::TYPE_FLOAT, UseTemplateData);

  return 0;
}

#include "cmtkSafeMain"
