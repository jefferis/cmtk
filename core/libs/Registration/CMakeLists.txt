##
##  Copyright 1997-2009 Torsten Rohlfing
##  Copyright 2004-2009 SRI International
##
##  This file is part of the Computational Morphometry Toolkit.
##
##  http://www.nitrc.org/projects/cmtk/
##
##  The Computational Morphometry Toolkit is free software: you can
##  redistribute it and/or modify it under the terms of the GNU General Public
##  License as published by the Free Software Foundation, either version 3 of
##  the License, or (at your option) any later version.
##
##  The Computational Morphometry Toolkit is distributed in the hope that it
##  will be useful, but WITHOUT ANY WARRANTY; without even the implied
##  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with the Computational Morphometry Toolkit.  If not, see
##  <http://www.gnu.org/licenses/>.
##
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

# Sources of non-templated classes.
SET(cmtkRegistration_SRCS
  cmtkAffineCongealingFunctional.cxx
  cmtkAffineCongealingFunctionalIO.cxx
  cmtkAffineGroupwiseRegistrationRMIFunctional.cxx
  cmtkAffineGroupwiseRegistrationRMIFunctionalIO.cxx
  cmtkAffineRegistrationCommandLine.cxx
  cmtkAffineRegistration.cxx
  cmtkBestDirectionOptimizer.cxx
  cmtkBestNeighbourOptimizer.cxx
  cmtkCongealingFunctionalBase.cxx
  cmtkCongealingFunctional.cxx
  cmtkDirectionSetOptimizer.cxx
  cmtkElasticRegistrationCommandLine.cxx
  cmtkElasticRegistration.cxx
  cmtkFunctionalAffine2D.cxx
  cmtkGroupwiseRegistrationFunctionalAffineInitializer.cxx
  cmtkGroupwiseRegistrationFunctionalBase.cxx
  cmtkGroupwiseRegistrationOutput.cxx
  cmtkGroupwiseRegistrationRMIFunctional.cxx
  cmtkImagePairAffineRegistrationFunctional.cxx
  cmtkImagePairSimilarityMeasure.cxx
  cmtkImagePairSimilarityMeasureMSD.cxx
  cmtkInverseInterpolationVolumeReconstructionBase.cxx
  cmtkIterativeDirectionOptimizer.cxx
  cmtkMakeInitialAffineTransformation.cxx
  cmtkMultiLevelOptimizer.cxx
  cmtkProtocolCallback.cxx
  cmtkReformatVolume.cxx
  cmtkReformatVolumeJacobian.cxx
  cmtkReformatVolumePlain.cxx
  cmtkRegistration2d2d.cxx
  cmtkRegistrationCallback.cxx
  cmtkRegistrationJointHistogram.cxx
  cmtkScalarImageSimilarity.cxx
  cmtkScalarImageSimilarityRMI.cxx
  cmtkSplineWarpCongealingFunctional.cxx
  cmtkSplineWarpCongealingFunctionalIO.cxx
  cmtkSplineWarpCongealingFunctionalStaticThreadStorage.cxx
  cmtkSplineWarpGroupwiseRegistrationRMIFunctional.cxx
  cmtkSplineWarpGroupwiseRegistrationRMIFunctionalIO.cxx
  cmtkSymmetricElasticFunctional.cxx
  cmtkSymmetryPlaneFunctional.cxx
  cmtkTypedArraySimilarity.cxx
  cmtkTypedArraySimilarityMemory.cxx
  cmtkTypedArraySimilarityRMI.cxx
  cmtkVolumeInjectionReconstruction.cxx
  cmtkVoxelMatchingAffineFunctional.cxx
  cmtkVoxelMatchingCorrRatio.cxx
  cmtkVoxelMatchingCrossCorrelation.cxx
  cmtkVoxelMatchingElasticFunctional.cxx
  cmtkVoxelMatchingFunctional.cxx
  cmtkVoxelMatchingMeanSquaredDifference.cxx
  cmtkVoxelMatchingMetric.cxx
  cmtkVoxelMatchingMetric_Type.cxx
  cmtkVoxelRegistration.cxx
  cmtkVoxelRegistrationImagePreprocessor.cxx
  cmtkXformList.cxx
  cmtkXformListEntry.cxx
)

ADD_LIBRARY(cmtkRegistration ${cmtkRegistration_SRCS})

TARGET_LINK_LIBRARIES(cmtkRegistration cmtkIO cmtkBase cmtkSystem cmtkNumerics )

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkRegistration PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkRegistration RUNTIME DESTINATION lib ARCHIVE DESTINATION lib COMPONENT runtime)
INSTALL(TARGETS cmtkRegistration LIBRARY DESTINATION lib ARCHIVE DESTINATION lib COMPONENT libraries)
FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION include/cmtk COMPONENT headers)
