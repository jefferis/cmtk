##
##  Copyright 1997-2009 Torsten Rohlfing
##
##  Copyright 2004-2010 SRI International
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
  cmtkAffineRegistrationCommandLine.cxx
  cmtkAffineRegistration.cxx
  cmtkBestDirectionOptimizer.cxx
  cmtkBestNeighbourOptimizer.cxx
  cmtkCongealingFunctional.cxx
  cmtkElasticRegistrationCommandLine.cxx
  cmtkElasticRegistration.cxx
  cmtkFunctionalAffine2D.cxx
  cmtkGroupwiseRegistrationFunctionalAffineInitializer.cxx
  cmtkGroupwiseRegistrationFunctionalBase.cxx
  cmtkGroupwiseRegistrationFunctionalXformTemplate.cxx
  cmtkGroupwiseRegistrationFunctionalXformTemplate_Affine_IO.cxx
  cmtkGroupwiseRegistrationFunctionalXformTemplate_SplineWarp.cxx
  cmtkGroupwiseRegistrationFunctionalXformTemplate_SplineWarp_IO.cxx
  cmtkGroupwiseRegistrationFunctionalXformTemplateBase.cxx
  cmtkGroupwiseRegistrationOutput.cxx
  cmtkGroupwiseRegistrationRMIFunctional.cxx
  cmtkImagePairRegistration.cxx
  cmtkImagePairRegistrationFunctional.cxx
  cmtkImagePairRegistrationImagePreprocessor.cxx
  cmtkImagePairAffineRegistration.cxx
  cmtkImagePairAffineRegistrationCommandLine.cxx
  cmtkImagePairAffineRegistrationFunctional.cxx
  cmtkImagePairNonrigidRegistration.cxx
  cmtkImagePairNonrigidRegistrationCommandLine.cxx
  cmtkImagePairNonrigidRegistrationFunctional.cxx
  cmtkImagePairSimilarityJointHistogram.cxx
  cmtkImagePairSimilarityMeasure.cxx
  cmtkImagePairSimilarityMeasureCR.cxx
  cmtkImagePairSimilarityMeasureMSD.cxx
  cmtkImagePairSimilarityMeasureNCC.cxx
  cmtkImagePairSimilarityMeasureRMS.cxx
  cmtkImagePairSymmetricNonrigidRegistrationFunctional.cxx
  cmtkImageSymmetryPlaneCommandLineBase.cxx
  cmtkImageSymmetryPlaneFunctional.cxx
  cmtkImageSymmetryPlaneFunctionalBase.cxx
  cmtkIterativeDirectionOptimizer.cxx
  cmtkMakeInitialAffineTransformation.cxx
  cmtkMultiChannelRegistrationFunctionalBase.cxx
  cmtkMultiLevelOptimizer.cxx
  cmtkProtocolCallback.cxx
  cmtkReformatVolume.cxx
  cmtkRegistration2d2d.cxx
  cmtkRegistrationCallback.cxx
  cmtkRegistrationJointHistogram.cxx
  cmtkScalarImageSimilarity.cxx
  cmtkScalarImageSimilarityRMI.cxx
  cmtkSplineWarpCongealingFunctional.cxx
  cmtkSplineWarpCongealingFunctionalStaticThreadStorage.cxx
  cmtkSplineWarpGroupwiseRegistrationRMIFunctional.cxx
  cmtkSymmetricElasticFunctional.cxx
  cmtkSymmetryPlaneFunctional.cxx
  cmtkTypedArraySimilarity.cxx
  cmtkTypedArraySimilarityMemory.cxx
  cmtkTypedArraySimilarityRMI.cxx
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
)

IF(CMTK_USE_SQLITE)
  SET(cmtkRegistration_SRCS ${cmtkRegistration_SRCS}
    cmtkImageXformDB.cxx
    cmtkRegistrationDatabaseUpdater.cxx)
ENDIF(CMTK_USE_SQLITE)

ADD_LIBRARY(cmtkRegistration ${cmtkRegistration_SRCS})

TARGET_LINK_LIBRARIES(cmtkRegistration cmtkIO cmtkBase cmtkSystem cmtkNumerics )

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkRegistration PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkRegistration
  RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT runtime
  LIBRARY DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT runtime
  ARCHIVE DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION ${CMTK_INSTALL_INCLUDE_DIR}/Registration COMPONENT headers)
