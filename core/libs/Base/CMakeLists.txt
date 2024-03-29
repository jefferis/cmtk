##
##  Copyright 1997-2010 Torsten Rohlfing
##
##  Copyright 2004-2014 SRI International
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
SET(cmtkBase_SRCS
  cmtkActiveDeformationModel.cxx
  cmtkActiveShapeModel.cxx
  cmtkAffineXform.cxx
  cmtkAffineXformUniformVolume.cxx
  cmtkAnatomicalOrientation.cxx
  cmtkAnatomicalOrientationBase.cxx
  cmtkAnatomicalOrientationPermutationMatrix.cxx
  cmtkBitVector.cxx
  cmtkCompatibilityMatrix4x4.cxx
  cmtkDataGrid.cxx
  cmtkDataGrid_Crop.cxx
  cmtkDataGridFilter.cxx
  cmtkDataGridConnectedComponents.cxx
  cmtkDataGridLocalCorrelation.cxx
  cmtkDataGridMorphologicalOperators.cxx
  cmtkDeformationField.cxx
  cmtkDeformationField_Jacobian.cxx
  cmtkDirectionSet.cxx
  cmtkFilterVolume.cxx
  cmtkFitAffineToLandmarks.cxx
  cmtkFitAffineToWarpXform.cxx
  cmtkFitAffineToXformList.cxx
  cmtkFitPolynomialToLandmarks.cxx
  cmtkFitRigidToLandmarks.cxx
  cmtkFitSplineWarpToDeformationField.cxx
  cmtkFitSplineWarpToLandmarks.cxx
  cmtkFitSplineWarpToXformList.cxx
  cmtkFitToXformListBase.cxx
  cmtkFunctional.cxx
  cmtkGeneralLinearModel.cxx
  cmtkHistogram.cxx
  cmtkHypothesisTests.cxx
  cmtkImageOperation.cxx
  cmtkImageOperationCropRegion.cxx
  cmtkImageOperationCropThreshold.cxx
  cmtkImageOperationDistanceMap.cxx
  cmtkImageOperationDownsample.cxx
  cmtkImageOperationHistogramEqualization.cxx
  cmtkImageOperationMapValues.cxx
  cmtkImageOperationRegionFilter.cxx
  cmtkImageOperationResampleIsotropic.cxx
  cmtkImageOperationScaleToRange.cxx
  cmtkImageOperationThreshold.cxx
  cmtkImageOperationOtsuThreshold.cxx
  cmtkImageOperationPruneHistogram.cxx
  cmtkImageOperationRevert.cxx
  cmtkJointHistogram.cxx
  cmtkJointHistogramBase.cxx
  cmtkLandmark.cxx
  cmtkLandmarkPair.cxx
  cmtkLandmarkList.cxx
  cmtkLandmarkPairList.cxx
  cmtkMagphanEMR051.cxx
  cmtkMagphanEMR051_GetPhantomImage.cxx
  cmtkMathFunctionWrappers.cxx
  cmtkMathUtil_LinAlg.cxx
  cmtkMathUtil_Random.cxx
  cmtkMathUtil_Statistics.cxx
  cmtkMatrix3x3.cxx
  cmtkMatrix4x4.cxx
  cmtkMetaInformationObject.cxx
  cmtkParametricPlane.cxx
  cmtkPolynomialXform.cxx
  cmtkScalarImage.cxx
  cmtkScalarImageGradientField.cxx
  cmtkSegmentationLabel.cxx
  cmtkSplineWarpXform.cxx
  cmtkSplineWarpXform_Inverse.cxx
  cmtkSplineWarpXform_Jacobian.cxx
  cmtkSplineWarpXform_Rigidity.cxx
  cmtkSplineWarpXformUniformVolume.cxx
  cmtkTransformChangeToSpaceAffine.cxx
  cmtkTransformChangeFromSpaceAffine.cxx
  cmtkTransformedVolumeAxes.cxx
  cmtkTypedArray.cxx
  cmtkTypedArray_Statistics.cxx
  cmtkTypedArrayFunctionHistogramMatching.cxx
  cmtkTypedArrayFunctionHistogramEqualization.cxx
  cmtkTypedArrayNoiseEstimatorNaiveGaussian.cxx
  cmtkTypes.cxx
  cmtkUniformDistanceMap.cxx
  cmtkUniformVolume.cxx
  cmtkUniformVolume_Crop.cxx
  cmtkUniformVolume_Differential.cxx
  cmtkUniformVolume_Resample.cxx
  cmtkUniformVolume_Space.cxx
  cmtkUniformVolumeGaussianFilter.cxx
  cmtkUniformVolumeLaplaceFilter.cxx
  cmtkUniformVolumeInterpolatorBase.cxx
  cmtkUniformVolumeInterpolatorPartialVolume.cxx
  cmtkUniformVolumeMorphologicalOperators.cxx
  cmtkUniformVolumePainter.cxx
  cmtkVolumeClipping.cxx
  cmtkVolume.cxx
  cmtkVolumeGridToGridLookup.cxx
  cmtkWarpXform.cxx
  cmtkXform.cxx
  cmtkXform_Inverse.cxx
  cmtkXformList.cxx
  cmtkXformListEntry.cxx
)

ADD_LIBRARY(cmtkBase ${cmtkBase_SRCS})

TARGET_LINK_LIBRARIES(cmtkBase cmtkSystem cmtkNumerics )

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkBase PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkBase
  RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT runtime
  LIBRARY DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT runtime
  ARCHIVE DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx}  DESTINATION ${CMTK_INSTALL_INCLUDE_DIR}/Base/ COMPONENT headers)
