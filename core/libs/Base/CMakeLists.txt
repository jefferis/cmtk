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
SET(cmtkBase_SRCS
  cmtkActiveDeformationModel.cxx
  cmtkActiveShapeModel.cxx
  cmtkAffineXform.cxx
  cmtkAffineXformUniformVolume.cxx
  cmtkAnatomicalOrientation.cxx
  cmtkAnatomicalOrientationPermutationMatrix.cxx
  cmtkBitVector.cxx
  cmtkControlPoint.cxx
  cmtkDataGrid.cxx
  cmtkDataGridFilter.cxx
  cmtkDataGridMorphologicalOperators.cxx
  cmtkDataTypeTraits.cxx
  cmtkDeformationField.cxx
  cmtkDeformationFieldJacobian.cxx
  cmtkDirectionSet.cxx
  cmtkFilterVolume.cxx
  cmtkFilterVolumeCoupe.cxx
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
  cmtkImageOperationMedianFilter.cxx
  cmtkImageOperationMedialSkeleton.cxx
  cmtkImageOperationScaleToRange.cxx
  cmtkImageOperationThreshold.cxx
  cmtkImageOperationRevert.cxx
  cmtkInfinitePlane.cxx
  cmtkIntersection.cxx
  cmtkJointHistogram.cxx
  cmtkJointHistogramBase.cxx
  cmtkLandmark.cxx
  cmtkLandmarkList.cxx
  cmtkMatchedLandmarkList.cxx
  cmtkMathUtilLinAlg.cxx
  cmtkMathUtilRandom.cxx
  cmtkMathUtilStatistics.cxx
  cmtkMatrix3x3.cxx
  cmtkMatrix4x4.cxx
  cmtkMetaInformationObject.cxx
  cmtkScalarImage.cxx
  cmtkSegmentationLabel.cxx
  cmtkSplineWarpXform.cxx
  cmtkSplineWarpXformInverse.cxx
  cmtkSplineWarpXformJacobian.cxx
  cmtkSplineWarpXformRigidity.cxx
  cmtkSplineWarpXformUniformVolume.cxx
  cmtkTransformChangeToSpaceAffine.cxx
  cmtkTransformChangeFromSpaceAffine.cxx
  cmtkTypedArray.cxx
  cmtkTypedArrayFunctionHistogramMatching.cxx
  cmtkTypedArrayFunctionHistogramEqualization.cxx
  cmtkTypedArrayNoiseEstimatorNaiveGaussian.cxx
  cmtkTypedArrayStatistics.cxx
  cmtkTypes.cxx
  cmtkUniformDistanceMap.cxx
  cmtkUniformVolume.cxx
  cmtkUniformVolumeCrop.cxx
  cmtkUniformVolumeDifferential.cxx
  cmtkUniformVolumeFilter.cxx
  cmtkUniformVolumeInterpolatorPartialVolume.cxx
  cmtkUniformVolumePainter.cxx
  cmtkUniformVolumeResample.cxx
  cmtkVector3D.cxx
  cmtkVolumeClipping.cxx
  cmtkVolume.cxx
  cmtkVolumeAxesHash.cxx
  cmtkVolumeGridToGridLookup.cxx
  cmtkWarpXform.cxx
  cmtkXform.cxx
)

ADD_LIBRARY(cmtkBase ${cmtkBase_SRCS})

TARGET_LINK_LIBRARIES(cmtkBase cmtkSystem cmtkNumerics )

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkBase PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkBase
  RUNTIME DESTINATION bin COMPONENT runtime
  LIBRARY DESTINATION lib COMPONENT runtime
  ARCHIVE DESTINATION lib COMPONENT libraries)

FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION include/cmtk COMPONENT headers)
