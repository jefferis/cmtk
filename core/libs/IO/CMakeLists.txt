##
##  Copyright 2004-2010 SRI International
##  Copyright 1997-2009 Torsten Rohlfing
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
SET(cmtkIO_SRCS
  cmtkAffineXformITKIO.cxx
  cmtkSplineWarpXformITKIO.cxx
  cmtkClassStreamAffineXform.cxx
  cmtkClassStreamHistogram.cxx
  cmtkClassStreamInfinitePlane.cxx
  cmtkClassStreamJointHistogram.cxx
  cmtkClassStreamLandmarkList.cxx
  cmtkClassStreamStudyList.cxx
  cmtkClassStreamWarpXform.cxx
  cmtkFileFormat.cxx
  cmtkGroupwiseRegistrationFunctionalIO.cxx
  cmtkImageInfo.cxx
  cmtkMountPoints.cxx
  cmtkPGM.cxx
  cmtkResourceFile.cxx
  cmtkScalarImageIO.cxx
  cmtkSegmentationLabelIO.cxx
  cmtkStudy.cxx
  cmtkStudyImageSet.cxx
  cmtkStudyList.cxx
  cmtkTypedStream.cxx
  cmtkTypedStreamStudylist.cxx
  cmtkVolumeFromFile.cxx
  cmtkVolumeFromFileAnalyze.cxx
  cmtkVolumeFromFileBioRad.cxx
  cmtkVolumeFromFileDICOM.cxx
  cmtkVolumeFromFileMetaFile.cxx
  cmtkVolumeFromFileNRRD.cxx
  cmtkVolumeFromFileNifti.cxx
  cmtkVolumeFromFileVanderbilt.cxx
  cmtkVolumeFromSlices.cxx
  cmtkVolumeFromStudy.cxx
  cmtkVolumeIO.cxx
  cmtkXformIO.cxx
  cmtkXformIONrrd.cxx
)

IF(CMTK_USE_DCMTK)
  SET(cmtkIO_SRCS ${cmtkIO_SRCS} cmtkDICOM.cxx)
ENDIF(CMTK_USE_DCMTK)

IF(CMTK_USE_MPI)
  SET(cmtkIO_SRCS ${cmtkIO_SRCS} cmtkUniformVolumeMPI.cxx)
ENDIF(CMTK_USE_MPI)

ADD_LIBRARY(cmtkIO ${cmtkIO_SRCS})

TARGET_LINK_LIBRARIES(cmtkIO cmtkBase cmtkSystem cmtkNumerics)

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkIO PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkIO RUNTIME DESTINATION lib ARCHIVE DESTINATION lib COMPONENT runtime)
INSTALL(TARGETS cmtkIO LIBRARY DESTINATION lib ARCHIVE DESTINATION lib COMPONENT libraries)
FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION include/cmtk COMPONENT headers)
