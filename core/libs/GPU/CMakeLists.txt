##
##  Copyright 1997-2011 Torsten Rohlfing
##
##  Copyright 2004-2011 SRI International
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

# Sources of non-specific classes.
SET(cmtkGPU_SRCS
  cmtkDeviceHistogram.cxx
  cmtkDeviceUniformVolume.cxx
  cmtkDeviceUniformVolumeArray.cxx
  cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice.cxx
  cmtkImagePairAffineRegistrationFunctionalDevice.cxx
  cmtkImageSymmetryPlaneFunctionalDevice.cxx
  cmtkSimpleLevelsetDevice.cxx
  )

IF(CMTK_USE_CUDA)
  SET(cmtkGPU_SRCS_CUDA
    cmtkDeviceImageConvolution_kernels.cu
    cmtkDeviceHistogram_kernels.cu
    cmtkDeviceMemoryCUDA.cxx
    cmtkDeviceArrayCUDA.cxx
    cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_kernels.cu
    cmtkImagePairAffineRegistrationFunctionalDevice_kernels.cu
    cmtkSimpleLevelsetDevice_kernels.cu
    cmtkSumReduction_kernel.cu
    )

  CUDA_ADD_LIBRARY(cmtkGPU ${cmtkGPU_SRCS_CUDA} ${cmtkGPU_SRCS})
ENDIF(CMTK_USE_CUDA)

TARGET_LINK_LIBRARIES(cmtkGPU cmtkSegmentation cmtkRegistration cmtkBase cmtkSystem cmtkNumerics)

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkGPU PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkGPU
  RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT runtime
  LIBRARY DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT runtime
  ARCHIVE DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION ${CMTK_INSTALL_INCLUDE_DIR}/GPU COMPONENT headers)
