##
##  Copyright 1997-2012 Torsten Rohlfing
##
##  Copyright 2004-2012 SRI International
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

SET(CMTK_LIBRARIES "cmtkSegmentation;cmtkRecon;cmtkRegistration;cmtkIO;cmtkBase;cmtkNumerics;cmtkSystem")
IF(CMTK_BUILD_UNSTABLE)
  SET(CMTK_LIBRARIES "cmtkUnstable;${CMTK_LIBRARIES}")
ENDIF(CMTK_BUILD_UNSTABLE)

SET(APPS avg_adm
  describe
  groupwise_init
  groupwise_affine
  groupwise_warp
  asegment
  average_affine
  average_images
  average_labels
  concat_affine
  convertx
  convert_warp
  destripe
  dof2mat
  epiunwarp
  fibxform
  film
  filter
  fit_affine_xform
  fit_spline_xform
  glm
  gmm
  gregxform
  hausdorff
  histogram
  imagemath
  jidb
  levelset
  lsba
  lmsba
  lvote
  make_initial_affine
  mat2dof
  mip
  mk_adni_phantom
  mk_analyze_hdr
  mk_nifti_hdr
  mk_phantom_3d
  mrbias
  overlap
  probe
  reformatx
  registration
  registrationx
  regress
  reorient
  sba
  sbai
  sequence
  similarity
  split
  statistics
  sympl
  symplx
  ttest
  unsplit
  volume_injection
  volume_reconstruction
  vtkxform
  warp
  warpx
  warp2ps
  xform2dfield
  xform2itk
  xform2scalar)

IF(CMTK_USE_CUDA)
  SET(CUDA_APPS 
    levelset_cuda
    mrbias_cuda
    symplx_cuda)
ENDIF(CMTK_USE_CUDA)

IF(CMTK_USE_SQLITE)
  LIST(APPEND APPS dbtool)
ENDIF(CMTK_USE_SQLITE)

IF(CMTK_ROOT_PATH_SRI24)
  LIST(APPEND APPS asegment_sri24)
ENDIF(CMTK_ROOT_PATH_SRI24)

IF(HAVE_STL_HASH_MAP)
  LIST(APPEND APPS mcaffine mcwarp)
ENDIF(HAVE_STL_HASH_MAP)

IF(CMTK_USE_DCMTK)
  LIST(APPEND APPS dcm2image)
ENDIF(CMTK_USE_DCMTK)

FOREACH(APPLICATION ${APPS})
  ADD_EXECUTABLE(${APPLICATION} ${APPLICATION}.cxx)
  TARGET_LINK_LIBRARIES(${APPLICATION} ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})

  IF(CMTK_BUILD_WRAPPER)
    INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION ${CMTK_INSTALL_LIB_DIR}/bin COMPONENT tools)
    SET(manpage cmtk-${APPLICATION}.1)
  ELSE(CMTK_BUILD_WRAPPER)    
    INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT tools)
    SET(manpage ${APPLICATION}.1)
  ENDIF(CMTK_BUILD_WRAPPER)
  
  IF(BUILD_MANPAGES)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${manpage}
      COMMAND ${APPLICATION} --man ">" ${manpage}
      DEPENDS ${APPLICATION})
    ADD_CUSTOM_TARGET(man_${APPLICATION} ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${manpage})
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${manpage} DESTINATION share/man/man1 COMPONENT tools)  
  ENDIF(BUILD_MANPAGES)
ENDFOREACH(APPLICATION ${APPS})

FOREACH(APPLICATION ${CUDA_APPS})
  CUDA_ADD_EXECUTABLE(${APPLICATION} ${APPLICATION}.cxx)
  TARGET_LINK_LIBRARIES(${APPLICATION} cmtkGPU ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})

  IF(CMTK_BUILD_WRAPPER)
    INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION ${CMTK_INSTALL_LIB_DIR}/bin COMPONENT tools)
  ELSE(CMTK_BUILD_WRAPPER)    
    INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT tools)
  ENDIF(CMTK_BUILD_WRAPPER)
ENDFOREACH(APPLICATION ${CUDA_APPS})

IF(CMTK_USE_CUDA)
  CUDA_ADD_EXECUTABLE(cudatest cudatest.cxx)

  IF(CMTK_BUILD_WRAPPER)
    INSTALL(TARGETS cudatest RUNTIME DESTINATION ${CMTK_INSTALL_LIB_DIR}/bin COMPONENT tools)
  ELSE(CMTK_BUILD_WRAPPER)    
    INSTALL(TARGETS cudatest RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT tools)
  ENDIF(CMTK_BUILD_WRAPPER)
ENDIF(CMTK_USE_CUDA)
