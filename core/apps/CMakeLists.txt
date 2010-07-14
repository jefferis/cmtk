##
##  Copyright 1997-2010 Torsten Rohlfing
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

SET(CMTK_LIBRARIES "cmtkSegmentation;cmtkRecon;cmtkRegistration;cmtkIO;cmtkBase;cmtkNumerics;cmtkSystem")

INCLUDE(${CMAKE_SOURCE_DIR}/libs/cmtkLibInitSrcs.cmake)

IF(CMTK_SINGLE_COMMAND_BINARY)
  SET(APPS cmtk)
ELSE(CMTK_SINGLE_COMMAND_BINARY)
  SET(APPS avg_adm
    congeal
    congeal_warp
    describe
    groupwise_rmi
    groupwise_init
    groupwise_rmi_warp
    )
  
  IF(CMTK_USE_CUDA)
    SET(CUDA_APPS 
      mrbias_cuda)
  ENDIF(CMTK_USE_CUDA)
  
  IF(CMTK_USE_OPENCL)
  ENDIF(CMTK_USE_OPENCL)
  
  IF(NOT CMTK_USE_MPI)
    SET(APPS ${APPS}
      asegment
      average_affine
      average_edt
      average_grey
      average_images
      average_labels
      concat_affine
      convert
      convertx
      convert_warp
      film
      filter
      glm
      gregxform
      histogram
      imagemath
      jidb
      levelset
      make_initial_affine
      mat2dof
      mip
      mk_analyze_hdr
      mk_phantom_3d
      mrbias
      overlap
      probe_xform
      randompxls
      randomwarp
      reformatx
      registration
      registrationx
      regress
      reorient
      runcheck
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
      xform2scalar)

    IF(CMTK_BUILD_SQLITE)
      SET(APPS ${APPS}
	dbtool)
    ENDIF(CMTK_BUILD_SQLITE)

    IF(CMTK_ROOT_PATH_SRI24)
      SET(APPS ${APPS} asegment_sri24)
    ENDIF(CMTK_ROOT_PATH_SRI24)

    IF(HAVE_STL_HASH_MAP)
      SET(APPS ${APPS} mcaffine mcwarp)
    ENDIF(HAVE_STL_HASH_MAP)

    IF(CMTK_USE_DCMTK)
      SET(APPS ${APPS} dcm2image)
    ENDIF(CMTK_USE_DCMTK)
    
  ENDIF(NOT CMTK_USE_MPI)
ENDIF(CMTK_SINGLE_COMMAND_BINARY)

FOREACH(APPLICATION ${APPS})
  ADD_EXECUTABLE(${APPLICATION} ${CMTK_LIB_INIT_SRCS} ${APPLICATION}.cxx)
  TARGET_LINK_LIBRARIES(${APPLICATION} ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})
  INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION bin COMPONENT applications)
ENDFOREACH(APPLICATION ${APPS})

FOREACH(APPLICATION ${CUDA_APPS})
  CUDA_ADD_EXECUTABLE(${APPLICATION} ${APPLICATION}.cxx)
  TARGET_LINK_LIBRARIES(${APPLICATION} cmtkGPU ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})
  INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION bin COMPONENT applications)
ENDFOREACH(APPLICATION ${CUDA_APPS})

IF(CMTK_USE_CUDA)
  CUDA_ADD_EXECUTABLE(cudatest cudatest.cxx)
  INSTALL(TARGETS cudatest RUNTIME DESTINATION bin COMPONENT applications)
ENDIF(CMTK_USE_CUDA)

IF(CMTK_USE_OPENCL)
  ADD_EXECUTABLE(cltest cltest.cxx)
  TARGET_LINK_LIBRARIES(cltest ${OPENCL_LIBRARIES})
  INSTALL(TARGETS cltest RUNTIME DESTINATION bin COMPONENT applications)
ENDIF(CMTK_USE_OPENCL)
