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

SET(CMTK_LIBRARIES "cmtkSegmentation;cmtkRegistration;cmtkIO;cmtkBase;cmtkNumerics;cmtkSystem")

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
      convert_warp
      cpmatch
      deformation_field
      film
      filter
      glm
      gregxform
      histogram
      imagemath
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
      regress
      reorient
      runcheck
      sequence
      similarity
      split
      statistics
      sympl
      ttest
      unsplit
      volume_injection
      volume_reconstruction
      warp2ps
      xform2dfield)

    IF(CMTK_ROOT_PATH_SRI24)
      SET(APPS ${APPS} asegment_sri24)
    ENDIF(CMTK_ROOT_PATH_SRI24)

    IF(HAVE_STL_HASH_MAP)
      SET(APPS ${APPS} mcaffine mcwarp)
    ENDIF(HAVE_STL_HASH_MAP)

    IF(CMTK_BUILD_UNSTABLE)
      SET(APPS ${APPS} jidb)
    ENDIF(CMTK_BUILD_UNSTABLE)  

    IF(CMTK_USE_DCMTK)
      SET(APPS ${APPS} dcm2image)
    ENDIF(CMTK_USE_DCMTK)
    
    IF(CMTK_USE_VTK)
      SET(VTK_APPS align_landmarks asegment sympl warp)
    ELSE(CMTK_USE_VTK)
      SET(APPS ${APPS} warp)
    ENDIF(CMTK_USE_VTK)
    
  ENDIF(NOT CMTK_USE_MPI)
ENDIF(CMTK_SINGLE_COMMAND_BINARY)

FOREACH(APPLICATION ${APPS})
  ADD_EXECUTABLE(${APPLICATION} ${CMTK_LIB_INIT_SRCS} ${APPLICATION}.cxx)
  TARGET_LINK_LIBRARIES(${APPLICATION} ${CMTK_LIBRARIES})
  INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION bin)
ENDFOREACH(APPLICATION ${APPS})

FOREACH(APPLICATION ${VTK_APPS})
  ADD_EXECUTABLE(${APPLICATION} ${CMTK_LIB_INIT_SRCS} ${APPLICATION}.cxx)
  INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION bin)
  TARGET_LINK_LIBRARIES(${APPLICATION} cmtkVTKWrapper ${CMTK_LIBRARIES} ${CMTK_VTK_LIBS})
ENDFOREACH(APPLICATION ${VTK_APPS})

