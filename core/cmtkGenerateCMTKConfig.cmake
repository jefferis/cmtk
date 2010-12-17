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

##
## This file borrows heavily from the analogous InsightToolkit file
##

# Generate the CMTKConfig.cmake file in the build tree.  Also configure
# one for installation.  The file tells external projects how to use
# CMTK.

#-----------------------------------------------------------------------------
# Settings common to the build and installation tree.

# The "use" file.
SET(CMTK_USE_FILE                  UseCMTK.cmake)

# The build settings file.
SET(CMTK_BUILD_SETTINGS_FILE       CMTKBuildSettings.cmake)


#-----------------------------------------------------------------------------
# Settings specific to the build tree.

# The library dependencies file.
SET(CMTK_LIBRARY_DEPENDS_FILE       CMTKLibraryDepends.cmake)

# Library directory.
SET(CMTK_LIBRARY_DIRS_CONFIG ${CMTK_LIBRARY_PATH})

# Binary directory.
SET(CMTK_BINARY_DIR_CONFIG ${CMTK_BINARY_DIR})

# Determine the include directories needed.
SET(CMTK_INCLUDE_DIRS_CONFIG 
  ${CMTK_INCLUDE_DIRS_BUILD_TREE} 
  ${CMTK_INCLUDE_DIRS_SYSTEM})

# Set data directory
SET(CMTK_DATA_ROOT_CONFIG ${CMTK_DATA_ROOT})

# Set DICOM dictionary paths for build AND install tree
SET(CMTK_DCMDICTPATH_CONFIG ${CMTK_LIBRARY_PATH})
SET(CMTK_DCMDICTPATH_INSTALL_CONFIG ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_LIB_DIR})

#-----------------------------------------------------------------------------
# Configure CMTKConfig.cmake for the build tree.
CONFIGURE_FILE(${CMTK_SOURCE_DIR}/CMTKConfig.cmake.in
  ${CMTK_BINARY_DIR}/CMTKConfig.cmake @ONLY IMMEDIATE)
CONFIGURE_FILE(${CMTK_SOURCE_DIR}/cmtkconfig.h.cmake
  ${CMTK_BINARY_DIR}/cmtkconfig.h @ONLY IMMEDIATE)

#-----------------------------------------------------------------------------
# Settings specific to the install tree.

# The library dependencies file.
SET(CMTK_LIBRARY_DEPENDS_FILE      CMTKLibraryDepends.cmake)

# Include directories.
SET(CMTK_INCLUDE_DIRS_CONFIG 
  \${CMTK_INSTALL_PREFIX}${CMTK_INSTALL_INCLUDE_DIR}
)

# List of CMTK libraries
SET(CMTK_LIBRARIES cmtkIO cmtkPipeline cmtkQt cmtkRegistration cmtkSegmentation cmtkRecon cmtkBase cmtkSystem cmtkNumerics)

IF(CMTK_BUILD_UNSTABLE)
  SET(CMTK_LIBRARIES cmtkUnstable ${CMTK_LIBRARIES})
ENDIF(CMTK_BUILD_UNSTABLE)

IF(CMTK_USE_CUDA)
  SET(CMTK_LIBRARIES cmtkGPU ${CMTK_LIBRARIES})
ENDIF(CMTK_USE_CUDA)

# Link directories.
SET(CMTK_LIBRARY_DIRS_CONFIG "\${CMTK_INSTALL_PREFIX}${CMTK_INSTALL_LIB_DIR}")

# Link directories.
# The install tree will use the directory where CMTKConfig.cmake is found, which
# happens to be "INSTALLATION/lib". That is, it is already the
# same directory where the libraries are installed. Therefore this variable
# must be empty here. See CMTKConfig.cmake.in for details on how this variable
# is used.
SET(CMTK_LIBRARY_DIRS_CONCUR "")  

# Binary directory.
SET(CMTK_BINARY_DIR_CONFIG ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_BIN_DIR})

#-----------------------------------------------------------------------------
# Configure CMTKConfig.cmake for the install tree.

# Construct the proper number of GET_FILENAME_COMPONENT(... PATH)
# calls to compute the installation prefix.
STRING(REGEX REPLACE "/" ";" CMTK_INSTALL_PACKAGE_DIR_COUNT "${CMTK_INSTALL_PACKAGE_DIR}")
SET(CMTK_CONFIG_CODE "
# Compute the installation prefix from this CMTKConfig.cmake file location.
GET_FILENAME_COMPONENT(CMTK_INSTALL_PREFIX \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)")
FOREACH(p ${CMTK_INSTALL_PACKAGE_DIR_COUNT})
  SET(CMTK_CONFIG_CODE
    "${CMTK_CONFIG_CODE}\nGET_FILENAME_COMPONENT(CMTK_INSTALL_PREFIX \"\${CMTK_INSTALL_PREFIX}\" PATH)"
    )
ENDFOREACH(p)


CONFIGURE_FILE(${CMTK_SOURCE_DIR}/CMTKConfig.cmake.in 
  ${CMTK_BINARY_DIR}/Install/CMTKConfig.cmake @ONLY IMMEDIATE)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/Install/CMTKConfig.cmake DESTINATION lib COMPONENT headers)

CONFIGURE_FILE(${CMTK_SOURCE_DIR}/cmtkconfig.h.cmake 
  ${CMTK_BINARY_DIR}/Install/cmtkconfig.h @ONLY IMMEDIATE)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/Install/cmtkconfig.h DESTINATION include COMPONENT headers)

