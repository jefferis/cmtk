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
SET(CMTK_LIBRARY_DIRS_CONFIG ${CMTK_BUILD_LIB_DIR})

# Determine the include directories needed.
SET(CMTK_INCLUDE_DIRS_CONFIG
  ${CMTK_INCLUDE_DIRS_BUILD_TREE}
  ${CMTK_INCLUDE_DIRS_SOURCE_TREE}
  ${CMTK_INCLUDE_DIRS_SYSTEM}
)

#-----------------------------------------------------------------------------
# Configure CMTKConfig.cmake for the build tree.
CONFIGURE_FILE(${CMTK_SOURCE_DIR}/CMTKConfig.cmake.in
               ${CMTK_BINARY_DIR}/CMTKConfig.cmake @ONLY IMMEDIATE)

#-----------------------------------------------------------------------------
# Settings specific to the install tree.

# The library dependencies file.
SET(CMTK_LIBRARY_DEPENDS_FILE      CMTKLibraryDepends.cmake)

# Include directories.
SET(CMTK_INCLUDE_DIRS_CONFIG
  ${CMTK_INCLUDE_DIRS_INSTALL_TREE}
  ${CMTK_INCLUDE_DIRS_SYSTEM}
)

# Link directories.
# The install tree will use the directory where CMTKConfig.cmake is found, which
# happens to be "INSTALLATION/lib/InsightToolkit". That is, it is already the
# same directory where the libraries are installed. Therefore this variable
# must be empty here. See CMTKConfig.cmake.in for details on how this variable
# is used.
SET(CMTK_LIBRARY_DIRS_CONFIG "")  

#-----------------------------------------------------------------------------
# Configure CMTKConfig.cmake for the install tree.
CONFIGURE_FILE(${CMTK_SOURCE_DIR}/CMTKConfig.cmake.in
               ${CMTK_BINARY_DIR}/Utilities/CMTKConfig.cmake @ONLY IMMEDIATE)
