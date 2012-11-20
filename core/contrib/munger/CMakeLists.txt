##
##  Copyright 2012 SRI International
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

FIND_PACKAGE(Perl)

IF(PERL_FOUND)

  ## Configure for the build tree
  SET(CMTK_BINARY_DIR_CONFIG ${CMAKE_BINARY_DIR}/bin)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/munger.in ${CMAKE_BINARY_DIR}/bin/munger @ONLY)

  ## Configure for the install tree
  SET(CMTK_BINARY_DIR_CONFIG ${CMAKE_INSTALL_PREFIX}/${CMTK_INSTALL_BIN_DIR})
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/munger.in ${CMAKE_BINARY_DIR}/Install/munger @ONLY)
  INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/Install/munger DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT tools)

ELSE(PERL_FOUND)

  MESSAGE(ERROR "Perl is needed for 'munger' script")

ENDIF(PERL_FOUND)
