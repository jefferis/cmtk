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


PROJECT(cmtkDCMTK)
INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(dcmtk_SRCS dcbytstr.cc dcerror.cc  dcobject.cc dctagkey.cc didocu.cc   ofconsol.cc
               dcdatset.cc dcfilefo.cc dcpcache.cc dctypes.cc  diutils.cc  oflist.cc
               dcdicent.cc dchashdi.cc dcsequen.cc dcvr.cc     djdecode.cc ofstd.cc
               dcdict.cc   dcitem.cc   dcstack.cc  dcvrui.cc   djutils.cc  ofstring.cc
               dcelem.cc   dclist.cc   dctag.cc    dcxfer.cc   ofcond.cc   ofthread.cc )

CONFIGURE_FILE(${DCMTK_INCLUDE_DIR}/config/cfunix.h.in ${CMAKE_CURRENT_BINARY_DIR}/cfunix.h)

ADD_LIBRARY(cmtkDCMTK ${dcmtk_SRCS})
TARGET_LINK_LIBRARIES(cmtkDCMTK)

INSTALL_TARGETS(/lib cmtkDCMTK)
INSTALL_FILES(/include/dcmtk "(\\.h)$")
INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_DIR}/config)
INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_DIR}/dcmdata)
INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_DIR}/dcmimgle)
INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_DIR}/dcmjpeg)
INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_DIR}/ofstd)

SET(CMTK_DCMTK_LIBRARIES cmtkDCMTK)



