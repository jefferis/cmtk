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

CHECK_INCLUDE_FILE_CXX (sstream HAVE_SSTREAM)
CHECK_INCLUDE_FILE_CXX (strstream HAVE_STRSTREAM)
CHECK_INCLUDE_FILE_CXX (sstream.h HAVE_SSTREAM_H)
CHECK_INCLUDE_FILE_CXX (strstrea.h HAVE_STRSTREA_H)
CHECK_INCLUDE_FILE_CXX (strstream.h HAVE_STRSTREAM_H)

SET(USE_STD_CXX_INCLUDES 1)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/dcmtk/config/cfunix.h.in ${CMAKE_CURRENT_BINARY_DIR}/cfunix.h)

INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/dcmtk/config
	${CMAKE_CURRENT_SOURCE_DIR}/dcmtk/ofstd
	${CMAKE_CURRENT_SOURCE_DIR}/dcmtk/dcmdata
	${CMAKE_CURRENT_SOURCE_DIR}/dcmtk/dcmimgle
)

INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(dcmtk_SRCS
  cmdlnarg.cc
  dcbytstr.cc
  dcchrstr.cc
  dccodec.cc
  dcdatset.cc
  dcdebug.cc
  dcdicdir.cc
  dcdicent.cc
  dcdict.cc
  dcdictzz.cc
  dcdirrec.cc
  dcelem.cc
  dcerror.cc
  dcfilefo.cc
  dchashdi.cc
  dcistrma.cc
  dcistrmf.cc
  dcitem.cc
  dclist.cc
  dcmetinf.cc
  dcobject.cc
  dcostrma.cc
  dcostrmf.cc
  dcpcache.cc
  dcpixel.cc
  dcpixseq.cc
  dcpxitem.cc
  dcsequen.cc
  dcstack.cc
  dcswap.cc
  dctag.cc
  dctagkey.cc
  dctypes.cc
  dcuid.cc
  dcvm.cc
  dcvr.cc
  dcvrae.cc
  dcvras.cc
  dcvrat.cc
  dcvrcs.cc
  dcvrda.cc
  dcvrds.cc
  dcvrdt.cc
  dcvrfd.cc
  dcvrfl.cc
  dcvris.cc
  dcvrlo.cc
  dcvrlt.cc
  dcvrobow.cc
  dcvrof.cc
  dcvrpn.cc
  dcvrpobw.cc
  dcvrsh.cc
  dcvrsl.cc
  dcvrss.cc
  dcvrst.cc
  dcvrtm.cc
  dcvrui.cc
  dcvrul.cc
  dcvrulup.cc
  dcvrus.cc
  dcvrut.cc
  dcxfer.cc
  didocu.cc
  diutils.cc
  ofcond.cc
  ofconsol.cc
  ofdate.cc
  ofdatime.cc
  oflist.cc
  ofstd.cc
  ofstring.cc
  ofthread.cc
  oftime.cc
)


ADD_LIBRARY(cmtkDCMTK ${dcmtk_SRCS})
TARGET_LINK_LIBRARIES(cmtkDCMTK)

INSTALL_TARGETS(/lib cmtkDCMTK)
INSTALL_FILES(/include/dcmtk "(\\.h)$")
