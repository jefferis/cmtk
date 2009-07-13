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

INCLUDE_DIRECTORIES(${DCMTK_BUNDLE_ROOT})
INCLUDE_DIRECTORIES(${DCMTK_BUNDLE_ROOT}/dcmtk)
INCLUDE_DIRECTORIES(${DCMTK_BUNDLE_ROOT}/dcmtk/config)
INCLUDE_DIRECTORIES(${DCMTK_BUNDLE_ROOT}/dcmtk/dcmdata)
INCLUDE_DIRECTORIES(${DCMTK_BUNDLE_ROOT}/dcmtk/dcmimgle)
INCLUDE_DIRECTORIES(${DCMTK_BUNDLE_ROOT}/dcmtk/dcmjpeg) 
INCLUDE_DIRECTORIES(${DCMTK_BUNDLE_ROOT}/dcmtk/ofstd)

INCLUDE_REGULAR_EXPRESSION("^.*.h$")

#SET(dcmtk_SRCS dcbytstr.cc dcerror.cc  dcobject.cc dctagkey.cc didocu.cc   ofconsol.cc
#               dcdatset.cc dcfilefo.cc dcpcache.cc dctypes.cc  diutils.cc  oflist.cc
#               dcdicent.cc dchashdi.cc dcsequen.cc dcvr.cc     djdecode.cc ofstd.cc
#               dcdict.cc   dcitem.cc   dcstack.cc  dcvrui.cc   djutils.cc  ofstring.cc
#               dcelem.cc   dclist.cc   dctag.cc    dcxfer.cc   ofcond.cc   ofthread.cc )

SET(dcmtk_SRCS cmdlnarg.cc  dcdebug.cc   dcelem.cc    dcistrmf.cc  dcostrma.cc
               dcpxitem.cc  dctagkey.cc  dcvras.cc  dcvrds.cc  dcvrlo.cc    dcvrpobw.cc
               dcvrtm.cc    dcvrut.cc    djcparam.cc  djdecpro.cc  ofconsol.cc  ofstring.cc
               dcbytstr.cc  dcdicdir.cc  dcerror.cc   dcitem.cc    dcostrmf.cc  dcsequen.cc
               dctypes.cc   dcvrat.cc  dcvrdt.cc  dcvrlt.cc    dcvrsh.cc    dcvrui.cc
               dcxfer.cc    djdecbas.cc  djdecsps.cc  ofdate.cc    ofthread.cc
               dcchrstr.cc  dcdicent.cc  dcfilefo.cc  dclist.cc    dcpcache.cc  dcstack.cc
               dcuid.cc     dcvr.cc    dcvrfd.cc  dcvrobow.cc  dcvrsl.cc    dcvrul.cc
               didocu.cc    djdecext.cc  djdecsv1.cc  ofdatime.cc  oftime.cc
               dccodec.cc   dcdict.cc    dchashdi.cc  dcmetinf.cc  dcpixel.cc   dcswap.cc
               dcvm.cc      dcvrcs.cc  dcvrfl.cc  dcvrof.cc    dcvrss.cc    dcvrulup.cc
               diutils.cc   djdeclol.cc  djutils.cc   oflist.cc
               dcdatset.cc  dcdirrec.cc  dcistrma.cc  dcobject.cc  dcpixseq.cc  dctag.cc
               dcvrae.cc    dcvrda.cc  dcvris.cc  dcvrpn.cc    dcvrst.cc    dcvrus.cc
               djcodecd.cc  djdecode.cc  ofcond.cc    ofstd.cc )


ADD_LIBRARY(cmtkDCMTK ${dcmtk_SRCS})
TARGET_LINK_LIBRARIES(cmtkDCMTK)

INSTALL_TARGETS(/lib cmtkDCMTK)
INSTALL_FILES(/include/dcmtk "(\\.h)$")

SET(CMTK_DCMTK_LIBRARIES cmtkDCMTK)



