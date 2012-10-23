# 
# NrrdIO: stand-alone code for basic nrrd functionality
# Copyright (C) 2011, 2010, 2009  University of Chicago
# Copyright (C) 2008, 2007, 2006, 2005  Gordon Kindlmann
# Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah
#
# This software is provided 'as-is', without any express or implied
# warranty.  In no event will the authors be held liable for any
# damages arising from the use of this software.
# 
# Permission is granted to anyone to use this software for any
# purpose, including commercial applications, and to alter it and
# redistribute it freely, subject to the following restrictions:
# 
# 1. The origin of this software must not be misrepresented; you must
#    not claim that you wrote the original software. If you use this
#    software in a product, an acknowledgment in the product
#    documentation would be appreciated but is not required.
#
#   2. Altered source versions must be plainly marked as such, and must
#    not be misrepresented as being the original software.
# 
# 3. This notice may not be removed or altered from any source distribution.
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
PROJECT(NrrdIO)
INCLUDE_REGULAR_EXPRESSION("^.*.h$")

#
# This CMake file configures the NrrdIO library build.  NrrdIO
# is used by Insight/Code/IO/itkNrrdIO for reading/writing 
# "Nearly Raw Raster Data" within the open-source Teem software 
# package. See http://teem.sourceforge.net for more information.
#

SET(nrrdio_SRCS comment.c enumsNrrd.c mop.c string.c 754.c defaultsNrrd.c
parseAir.c dio.c format.c parseNrrd.c formatEPS.c encoding.c
formatNRRD.c encodingAscii.c formatPNG.c encodingBzip2.c
formatPNM.c accessors.c encodingGzip.c formatText.c 
array.c encodingHex.c formatVTK.c read.c arraysNrrd.c encodingRaw.c 
gzio.c reorder.c write.c axis.c endianAir.c  keyvalue.c
biffbiff.c biffmsg.c endianNrrd.c methodsNrrd.c  sane.c enum.c
miscAir.c simple.c )

# Turn on TEEM_BUILD so that the proper dll export def's are 
# used on windows builds.
ADD_DEFINITIONS(-DTEEM_BUILD=1)

#The QNANHIBIT variable is configured by the root level CMakeLists.txt
IF(QNANHIBIT)
  ADD_DEFINITIONS(-DTEEM_QNANHIBIT=1)
ELSE(QNANHIBIT)
  ADD_DEFINITIONS(-DTEEM_QNANHIBIT=0)
ENDIF(QNANHIBIT)

#DirectIO is the fast way to do multi-gigabyte I/O and currently only available
#for SGI platforms.  Use of DirectIO is enabled manually for now.
#OPTION(USE_DIRECTIO "Use DirectIO for Nrrd file IO.  Only valid on SGI systems." 0)
#MARK_AS_ADVANCED(USE_DIRECTIO)
#IF(USE_DIRECTIO)
#  ADD_DEFINITIONS(-DTEEM_DIO=1)
#ELSE(USE_DIRECTIO)
  ADD_DEFINITIONS(-DTEEM_DIO=0)
#ENDIF(USE_DIRECTIO)

# Possibly turn on usage of zlib compression (requires linking with libz)
# (i.e., programs compiled with ITKNrrdIO must also be compiled with zlib)
#ADD_DEFINITIONS(-DTEEM_ZLIB=1)

ADD_LIBRARY(NrrdIO ${nrrdio_SRCS} )

## These are ITK-specific
#TARGET_LINK_LIBRARIES(NrrdIO ${ITK_ZLIB_LIBRARIES} )
#INSTALL_TARGETS(/lib/InsightToolkit ITKNrrdIO)
#INSTALL_FILES(/include/InsightToolkit/Utilities/NrrdIO "(\\.h)$")

