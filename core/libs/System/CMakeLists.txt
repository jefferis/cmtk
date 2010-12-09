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

# Sources of non-templated classes.
SET(cmtkSystem_SRCS
  cmtkCommandLine.cxx
  cmtkCommandLineCallback.cxx
  cmtkCommandLineItem.cxx
  cmtkCommandLineKeyToAction.cxx
  cmtkCommandLineKeyToActionEnum.cxx
  cmtkCommandLineKeyToActionSingle.cxx
  cmtkCommandLineXML.cxx
  cmtkCommandLineNonOptionParameter.cxx
  cmtkCommandLineNonOptionParameterVector.cxx
  cmtkCompressedStream.cxx
  cmtkCompressedStreamFile.cxx
  cmtkCompressedStreamPipe.cxx
  cmtkCompressedStreamReaderBase.cxx
  cmtkCompressedStreamZlib.cxx
  cmtkConsole.cxx
  cmtkFileUtils.cxx
  cmtkMemory.cxx
  cmtkMountPoints.cxx
  cmtkProgress.cxx	
  cmtkProgressConsole.cxx	
  cmtkStackBacktrace
  cmtkStrUtility.cxx
  cmtkThreads.cxx
  cmtkThreadPool.cxx
  cmtkThreadSemaphore.cxx
  cmtkTimers.cxx
)

IF(CMTK_USE_BZIP2)
  SET(cmtkSystem_SRCS
    ${cmtkSystem_SRCS}
    cmtkCompressedStreamBZip2.cxx
    )
ENDIF(CMTK_USE_BZIP2)

IF(CMTK_USE_LZMA)
  SET(cmtkSystem_SRCS
    ${cmtkSystem_SRCS}
    cmtkCompressedStreamLZMA.cxx
    )
ENDIF(CMTK_USE_LZMA)

ADD_LIBRARY(cmtkSystem ${cmtkSystem_SRCS})

TARGET_LINK_LIBRARIES(cmtkSystem ${CMTK_ZLIB_LIB} ${CMTK_BZIP2_LIBS} ${CMTK_LZMA_LIBS} mxml)

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkSystem PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkSystem
  RUNTIME DESTINATION bin COMPONENT runtime
  LIBRARY DESTINATION lib COMPONENT runtime
  ARCHIVE DESTINATION lib COMPONENT libraries)

FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION include/cmtk/System COMPONENT headers)
