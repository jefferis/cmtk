PROJECT(MiniXML)
INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(MXML_VERSION_MAJOR 2)
SET(MXML_VERSION_MINOR 6)
SET(MXML_VERSION "Mini-XML v2.6")

SET(mxml_SRCS mxml-attr.c  
	      mxml-entity.c
	      mxml-index.c
	      mxml-private.c
	      mxml-set.c
	      mxml-file.c
	      mxml-node.c
	      mxml-search.c
	      mxml-string.c)

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CheckCInline.cmake)

CHECK_TYPE_SIZE("long long" LONG_LONG_SIZE)
IF(LONG_LONG_SIZE)
	SET(HAVE_LONG_LONG 1)
ENDIF(LONG_LONG_SIZE)

INCLUDE(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(strdup HAVE_STRDUP)
CHECK_FUNCTION_EXISTS(snprintf HAVE_SNPRINTF)
CHECK_FUNCTION_EXISTS(vsnprintf HAVE_VSNPRINTF)
CHECK_INCLUDE_FILES(pthread.h HAVE_PTHREAD_H)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_LIBRARY(mxml ${mxml_SRCS})
TARGET_LINK_LIBRARIES(mxml)

INSTALL_TARGETS(/lib mxml)
INSTALL_FILES(/include/mxml "(\\.h)$")

SET(MXML_LIBRARIES mxml)
