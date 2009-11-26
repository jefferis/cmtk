PROJECT(cmtkZlib)
INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(zlib_SRCS adler32.c   crc32.c    gzio.c     inffast.c  inftrees.c
	      uncompr.c compress.c  deflate.c  infback.c  inflate.c  trees.c
	      zutil.c)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/zconf.in.h ${CMAKE_CURRENT_BINARY_DIR}/zconf.h)

ADD_LIBRARY(cmtkZlib ${zlib_SRCS})

INSTALL(TARGETS cmtkZlib RUNTIME DESTINATION lib ARCHIVE DESTINATION lib COMPONENT runtime)
INSTALL(TARGETS cmtkZlib LIBRARY DESTINATION lib ARCHIVE DESTINATION lib COMPONENT libraries)
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${files} DESTINATION include/cmtk/zlib COMPONENT headers)

SET(ZLIB_LIBRARY "cmtkZlib" CACHE STRING "Library to link for libz support" FORCE)
