PROJECT(cmtkZlib)
INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(zlib_SRCS adler32.c   crc32.c    gzio.c     inffast.c  inftrees.c
	      uncompr.c compress.c  deflate.c  infback.c  inflate.c  trees.c
	      zutil.c)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/zconf.in.h ${CMAKE_CURRENT_BINARY_DIR}/zconf.h)

ADD_LIBRARY(cmtkZlib ${zlib_SRCS})
TARGET_LINK_LIBRARIES(cmtkZlib)

INSTALL_TARGETS(/lib cmtkZlib)
INSTALL_FILES(/include/zlib "(\\.h)$")

SET(ZLIB_LIBRARIES cmtkZlib)
