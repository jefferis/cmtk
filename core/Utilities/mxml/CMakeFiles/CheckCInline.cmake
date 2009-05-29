# This handy test is from Jack Kelly on the cmake email list. 
#   http://www.cmake.org/Wiki/CMakeTestInline

# Inspired from /usr/share/autoconf/autoconf/c.m4
FOREACH(KEYWORD "inline" "__inline__" "__inline")
   IF(NOT DEFINED C_INLINE)
     TRY_COMPILE(C_HAS_KEYWORD "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/CheckCInline.c" COMPILE_DEFINITIONS "-Dinline=${KEYWORD}")
     IF(C_HAS_KEYWORD)
       SET(C_INLINE ${KEYWORD})
     ENDIF(C_HAS_KEYWORD)
   ENDIF(NOT DEFINED C_INLINE)
ENDFOREACH(KEYWORD)
