##
##  Copyright 2011 SRI International
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

# ########## Project setup ##########
PROJECT(numdiff-5.2.1)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0005 OLD) ## explicit escaping of compiler definition strings
ENDIF(COMMAND cmake_policy)

SET(CMAKE_CONFIG_PACKAGE_BUGREPORT "<ivprimi at libero dot it>")
SET(CMAKE_CONFIG_PACKAGE_NAME "Numdiff")
SET(CMAKE_CONFIG_PACKAGE_STRING "Numdiff 5.2")
SET(CMAKE_CONFIG_PACKAGE_TARNAME "numdiff")
SET(CMAKE_CONFIG_PACKAGE_VERSION "5.2")

ADD_DEFINITIONS(-DPACKAGE=\\\"numdiff\\\")

# ######### General setup ##########
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# ########## numdiff-5.2.1 executable ##########
# Sources:
SET(numdiff_SRCS
  analyze.c arith.c cmpbuf.c cmpfns.c error.c exitfail.c flags.c getopt.c getopt1.c inout.c io.c main.c numutil.c options.c setmode.c side.c util.c xalloc-die.c xmalloc.c)

INCLUDE(FindCygwin)
SET(POSSIBLE_BIN_PATHS ${CYGWIN_INSTALL_PATH}/bin)

FIND_PROGRAM(FIND_EXECUTABLE NAMES find PATHS ${POSSIBLE_BIN_PATHS})

#-----------------------------------------------------------------------------
# Check for system headers.
INCLUDE (CheckIncludeFiles)

CHECK_INCLUDE_FILES(float.h HAVE_FLOAT_H)
CHECK_INCLUDE_FILES(stdarg.h HAVE_STDARG_H)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)

# see https://www.gnu.org/s/hello/manual/autoconf/Particular-Headers.html
IF(HAVE_FLOAT_H AND HAVE_STDARG_H AND HAVE_STDLIB_H AND HAVE_STRING_H)
  SET(STDC_HEADERS 1)
  ADD_DEFINITIONS(-DSTDC_HEADERS)
ENDIF(HAVE_FLOAT_H AND HAVE_STDARG_H AND HAVE_STDLIB_H AND HAVE_STRING_H)

CHECK_INCLUDE_FILES(dirent.h HAVE_DIRENT_H)
CHECK_INCLUDE_FILES(stat.h HAVE_STAT_H)
CHECK_INCLUDE_FILES(fcntl.h HAVE_FCNTL_H)
CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(libintl.h HAVE_LIBINTL_H)
CHECK_INCLUDE_FILES(locale.h HAVE_LOCALE_H)
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(stdbool.h HAVE_STDBOOL_H)
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES(stddef.h HAVE_STDDEF_H)
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(time.h HAVE_TIME_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(vfork.h HAVE_VFORK_H)
CHECK_INCLUDE_FILES(wchar.h HAVE_WCHAR_H)
CHECK_INCLUDE_FILES(wctype.h HAVE_WCTYPE_H)

CHECK_INCLUDE_FILES (sys/dir.h HAVE_SYS_DIR_H)
CHECK_INCLUDE_FILES (sys/file.h HAVE_SYS_FILE_H)
CHECK_INCLUDE_FILES (sys/file.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES (sys/wait.h HAVE_SYS_WAIT_H)

#AC_C_CONST -- obsolete according to https://www.gnu.org/software/hello/manual/autoconf/C-Compiler.html#index-AC_005fC_005fCONST-885
SET(CMAKE_CONFIG_CONST "const")

#AC_C_INLINE
SET(CMAKE_CONFIG_INLINE "${C_INLINE}")

# Check for variable arrays support in C compiler
INCLUDE(CheckCSourceCompiles)
CHECK_C_SOURCE_COMPILES("int main() { int size = 10; char array[size]; }" HAVE_C_VARARRAYS)

# Check for types and substitute replacement types if they don't exist
INCLUDE(CheckTypeSize)

CHECK_TYPE_SIZE(pid_t PID_T_SIZE)
IF(PID_T_SIZE)
  SET(CMAKE_CONFIG_PID_T pid_t)
ELSE(PID_T_SIZE)
  SET(CMAKE_CONFIG_PID_T int)
ENDIF(PID_T_SIZE)

CHECK_TYPE_SIZE(ptrdiff_t PTRDIFF_T_SIZE)
IF(PTRDIFF_T_SIZE)
  SET(CMAKE_CONFIG_PTRDIFF_T ptrdiff_t)
ELSE(PTRDIFF_T_SIZE)
  SET(CMAKE_CONFIG_PTRDIFF_T int)
ENDIF(PTRDIFF_T_SIZE)

CHECK_TYPE_SIZE(ssize_t SSIZE_T_SIZE)
IF(SSIZE_T_SIZE)
  SET(CMAKE_CONFIG_SSIZE_T ssize_t)
ELSE(SSIZE_T_SIZE)
  SET(CMAKE_CONFIG_SSIZE_T int)
ENDIF(SSIZE_T_SIZE)

CHECK_TYPE_SIZE(size_t SIZE_T_SIZE)
IF(SIZE_T_SIZE)
  SET(CMAKE_CONFIG_SIZE_T size_t)
ELSE(SIZE_T_SIZE)
  SET(CMAKE_CONFIG_SIZE_T unsigned int)
ENDIF(SIZE_T_SIZE)

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMake/CheckStructMember.cmake)
CHECK_STRUCT_MEMBER("struct stat" st_blksize "sys/stat.h" HAVE_STRUCT_STAT_ST_BLKSIZE)
CHECK_STRUCT_MEMBER("struct stat" st_rdev "sys/stat.h" HAVE_STRUCT_STAT_ST_RDEV)
CHECK_STRUCT_MEMBER("struct stat" st_mtim_nsed "sys/stat.h" HAVE_STRUCT_STAT_ST_MTIM_NSEC)

# --------------------------------------
# Checks for libraries
# --------------------------------------

INCLUDE(CheckLibraryExists)
CHECK_INCLUDE_FILES(gmp.h HAVE_GMP_H)
CHECK_LIBRARY_EXISTS(gmp  __gmpz_init "" HAVE_LIBGMP)

# --------------------------------------
# Checks for library functions
# --------------------------------------

INCLUDE(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(vsprintf HAVE_VSPRINTF)
IF(NOT HAVE_VSPRINTF)
  MESSAGE(FATAL_ERROR "Declaration of vsprintf() not present!")
ENDIF(NOT HAVE_VSPRINTF)

CHECK_FUNCTION_EXISTS(memset HAVE_MEMSET)
IF(NOT HAVE_MEMSET)
  MESSAGE(FATAL_ERROR "Declaration of memset() not present!")
ENDIF(NOT HAVE_MEMSET)

CHECK_FUNCTION_EXISTS(memchr HAVE_MEMCHR)
IF(NOT HAVE_MEMCHR)
  MESSAGE(FATAL_ERROR "Declaration of memchr() not present!")
ENDIF(NOT HAVE_MEMCHR)

CHECK_FUNCTION_EXISTS(strerror HAVE_STRERROR)
IF(NOT HAVE_STRERROR)
  MESSAGE(FATAL_ERROR "Declaration of strerror() not present!")
ENDIF(NOT HAVE_STRERROR)

CHECK_FUNCTION_EXISTS(strspn HAVE_STRSPN)
IF(NOT HAVE_STRSPN)
  MESSAGE(FATAL_ERROR "Declaration of strspn() not present!")
ENDIF(NOT HAVE_STRSPN)

CHECK_FUNCTION_EXISTS(strcspn HAVE_STRCSPN)
IF(NOT HAVE_STRCSPN)
  MESSAGE(FATAL_ERROR "Declaration of strcspn() not present!")
ENDIF(NOT HAVE_STRCSPN)

CHECK_FUNCTION_EXISTS(fork HAVE_FORK)
IF(HAVE_FORK)
  SET(HAVE_WORKING_FORK 1)
  SET(CMAKE_CONFIG_VFORK fork)
ENDIF(HAVE_FORK)

CHECK_FUNCTION_EXISTS(strerror_r HAVE_STRERROR_R)
CHECK_FUNCTION_EXISTS(localeconv HAVE_LOCALECONV)
CHECK_FUNCTION_EXISTS(gettimeofday HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(isascii HAVE_ISASCII)
CHECK_FUNCTION_EXISTS(setlocale HAVE_SETLOCALE)
CHECK_FUNCTION_EXISTS(strchr HAVE_STRCHR)
CHECK_FUNCTION_EXISTS(strftime HAVE_STRFTIME)
CHECK_FUNCTION_EXISTS(waitpid HAVE_WAITPID)

CHECK_INCLUDE_FILES (alloca.h HAVE_ALLOCA_H)
CHECK_FUNCTION_EXISTS(alloca HAVE_ALLOCA)

CHECK_FUNCTION_EXISTS(bcopy HAVE_BCOPY)
CHECK_FUNCTION_EXISTS(btowc HAVE_BTOWC)
CHECK_FUNCTION_EXISTS(bzero HAVE_BZERO)
CHECK_FUNCTION_EXISTS(clock_gettime HAVE_CLOCK_GETTIME)

# setmode DOS ?
##AC_FUNC_SETMODE_DOS

## Check for clock_gettime and its library.  Solaris puts it in -lrt or
## -lposix4, but we don't want to link that library unless we have to.
#diff_saved_libs=$LIBS
#AC_SEARCH_LIBS(clock_gettime, [rt posix4])
#case $ac_cv_search_clock_gettime in
#no | 'none required')
#  LIB_CLOCK_GETTIME=;;
#*)
#  LIB_CLOCK_GETTIME=$ac_cv_search_clock_gettime;;
#esac
#AC_SUBST([LIB_CLOCK_GETTIME])


OPTION(NUMDIFF_USE_UNLOCKED_IO "Define if you want getc etc. to use unlocked I/O if available. Unlocked I/O can improve performance in unithreaded apps, but it is not safe for multithreaded apps." ON )

CHECK_FUNCTION_EXISTS(clearerr_unlocked HAVE_DECL_CLEARERR_UNLOCKED)
CHECK_FUNCTION_EXISTS(feof_unlocked HAVE_DECL_FEOF_UNLOCKED)
CHECK_FUNCTION_EXISTS(ferror_unlocked HAVE_DECL_FERROR_UNLOCKED)
CHECK_FUNCTION_EXISTS(fflush_unlocked HAVE_DECL_FFLUSH_UNLOCKED)
CHECK_FUNCTION_EXISTS(fgets_unlocked HAVE_DECL_FGETS_UNLOCKED)
CHECK_FUNCTION_EXISTS(fputc_unlocked HAVE_DECL_FPUTC_UNLOCKED)
CHECK_FUNCTION_EXISTS(fputs_unlocked HAVE_DECL_FPUTS_UNLOCKED)
CHECK_FUNCTION_EXISTS(fread_unlocked HAVE_DECL_FREAD_UNLOCKED)
CHECK_FUNCTION_EXISTS(fwrite_unlocked HAVE_DECL_FWRITE_UNLOCKED)
CHECK_FUNCTION_EXISTS(getc_unlocked HAVE_DECL_GETC_UNLOCKED)
CHECK_FUNCTION_EXISTS(getchar_unlocked HAVE_DECL_GETCHAR_UNLOCKED)
CHECK_FUNCTION_EXISTS(putc_unlocked HAVE_DECL_PUTC_UNLOCKED)
CHECK_FUNCTION_EXISTS(putchar_unlocked HAVE_DECL_PUTCHAR_UNLOCKED)
CHECK_FUNCTION_EXISTS(strerror HAVE_DECL_STRERROR)

IF(HAVE_LIBINTL_H)
  OPTION(NUMDIFF_USE_NLS "Turn on Natural Language Support" ON)

  IF(NUMDIFF_USE_NLS)
    # Check if we have everything we need for NLS:
    # on some systems, libintl is separate, on others part of libc.
    CHECK_FUNCTION_EXISTS(_libintl_gettext HAVE_LIBINTL_GETTEXT_STD)
    IF(NOT HAVE_LIBINTL_GETTEXT_STD)
      FIND_LIBRARY(LIBINTL_LIBRARY intl)
      IF(LIBINTL_LIBRARY_FOUND)
	LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${LIBINTL_LIBRARY})
	CHECK_FUNCTION_EXISTS(_libintl_gettext HAVE_LIBINTL_GETTEXT)
	IF(NOT HAVE_LIBINTL_GETTEXT)
	  SET(NUMDIFF_USE_NLS OFF CACHE "Turn on Natural Language Support")
	ENDIF(NOT HAVE_LIBINTL_GETTEXT)
      ENDIF(LIBINTL_LIBRARY_FOUND)	
    ENDIF(NOT HAVE_LIBINTL_GETTEXT_STD)
  ENDIF(NUMDIFF_USE_NLS)

  IF(NUMDIFF_USE_NLS)
    ADD_DEFINITIONS(-DENABLE_NLS)
  ENDIF(NUMDIFF_USE_NLS)

ENDIF(HAVE_LIBINTL_H)

IF(HAVE_GMP_H AND HAVE_LIBGMP)
  OPTION(NUMDIFF_USE_GMP "Turn on use of GNU MP" ON)
  IF(NUMDIFF_USE_GMP)
    FIND_LIBRARY(LIBGMP_LIBRARY gmp)
  ELSE(NUMDIFF_USE_GMP)
    ADD_DEFINITIONS(-DDISABLE_GMP)
  ENDIF(NUMDIFF_USE_GMP)
ELSE(HAVE_GMP_H AND HAVE_LIBGMP)
  ADD_DEFINITIONS(-DDISABLE_GMP)
ENDIF(HAVE_GMP_H AND HAVE_LIBGMP)

CONFIGURE_FILE(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

# actual target:
ADD_EXECUTABLE(numdiff ${numdiff_SRCS})
TARGET_LINK_LIBRARIES(numdiff ${LIBGMP_LIBRARY} ${LIBINTL_LIBRARY})

# add install target:
INSTALL(TARGETS numdiff DESTINATION bin)
