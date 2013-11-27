/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

/* cmtkconfig.h.cmake  */

#ifndef __cmtkconfig_h_included__
#define __cmtkconfig_h_included__

#define CMTK_VERSION_MAJOR @CMTK_VERSION_MAJOR@
#define CMTK_VERSION_MINOR @CMTK_VERSION_MINOR@
#define CMTK_VERSION_PATCH "@CMTK_VERSION_PATCH@"
#define CMTK_VERSION_STRING "@CMTK_VERSION_STRING@"

#define CMTK_ROOT_PATH_SRI24 "@CMTK_ROOT_PATH_SRI24@"

#define CMTK_BINARY_DIR "@CMTK_BINARY_DIR_CONFIG@"
#define CMTK_DATADIR "@CMTK_DATA_ROOT_CONFIG@/testing/inputs"

#define CMTK_DCMDICTPATH "@CMTK_DCMDICTPATH_CONFIG@"
#define CMTK_DCMDICTPATH_INSTALL "@CMTK_DCMDICTPATH_INSTALL_CONFIG@"

// Unless in "DEBUG" build, turn off AlgLib assertions
#ifndef DEBUG
#define NO_AP_ASSERT 1
#endif

//
// Configuration options
//

#cmakedefine CMTK_BUILD_UNSTABLE 1
#cmakedefine CMTK_BUILD_STACKTRACE 1
#cmakedefine CMTK_BUILD_RTRACKER 1
#cmakedefine CMTK_BUILD_DEMO 1
#cmakedefine CMTK_BUILD_NRRD 1

#cmakedefine CMTK_USE_SMP 1
#cmakedefine CMTK_USE_PTHREADS 1
#cmakedefine CMTK_USE_CUDA 1
#cmakedefine CMTK_USE_BZIP2 1
#cmakedefine CMTK_USE_LZMA 1
#cmakedefine CMTK_USE_DCMTK 1
#cmakedefine CMTK_USE_SQLITE 1

#cmakedefine CMTK_USE_FFTW_FOUND 1

#cmakedefine CMTK_COORDINATES_DOUBLE 1
#ifndef CMTK_COORDINATES_DOUBLE
#  define CMTK_COORDINATES_FLOAT 1
#endif

#cmakedefine CMTK_DATA_DOUBLE 1
#ifndef CMTK_DATA_DOUBLE
#  define CMTK_DATA_FLOAT 1
#endif

#cmakedefine CMTK_NUMERICS_DOUBLE 1
#ifndef CMTK_NUMERICS_DOUBLE
#  define CMTK_NUMERICS_FLOAT 1
#endif

#cmakedefine CMTK_COMPILER_VAR_AUTO_ARRAYSIZE 1

#cmakedefine HAVE_DIRENT_H 1
#cmakedefine HAVE_EXECINFO_H 1
#cmakedefine HAVE_FCNTL_H 1
#cmakedefine HAVE_IEEEFP_H 1
#cmakedefine HAVE_INTTYPES_H 1
#cmakedefine HAVE_MALLOC_H 1
#cmakedefine HAVE_PTHREAD_H 1
#cmakedefine HAVE_STDINT_H 1
#cmakedefine HAVE_TERMIOS_H 1
#cmakedefine HAVE_UNISTD_H 1
#cmakedefine HAVE_VALUES_H 1

#cmakedefine HAVE_SYS_IOCTL_H 1
#cmakedefine HAVE_SYS_PROCFS_H 1
#cmakedefine HAVE_SYS_STAT_H 1
#cmakedefine HAVE_SYS_TIMES_H 1
#cmakedefine HAVE_SYS_TIME_H 1
#cmakedefine HAVE_SYS_TYPES_H 1
#cmakedefine HAVE_SYS_UTSNAME_H 1

#cmakedefine HAVE_HASH_MAP 1
#cmakedefine HAVE_HASH_MAP_H 1

#cmakedefine HAVE_UNORDERED_MAP 1
#cmakedefine HAVE_UNORDERED_MAP_TR1 1

#cmakedefine WORDS_BIGENDIAN 1

/* Use stat64 on systems where it is available and stat is not 64bit aware. */
#cmakedefine CMTK_USE_STAT64 1

// Flag for Grand Central Dispatch
#cmakedefine CMTK_USE_GCD 1

/* The size of a `char', as computed by sizeof. */
#cmakedefine SIZEOF_CHAR @CMAKE_SIZEOF_CHAR@

/* The size of a `double', as computed by sizeof. */
#cmakedefine SIZEOF_DOUBLE @CMAKE_SIZEOF_DOUBLE@

/* The size of a `float', as computed by sizeof. */
#cmakedefine SIZEOF_FLOAT @CMAKE_SIZEOF_FLOAT@

/* The size of a `int', as computed by sizeof. */
#cmakedefine SIZEOF_INT @CMAKE_SIZEOF_INT@

/* The size of a `long', as computed by sizeof. */
#cmakedefine SIZEOF_LONG @CMAKE_SIZEOF_LONG@

/* The size of a `short', as computed by sizeof. */
#cmakedefine SIZEOF_SHORT @CMAKE_SIZEOF_SHORT@

/* The size of a `void *', as computed by sizeof. */
#cmakedefine SIZEOF_VOID_P @CMAKE_SIZEOF_VOID_P@

/// Macro to prevent warnings from unused function arguments.
#define UNUSED(a) ((void)a)

#ifdef _MSC_VER
// disable warnings about insecure functions (we want to be portable)
#  define _CRT_SECURE_NO_DEPRECATE
// disable warnings about unknown (i.e., gcc) pragmas
#  pragma warning ( disable: 4068 )
// disable warnings about unimplemented "throw()" declaration
#pragma warning(disable: 4290)
// enable POSIX compliance
#  define _POSIX_
#  define NOMINMAX
#  define snprintf _snprintf
#  define strdup _strdup
#  define isnan(a) (!_finite(a))
#  define finite(a) _finite(a)
#  define random rand
#  define srandom srand

#  define CMTK_PATH_SEPARATOR '\\'
#  define CMTK_PATH_SEPARATOR_STR "\\"

#else

#  define CMTK_PATH_SEPARATOR '/'
#  define CMTK_PATH_SEPARATOR_STR "/"

#endif //#ifdef _MSC_VER

#endif // #ifndef __cmtkconfig_h_included__
