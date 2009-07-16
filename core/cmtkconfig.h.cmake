/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#cmakedefine HAVE_DIRENT_H 1
#cmakedefine HAVE_EXECINFO_H 1
#cmakedefine HAVE_FCNTL_H 1
#cmakedefine HAVE_FLOAT_H 1
#cmakedefine HAVE_IEEEFP_H 1
#cmakedefine HAVE_LIMITS_H 1
#cmakedefine HAVE_MALLOC_H 1
#cmakedefine HAVE_PTHREAD_H 1
#cmakedefine HAVE_SIGNAL_H 1
#cmakedefine HAVE_STDARG_H 1
#cmakedefine HAVE_STDDEF_H 1
#cmakedefine HAVE_STDINT_H 1
#cmakedefine HAVE_STDIO_H 1
#cmakedefine HAVE_STDLIB_H 1
#cmakedefine HAVE_STRING_H 1
#cmakedefine HAVE_TIME_H 1
#cmakedefine HAVE_UCONTEXT_H 1
#cmakedefine HAVE_UNISTD_H 1
#cmakedefine HAVE_VALUES_H 1
#cmakedefine HAVE_VARARGS_H 1

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

#ifdef _MSC_VER
// disable warnings about insecure functions (we want to be portable)
#  define _CRT_SECURE_NO_DEPRECATE
// enable POSIX compliance
#  define _POSIX_
#  define NOMINMAX
#  define snprintf _snprintf
#  define strdup _strdup
#  define isnan(a) (!_finite(a))
#  define finite(a) _finite(a)
#  define random rand
#  define srandom srand
#endif //#ifdef _MSC_VER

#endif // #ifndef __cmtkconfig_h_included__
