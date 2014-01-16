PROJECT(D4CMTK)

##
##  Copyright 2004-2012, 2014 SRI International
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

SET(DCMTK_PREFIX ${CMAKE_BINARY_DIR})
SET(DCM_DICT_DEFAULT_PATH "/dicom.dic;${DCMTK_PREFIX}/lib/dicom.dic;${CMTK_INSTALL_DATA_DIR}/dicom.dic")

SET(WITH_ZLIB 1)

IF(CMTK_USE_SMP)
  SET(WITH_THREADS 1)
ENDIF(CMTK_USE_SMP)

INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityC.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityCXX.cmake)

CHECK_TYPE_SIZE(int64_t HAVE_INT64_T)
CHECK_TYPE_SIZE(uint64_t HAVE_UINT64_T)

CHECK_INCLUDE_FILE_CXX (sstream HAVE_SSTREAM)
CHECK_INCLUDE_FILE_CXX (strstream HAVE_STRSTREAM)
CHECK_INCLUDE_FILE_CXX (sstream.h HAVE_SSTREAM_H)
CHECK_INCLUDE_FILE_CXX (strstrea.h HAVE_STRSTREA_H)
CHECK_INCLUDE_FILE_CXX (strstream.h HAVE_STRSTREAM_H)
CHECK_INCLUDE_FILE_CXX (alloca.h HAVE_ALLOCA_H)
CHECK_INCLUDE_FILE_CXX (arpa/inet.h HAVE_ARPA_INET_H)
CHECK_INCLUDE_FILE_CXX (assert.h HAVE_ASSERT_H)
CHECK_INCLUDE_FILE_CXX (ctype.h HAVE_CTYPE_H)
CHECK_INCLUDE_FILE_CXX (dirent.h HAVE_DIRENT_H)
CHECK_INCLUDE_FILE_CXX (errno.h HAVE_ERRNO_H)
CHECK_INCLUDE_FILE_CXX (fcntl.h HAVE_FCNTL_H)
CHECK_INCLUDE_FILE_CXX (float.h HAVE_FLOAT_H)
CHECK_INCLUDE_FILE_CXX (fnmatch.h HAVE_FNMATCH_H)
CHECK_INCLUDE_FILE_CXX (fstream HAVE_FSTREAM)
CHECK_INCLUDE_FILE_CXX (fstream.h HAVE_FSTREAM_H)
CHECK_INCLUDE_FILE_CXX (grp.h HAVE_GRP_H)
CHECK_INCLUDE_FILE_CXX (ieeefp.h HAVE_IEEEFP_H)
CHECK_INCLUDE_FILE_CXX (inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE_CXX (iomanip HAVE_IOMANIP)
CHECK_INCLUDE_FILE_CXX (iomanip.h HAVE_IOMANIP_H)
CHECK_INCLUDE_FILE_CXX (iostream HAVE_IOSTREAM)
CHECK_INCLUDE_FILE_CXX (iostream.h HAVE_IOSTREAM_H)
CHECK_INCLUDE_FILE_CXX (io.h HAVE_IO_H)
CHECK_INCLUDE_FILE_CXX (iso646.h HAVE_ISO646_H)
CHECK_INCLUDE_FILE_CXX (libc.h HAVE_LIBC_H)
IF(HAVE_LIBC_H)
  SET(INCLUDE_LIBC_H_AS_CXX 1)
ELSE(HAVE_LIBC_H)
  CHECK_INCLUDE_FILE (libc.h HAVE_LIBC_H)
ENDIF(HAVE_LIBC_H)

CHECK_INCLUDE_FILE_CXX (libpng/png.h HAVE_LIBPNG_PNG_H)
CHECK_INCLUDE_FILE_CXX (limits.h HAVE_LIMITS_H)
CHECK_INCLUDE_FILE_CXX (locale.h HAVE_LOCALE_H)
CHECK_INCLUDE_FILE_CXX (malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILE_CXX (math.h HAVE_MATH_H)
CHECK_INCLUDE_FILE_CXX (memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILE_CXX (netdb.h HAVE_NETDB_H)
CHECK_INCLUDE_FILE_CXX (netinet/in.h HAVE_NETINET_IN_H)
CHECK_INCLUDE_FILE_CXX (netinet/in_systm.h HAVE_NETINET_IN_SYSTM_H)
CHECK_INCLUDE_FILE_CXX (netinet/tcp.h HAVE_NETINET_TCP_H)
CHECK_INCLUDE_FILE_CXX (new.h HAVE_NEW_H)
CHECK_INCLUDE_FILE_CXX (pthread.h HAVE_PTHREAD_H)
CHECK_INCLUDE_FILE_CXX (pwd.h HAVE_PWD_H)
CHECK_INCLUDE_FILE_CXX (semaphore.h HAVE_SEMAPHORE_H)
CHECK_INCLUDE_FILE_CXX (setjmp.h HAVE_SETJMP_H)
CHECK_INCLUDE_FILE_CXX (signal.h HAVE_SIGNAL_H)
CHECK_INCLUDE_FILE_CXX (sstream.h HAVE_SSTREAM_H)
CHECK_INCLUDE_FILE_CXX (stat.h HAVE_STAT_H)
CHECK_INCLUDE_FILE_CXX (stdarg.h HAVE_STDARG_H)
CHECK_INCLUDE_FILE_CXX (stdbool.h HAVE_STDBOOL_H)
CHECK_INCLUDE_FILE_CXX (stddef.h HAVE_STDDEF_H)
CHECK_INCLUDE_FILE_CXX (stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILE_CXX (stdio.h HAVE_STDIO_H)
CHECK_INCLUDE_FILE_CXX (stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILE_CXX (strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILE_CXX (string.h HAVE_STRING_H)
CHECK_INCLUDE_FILE_CXX (synch.h HAVE_SYNCH_H)
CHECK_INCLUDE_FILE_CXX (sys/errno.h HAVE_SYS_ERRNO_H)
CHECK_INCLUDE_FILE_CXX (sys/file.h HAVE_SYS_FILE_H)
CHECK_INCLUDE_FILE_CXX (sys/param.h HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILE_CXX (sys/resource.h HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILE_CXX (sys/select.h HAVE_SYS_SELECT_H)
CHECK_INCLUDE_FILE_CXX (sys/socket.h HAVE_SYS_SOCKET_H)
CHECK_INCLUDE_FILE_CXX (sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE_CXX (sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE_CXX (sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE_CXX (sys/utime.h HAVE_SYS_UTIME_H)
CHECK_INCLUDE_FILE_CXX (sys/utsname.h HAVE_SYS_UTSNAME_H)
CHECK_INCLUDE_FILE_CXX (sys/wait.h HAVE_SYS_WAIT_H)
CHECK_INCLUDE_FILE_CXX (thread.h HAVE_THREAD_H)
CHECK_INCLUDE_FILE_CXX (time.h HAVE_TIME_H)
CHECK_INCLUDE_FILE_CXX (unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILE_CXX (unix.h HAVE_UNIX_H)
CHECK_INCLUDE_FILE_CXX (utime.h HAVE_UTIME_H)
CHECK_INCLUDE_FILE_CXX (wctype.h HAVE_WCTYPE_H)

INCLUDE(CheckFunctionExists)
CHECK_FUNCTION_EXISTS (accept HAVE_ACCEPT)
CHECK_FUNCTION_EXISTS (access HAVE_ACCESS)
CHECK_FUNCTION_EXISTS (bcmp HAVE_BCMP)
CHECK_FUNCTION_EXISTS (bcopy HAVE_BCOPY)
CHECK_FUNCTION_EXISTS (bind HAVE_BIND)
CHECK_FUNCTION_EXISTS (bzero HAVE_BZERO)
CHECK_FUNCTION_EXISTS (connect HAVE_CONNECT)
CHECK_FUNCTION_EXISTS (cuserid HAVE_CUSERID)
CHECK_FUNCTION_EXISTS (finite HAVE_FINITE)
CHECK_FUNCTION_EXISTS (flock HAVE_FLOCK)
CHECK_FUNCTION_EXISTS (fork HAVE_FORK)
CHECK_FUNCTION_EXISTS (getenv HAVE_GETENV)
CHECK_FUNCTION_EXISTS (geteuid HAVE_GETEUID)
CHECK_FUNCTION_EXISTS (getgrnam HAVE_GETGRNAM)
CHECK_FUNCTION_EXISTS (gethostbyname HAVE_GETHOSTBYNAME)
CHECK_FUNCTION_EXISTS (gethostid HAVE_GETHOSTID)
CHECK_FUNCTION_EXISTS (gethostname HAVE_GETHOSTNAME)
CHECK_FUNCTION_EXISTS (getlogin HAVE_GETLOGIN)
CHECK_FUNCTION_EXISTS (getpid HAVE_GETPID)
CHECK_FUNCTION_EXISTS (getpwnam HAVE_GETPWNAM)
CHECK_FUNCTION_EXISTS (getsockname HAVE_GETSOCKNAME)
CHECK_FUNCTION_EXISTS (getsockopt HAVE_GETSOCKOPT)
CHECK_FUNCTION_EXISTS (getuid HAVE_GETUID)
CHECK_FUNCTION_EXISTS (index HAVE_INDEX)
CHECK_FUNCTION_EXISTS (isinf HAVE_ISINF)
CHECK_FUNCTION_EXISTS (isnan HAVE_ISNAN)
CHECK_FUNCTION_EXISTS (itoa HAVE_ITOA)
CHECK_FUNCTION_EXISTS (listen HAVE_LISTEN)
CHECK_FUNCTION_EXISTS (lockf HAVE_LOCKF)
CHECK_FUNCTION_EXISTS (malloc_debug HAVE_MALLOC_DEBUG)
CHECK_FUNCTION_EXISTS (memcmp HAVE_MEMCMP)
CHECK_FUNCTION_EXISTS (memcpy HAVE_MEMCPY)
CHECK_FUNCTION_EXISTS (memmove HAVE_MEMMOVE)
CHECK_FUNCTION_EXISTS (memset HAVE_MEMSET)
CHECK_FUNCTION_EXISTS (mkstemp HAVE_MKSTEMP)
CHECK_FUNCTION_EXISTS (mktemp HAVE_MKTEMP)
CHECK_FUNCTION_EXISTS (rindex HAVE_RINDEX)
CHECK_FUNCTION_EXISTS (select HAVE_SELECT)
CHECK_FUNCTION_EXISTS (setsockopt HAVE_SETSOCKOPT)
CHECK_FUNCTION_EXISTS (setuid HAVE_SETUID)
CHECK_FUNCTION_EXISTS (sleep HAVE_SLEEP)
CHECK_FUNCTION_EXISTS (socket HAVE_SOCKET)
CHECK_FUNCTION_EXISTS (stat HAVE_STAT)
CHECK_FUNCTION_EXISTS (strchr HAVE_STRCHR)
CHECK_FUNCTION_EXISTS (strdup HAVE_STRDUP)
CHECK_FUNCTION_EXISTS (strerror HAVE_STRERROR)
CHECK_FUNCTION_EXISTS (strlcat HAVE_STRLCAT)
CHECK_FUNCTION_EXISTS (strlcpy HAVE_STRLCPY)
CHECK_FUNCTION_EXISTS (strstr HAVE_STRSTR)
CHECK_FUNCTION_EXISTS (strtoul HAVE_STRTOUL)
CHECK_FUNCTION_EXISTS (sysinfo HAVE_SYSINFO)
CHECK_FUNCTION_EXISTS (tempnam HAVE_TEMPNAM)
CHECK_FUNCTION_EXISTS (tmpnam HAVE_TMPNAM)
CHECK_FUNCTION_EXISTS (uname HAVE_UNAME)
CHECK_FUNCTION_EXISTS (usleep HAVE_USLEEP)
CHECK_FUNCTION_EXISTS (vprintf HAVE_VPRINTF)
CHECK_FUNCTION_EXISTS (waitpid HAVE_WAITPID)
CHECK_FUNCTION_EXISTS (_findfirst HAVE__FINDFIRST)

CHECK_FUNCTION_EXISTS(accept HAVE_PROTOTYPE_ACCEPT)
CHECK_FUNCTION_EXISTS(bind HAVE_PROTOTYPE_BIND)
CHECK_FUNCTION_EXISTS(bzero HAVE_PROTOTYPE_BZERO)
CHECK_FUNCTION_EXISTS(connect HAVE_PROTOTYPE_CONNECT)
CHECK_FUNCTION_EXISTS(finite HAVE_PROTOTYPE_FINITE)
CHECK_FUNCTION_EXISTS(flock HAVE_PROTOTYPE_FLOCK)
CHECK_FUNCTION_EXISTS(gethostbyname HAVE_PROTOTYPE_GETHOSTBYNAME)
CHECK_FUNCTION_EXISTS(gethostid HAVE_PROTOTYPE_GETHOSTID)
CHECK_FUNCTION_EXISTS(gethostname HAVE_PROTOTYPE_GETHOSTNAME)
CHECK_FUNCTION_EXISTS(getsockname HAVE_PROTOTYPE_GETSOCKNAME)
CHECK_FUNCTION_EXISTS(getsockopt HAVE_PROTOTYPE_GETSOCKOPT)
CHECK_FUNCTION_EXISTS(gettimeofday HAVE_PROTOTYPE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(isinf HAVE_PROTOTYPE_ISINF)
CHECK_FUNCTION_EXISTS(isnan HAVE_PROTOTYPE_ISNAN)
CHECK_FUNCTION_EXISTS(listen HAVE_PROTOTYPE_LISTEN)
CHECK_FUNCTION_EXISTS(mkstemp HAVE_PROTOTYPE_MKSTEMP)
CHECK_FUNCTION_EXISTS(mktemp HAVE_PROTOTYPE_MKTEMP)
CHECK_FUNCTION_EXISTS(select HAVE_PROTOTYPE_SELECT)
CHECK_FUNCTION_EXISTS(setsockopt HAVE_PROTOTYPE_SETSOCKOPT)
CHECK_FUNCTION_EXISTS(socket HAVE_PROTOTYPE_SOCKET)
CHECK_FUNCTION_EXISTS(std::_vfprintf HAVE_PROTOTYPE_STD__VFPRINTF)
CHECK_FUNCTION_EXISTS(strcasecmp HAVE_PROTOTYPE_STRCASECMP)
CHECK_FUNCTION_EXISTS(stdncasecmp HAVE_PROTOTYPE_STRNCASECMP)
CHECK_FUNCTION_EXISTS(usleep HAVE_PROTOTYPE_USLEEP)
CHECK_FUNCTION_EXISTS(wait3 HAVE_PROTOTYPE_WAIT3)
CHECK_FUNCTION_EXISTS(waitpid HAVE_PROTOTYPE_WAITPID)
CHECK_FUNCTION_EXISTS(_stricmp HAVE_PROTOTYPE__STRICMP)

#
# The following from ITK's DCMTK git repository
#
INCLUDE(CMake/CheckFunctionWithHeaderExists.cmake)
INCLUDE(CMake/dcmtkTryCompile.cmake)

  SET(HEADERS)

  IF(HAVE_IO_H)
    SET(HEADERS ${HEADERS} io.h)
  ENDIF(HAVE_IO_H)

  IF(HAVE_IOSTREAM_H)
    SET(HEADERS ${HEADERS} iostream.h)
  ENDIF(HAVE_IOSTREAM_H)

  IF(HAVE_MATH_H)
    SET(HEADERS ${HEADERS} math.h)
  ENDIF(HAVE_MATH_H)

  IF(HAVE_LIBC_H)
    SET(HEADERS ${HEADERS} libc.h)
  ENDIF(HAVE_LIBC_H)

  IF(HAVE_PTHREAD_H)
    SET(HEADERS ${HEADERS} pthread.h)
  ENDIF(HAVE_PTHREAD_H)

  IF(HAVE_UNISTD_H)
    SET(HEADERS ${HEADERS} unistd.h)
  ENDIF(HAVE_UNISTD_H)

  IF(HAVE_STDLIB_H)
    SET(HEADERS ${HEADERS} stdlib.h)
  ENDIF(HAVE_STDLIB_H)

  IF(HAVE_NETDB_H)
    SET(HEADERS ${HEADERS} netdb.h)
  ENDIF(HAVE_NETDB_H)

  IF(HAVE_SYS_FILE_H)
    SET(HEADERS ${HEADERS} sys/file.h)
  ENDIF(HAVE_SYS_FILE_H)

  IF(HAVE_STRING_H)
    SET(HEADERS ${HEADERS} string.h)
  ENDIF(HAVE_STRING_H)

  IF(HAVE_STRINGS_H)
    SET(HEADERS ${HEADERS} strings.h)
  ENDIF(HAVE_STRINGS_H)

  IF(HAVE_SYS_WAIT_H)
    SET(HEADERS ${HEADERS} sys/wait.h)
  ENDIF(HAVE_SYS_WAIT_H)

  IF(HAVE_SYS_TIME_H)
    SET(HEADERS ${HEADERS} sys/time.h)
  ENDIF(HAVE_SYS_TIME_H)

  IF(HAVE_SYS_RESOURCE_H)
    SET(HEADERS ${HEADERS} sys/resource.h)
  ENDIF(HAVE_SYS_RESOURCE_H)

  IF(HAVE_SYS_TYPES_H)
    SET(HEADERS ${HEADERS} sys/types.h)
  ENDIF(HAVE_SYS_TYPES_H)

  IF(HAVE_SYS_SOCKET_H)
    SET(HEADERS ${HEADERS} sys/socket.h)
  ENDIF(HAVE_SYS_SOCKET_H)

  IF(HAVE_STDARG_H)
    SET(HEADERS ${HEADERS} stdarg.h)
  ENDIF(HAVE_STDARG_H)

  IF(HAVE_STDIO_H)
    SET(HEADERS ${HEADERS} stdio.h)
  ENDIF(HAVE_STDIO_H)

  IF(HAVE_SYS_SELECT_H)
    SET(HEADERS ${HEADERS} sys/select.h)
  ENDIF(HAVE_SYS_SELECT_H)

  IF(HAVE_WINDOWS_H)
    SET(HEADERS ${HEADERS} windows.h)
  ENDIF(HAVE_WINDOWS_H)

  IF(HAVE_WINSOCK_H)
    SET(HEADERS ${HEADERS} winsock.h)
    SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ws2_32 netapi32 wsock32)
  ENDIF(HAVE_WINSOCK_H)

  # std::vsnprintf and std::vsnprintf need the C++ version of the headers.
  # We just assume they exist when the C version was found
  SET(CXXHEADERS)

  IF(HAVE_CSTDIO)
    SET(CXXHEADERS ${CXXHEADERS} cstdio)
  ENDIF(HAVE_CSTDIO)

  IF(HAVE_CSTDARG)
    SET(CXXHEADERS ${CXXHEADERS} cstdarg)
  ENDIF(HAVE_CSTDARG)

  CHECK_FUNCTIONWITHHEADER_EXISTS(isinf "${HEADERS}" HAVE_PROTOTYPE_ISINF)
  CHECK_FUNCTIONWITHHEADER_EXISTS(isnan "${HEADERS}" HAVE_PROTOTYPE_ISNAN)
  CHECK_FUNCTIONWITHHEADER_EXISTS(finite "${HEADERS}" HAVE_PROTOTYPE_FINITE)
  CHECK_FUNCTIONWITHHEADER_EXISTS(flock "${HEADERS}" HAVE_PROTOTYPE_FLOCK)
  CHECK_FUNCTIONWITHHEADER_EXISTS(gethostbyname "${HEADERS}" HAVE_PROTOTYPE_GETHOSTBYNAME)
  CHECK_FUNCTIONWITHHEADER_EXISTS(gethostbyname_r "${HEADERS}" HAVE_PROTOTYPE_GETHOSTBYNAME_R)
  CHECK_FUNCTIONWITHHEADER_EXISTS(gethostid "${HEADERS}" HAVE_PROTOTYPE_GETHOSTID)
  CHECK_FUNCTIONWITHHEADER_EXISTS(bzero "${HEADERS}" HAVE_PROTOTYPE_BZERO)
  CHECK_FUNCTIONWITHHEADER_EXISTS(gethostname "${HEADERS}" HAVE_PROTOTYPE_GETHOSTNAME)
  CHECK_FUNCTIONWITHHEADER_EXISTS(waitpid "${HEADERS}" HAVE_PROTOTYPE_WAITPID)
  CHECK_FUNCTIONWITHHEADER_EXISTS(wait3 "${HEADERS}" HAVE_PROTOTYPE_WAIT3)
  CHECK_FUNCTIONWITHHEADER_EXISTS(usleep "${HEADERS}" HAVE_PROTOTYPE_USLEEP)
  CHECK_FUNCTIONWITHHEADER_EXISTS(accept "${HEADERS}" HAVE_PROTOTYPE_ACCEPT)
  CHECK_FUNCTIONWITHHEADER_EXISTS(bind "${HEADERS}" HAVE_PROTOTYPE_BIND)
  CHECK_FUNCTIONWITHHEADER_EXISTS(getsockname "${HEADERS}" HAVE_PROTOTYPE_GETSOCKNAME)
  CHECK_FUNCTIONWITHHEADER_EXISTS(getsockopt "${HEADERS}" HAVE_PROTOTYPE_GETSOCKOPT)
  CHECK_FUNCTIONWITHHEADER_EXISTS(getsockname "${HEADERS}" HAVE_PROTOTYPE_GETSOCKNAME)
  CHECK_FUNCTIONWITHHEADER_EXISTS(setsockopt "${HEADERS}" HAVE_PROTOTYPE_SETSOCKOPT)
  CHECK_FUNCTIONWITHHEADER_EXISTS(socket "${HEADERS}" HAVE_PROTOTYPE_SOCKET)
  CHECK_FUNCTIONWITHHEADER_EXISTS(listen "${HEADERS}" HAVE_PROTOTYPE_LISTEN)
  CHECK_FUNCTIONWITHHEADER_EXISTS(std::vfprintf "${CXXHEADERS}" HAVE_PROTOTYPE_STD__VFPRINTF)
  CHECK_FUNCTIONWITHHEADER_EXISTS(std::vsnprintf "${CXXHEADERS}" HAVE_PROTOTYPE_STD__VSNPRINTF)
  CHECK_FUNCTIONWITHHEADER_EXISTS(_stricmp "${HEADERS}" HAVE_PROTOTYPE__STRICMP)
  # "definition" is an (exchangeable) identifier that is needed for successful compile test
  CHECK_FUNCTIONWITHHEADER_EXISTS("socklen_t definition" "${HEADERS}" HAVE_DECLARATION_SOCKLEN_T)
  CHECK_FUNCTIONWITHHEADER_EXISTS(gettimeofday "${HEADERS}" HAVE_PROTOTYPE_GETTIMEOFDAY)
  CHECK_FUNCTIONWITHHEADER_EXISTS(connect "${HEADERS}" HAVE_PROTOTYPE_CONNECT)
  CHECK_FUNCTIONWITHHEADER_EXISTS(mkstemp "${HEADERS}" HAVE_PROTOTYPE_MKSTEMP)
  CHECK_FUNCTIONWITHHEADER_EXISTS(mktemp "${HEADERS}" HAVE_PROTOTYPE_MKTEMP)
  CHECK_FUNCTIONWITHHEADER_EXISTS(select "${HEADERS}" HAVE_PROTOTYPE_SELECT)
  CHECK_FUNCTIONWITHHEADER_EXISTS(strcasecmp "${HEADERS}" HAVE_PROTOTYPE_STRCASECMP)
  CHECK_FUNCTIONWITHHEADER_EXISTS(strncasecmp "${HEADERS}" HAVE_PROTOTYPE_STRNCASECMP)
  CHECK_FUNCTIONWITHHEADER_EXISTS(strerror_r "${HEADERS}" HAVE_PROTOTYPE_STRERROR_R)
  # "definition" is an (exchangeable) identifier that is needed for successful compile test
  CHECK_FUNCTIONWITHHEADER_EXISTS("std::ios_base::openmode definition" "${HEADERS}" HAVE_DECLARATION_STD__IOS_BASE__OPENMODE)
  CHECK_FUNCTIONWITHHEADER_EXISTS(pthread_rwlock_init "${HEADERS}" HAVE_PTHREAD_RWLOCK)
  CHECK_FUNCTIONWITHHEADER_EXISTS("__sync_add_and_fetch((int*)0,0)" "${HEADERS}" HAVE_SYNC_ADD_AND_FETCH)
  CHECK_FUNCTIONWITHHEADER_EXISTS("__sync_sub_and_fetch((int*)0,0)" "${HEADERS}" HAVE_SYNC_SUB_AND_FETCH)
  CHECK_FUNCTIONWITHHEADER_EXISTS("InterlockedIncrement((long*)0)" "${HEADERS}" HAVE_INTERLOCKED_INCREMENT)
  CHECK_FUNCTIONWITHHEADER_EXISTS("InterlockedDecrement((long*)0)" "${HEADERS}" HAVE_INTERLOCKED_DECREMENT)
  CHECK_FUNCTIONWITHHEADER_EXISTS(nanosleep "${HEADERS}" HAVE_PROTOTYPE_NANOSLEEP)

# Check for thread type
IF (HAVE_WINDOWS_H)
    SET(HAVE_INT_TYPE_PTHREAD_T 1)
ELSE (HAVE_WINDOWS_H)
    DCMTK_TRY_COMPILE(HAVE_INT_TYPE_PTHREAD_T "pthread_t is an integer type"
        "// test to see if pthread_t is a pointer type or not

#include <pthread.h>

int main ()
{
  pthread_t p;
  unsigned long l = p;
  return 0;
}")
    IF (NOT HAVE_INT_TYPE_PTHREAD_T)
        SET(HAVE_POINTER_TYPE_PTHREAD_T 1 CACHE INTERNAL "Set if pthread_t is a pointer type")
    ELSE (NOT HAVE_INT_TYPE_PTHREAD_T)
        SET(HAVE_POINTER_TYPE_PTHREAD_T 0 CACHE INTERNAL "Set if pthread_t is a pointer type")
    ENDIF (NOT HAVE_INT_TYPE_PTHREAD_T)
ENDIF(HAVE_WINDOWS_H)

# Check if typename works properly. Only MSC6 really fails here.
DCMTK_TRY_COMPILE(HAVE_TYPENAME "typename works correctly"
    "
template<typename type>
type identical(type arg)
{
    return arg;
}

template<typename type>
typename type::a foo(typename type::a arg)
{
    // This one is the usage of \"typename\" that we are checking for.
    // MSC6 refuses this one incorrectly.
    return identical<typename type::a>(arg);
}

class test {
public:
    typedef int a;
};

int main()
{
    return foo<test>(0);
}")

# Check if ENAMETOOLONG is defined.
DCMTK_TRY_COMPILE(HAVE_ENAMETOOLONG "ENAMETOOLONG is defined"
    "#include <errno.h>

int main()
{
    int value = ENAMETOOLONG;
    return 0;
}")

# Check if strerror_r returns a char* is defined.
DCMTK_TRY_COMPILE(HAVE_INT_STRERROR_R "strerror_r returns an int"
    "#include <string.h>

int main()
{
    char *buf = 0;
    int i = strerror_r(0, buf, 100);
    return i;
}")
IF(HAVE_INT_STRERROR_R)
  SET(HAVE_CHARP_STRERROR_R 0 CACHE INTERNAL "Set if strerror_r() returns a char*")
ELSE(HAVE_INT_STRERROR_R)
  SET(HAVE_CHARP_STRERROR_R 1 CACHE INTERNAL "Set if strerror_r() returns a char*")
ENDIF(HAVE_INT_STRERROR_R)

# Check if variable length arrays are supported.
DCMTK_TRY_COMPILE(HAVE_VLA "variable length arrays are supported"
    "int main()
{
    int n = 42;
    int foo[n];
    return 0;
}")

# Check if std::ios::nocreate exists
DCMTK_TRY_COMPILE(HAVE_IOS_NOCREATE "std::ios::nocreate exists"
    "#include <fstream>

using namespace std;

int main()
{
  std::ifstream file(\"name\", std::ios::nocreate);
  return 0;
}")

# Check for explicit large file support
DCMTK_TRY_COMPILE(_LARGEFILE64_SOURCE "explicit large file support is available"
    "#define _LARGEFILE64_SOURCE
#include <stdio.h>

using namespace std;

int main()
{
  FILE *f = fopen64(\"name\", \"r\");
  return 0;
}")

IF(WIN32)
  # If someone can tell me how to convince TRY_COMPILE to link against winsock,
  # we could use tests for these. Until then, here is what would be the result:
  SET(HAVE_INTP_ACCEPT 1 CACHE INTERNAL "Set if socket functions accept an int* argument")
  SET(HAVE_INTP_GETSOCKOPT 1 CACHE INTERNAL "Set if socket functions accept an int* argument")
  SET(HAVE_INTP_SELECT 0 CACHE INTERNAL "Set if select() accepts an int* argument")
ELSE(WIN32)
  # Check if socket functions accept an int*
  DCMTK_TRY_COMPILE(HAVE_INTP_SOCKET, "socket functions accept an int* argument"
      "
#ifdef __cplusplus
extern \"C\" {
#endif
#ifdef _WIN32
/* Windows is pure evil */
#include <windows.h>
#else
#include <sys/socket.h>
#endif
#ifdef __cplusplus
}
#endif

int main()
{
    int i;
    struct sockaddr *addr = 0;
    int addrlen = 0;
    int optlen = 0;

    i = accept(1, addr, &addrlen);
    i = getsockopt(0, 0, 0, 0, &optlen);

    return 0;
}")
  IF(HAVE_INTP_SOCKET)
    SET(HAVE_INTP_ACCEPT 1 CACHE INTERNAL "Set if socket functions accept an int* argument")
    SET(HAVE_INTP_GETSOCKOPT 1 CACHE INTERNAL "Set if socket functions accept an int* argument")
  ELSE(HAVE_INTP_SOCKET)
    SET(HAVE_INTP_ACCEPT 0 CACHE INTERNAL "Set if socket functions accept an int* argument")
    SET(HAVE_INTP_GETSOCKOPT 0 CACHE INTERNAL "Set if socket functions accept an int* argument")
  ENDIF(HAVE_INTP_SOCKET)

  # Check if select() accepts an int*
  DCMTK_TRY_COMPILE(HAVE_INTP_SELECT "select() accepts an int* argument"
      "
#ifdef __cplusplus
extern \"C\" {
#endif
#ifdef _WIN32
/* Windows is pure evil */
#include <windows.h>
#else
#include <sys/socket.h>
#endif
#ifdef __cplusplus
}
#endif

int main()
{
    int i;
    int fds = 0;

    i = select(1, &fds, &fds, &fds, 0);

    return 0;
}
")
ENDIF(WIN32)

# Windows being windows, it lies about its processor type to 32 bit binaries
SET(SYSTEM_PROCESSOR "$ENV{PROCESSOR_ARCHITEW6432}")
IF(NOT SYSTEM_PROCESSOR)
  SET(SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}")
ENDIF(NOT SYSTEM_PROCESSOR)
# CMake doesn't provide a configure-style system type string
SET(CANONICAL_HOST_TYPE "${SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}")
UNSET(SYSTEM_PROCESSOR)

#
# End of stuff from ITK's git repository
#

#Include(CheckLibraryExists)
#CHECK_LIBRARY_EXISTS(iostream HAVE_LIBIOSTREAM)
#CHECK_LIBRARY_EXISTS(nsl HAVE_LIBNSL)
#CHECK_LIBRARY_EXISTS(socket HAVE_LIBSOCKET)

SET(CMAKE_EXTRA_INCLUDE_FILES sys/types.h)
CHECK_TYPE_SIZE("pid_t" SYS_TYPES_PID_T)
CHECK_TYPE_SIZE("size_t" SYS_TYPES_SIZE_T)
CHECK_TYPE_SIZE("ssize_t" SYS_TYPES_SSIZE_T)
CHECK_TYPE_SIZE("uchar" SYS_TYPES_UCHAR)
CHECK_TYPE_SIZE("uint" SYS_TYPES_UINT)
CHECK_TYPE_SIZE("ulong" SYS_TYPES_ULONG)
SET(CMAKE_EXTRA_INCLUDE_FILES) 

Include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(DIR ndir.h HAVE_NDIR_H)
CHECK_SYMBOL_EXISTS(DIR sys/dir.h HAVE_SYS_DIR_H)
CHECK_SYMBOL_EXISTS(DIR sys/ndir.h HAVE_SYS_NDIR_H)
IF(HAVE_IOSTREAM_H)
  CHECK_SYMBOL_EXISTS(ios::nocreate iostream.h HAVE_IOS_NOCREATE)
ENDIF(HAVE_IOSTREAM_H)

SET(USE_STD_CXX_INCLUDES 1)

SET(HAVE_STD_NAMESPACE ${CMAKE_STD_NAMESPACE})
SET(HAVE_CLASS_TEMPLATE 1)
SET(HAVE_FUNCTION_TEMPLATE 1)
SET(HAVE_EXPLICIT_TEMPLATE_SPECIALIZATION 1)
SET(HAVE_CXX_BOOL 1)

## Copy variables from CMake's internal names to their names in configurable header 
SET(SIZEOF_CHAR ${CMAKE_SIZEOF_CHAR})
SET(SIZEOF_DOUBLE ${CMAKE_SIZEOF_DOUBLE})
SET(SIZEOF_FLOAT ${CMAKE_SIZEOF_FLOAT})
SET(SIZEOF_INT ${CMAKE_SIZEOF_INT})
SET(SIZEOF_LONG ${CMAKE_SIZEOF_LONG})
SET(SIZEOF_SHORT ${CMAKE_SIZEOF_SHORT})
SET(SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})

## Created configured configuration header
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/dcmtk/config/cfunix.h.in ${CMAKE_CURRENT_BINARY_DIR}/cfunix.h @ONLY)

INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${ZLIB_INCLUDE_DIRS})

INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(dcmtk_SRCS
  cmdlnarg.cc
  dcbytstr.cc
  dcchrstr.cc
  dccodec.cc
  dcdatset.cc
  dcdebug.cc
  dcdicdir.cc
  dcdicent.cc
  dcdict.cc
  dcdictzz.cc
  dcdirrec.cc
  dcelem.cc
  dcerror.cc
  dcfilefo.cc
  dchashdi.cc
  dcistrma.cc
  dcistrmf.cc
  dcistrmz.cc
  dcitem.cc
  dclist.cc
  dcmetinf.cc
  dcobject.cc
  dcostrma.cc
  dcostrmf.cc
  dcostrmz.cc
  dcpcache.cc
  dcpixel.cc
  dcpixseq.cc
  dcpxitem.cc
  dcsequen.cc
  dcstack.cc
  dcswap.cc
  dctag.cc
  dctagkey.cc
  dctypes.cc
  dcuid.cc
  dcvm.cc
  dcvr.cc
  dcvrae.cc
  dcvras.cc
  dcvrat.cc
  dcvrcs.cc
  dcvrda.cc
  dcvrds.cc
  dcvrdt.cc
  dcvrfd.cc
  dcvrfl.cc
  dcvris.cc
  dcvrlo.cc
  dcvrlt.cc
  dcvrobow.cc
  dcvrof.cc
  dcvrpn.cc
  dcvrpobw.cc
  dcvrsh.cc
  dcvrsl.cc
  dcvrss.cc
  dcvrst.cc
  dcvrtm.cc
  dcvrui.cc
  dcvrul.cc
  dcvrulup.cc
  dcvrus.cc
  dcvrut.cc
  dcxfer.cc
  didocu.cc
  diutils.cc
  ofcond.cc
  ofconsol.cc
  ofcrc32.cc
  ofdate.cc
  ofdatime.cc
  oflist.cc
  ofstd.cc
  ofstring.cc
  ofthread.cc
  oftime.cc
  # files added in dcmtk 3.6.0
  #  oflog:
  apndimpl.cc
  appender.cc
  config.cc
  consap.cc
  factory.cc
  fileap.cc
  filter.cc
  globinit.cc
  hierarchy.cc
  hierlock.cc
  layout.cc
  lloguser.cc
  logevent.cc
  logger.cc
  logimpl.cc
  loglevel.cc
  loglog.cc
  ndc.cc
  ntelogap.cc
  nullap.cc
  objreg.cc
  oflog.cc
  patlay.cc
  pointer.cc
  property.cc
  rootlog.cc
  sleep.cc
  sockbuff.cc
  socketap.cc
  socket.cc
  strhelp.cc
  syncprims.cc
  syslogap.cc
  threads.cc
  timehelp.cc
  # ofstd:
  ofcmdln.cc
  ofconapp.cc
  #  dcmdata:
  dcwcache.cc
  vrscan.cc
  vrscanl.c
)

IF(WIN32)   # special handling for Windows systems

  IF(MINGW)
    # Avoid auto-importing warnings on MinGW
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-auto-import")
  ELSE(MINGW)
    # On Visual Studio 8 MS deprecated C. This removes all 1.276E1265 security warnings.
    IF(NOT BORLAND)
      IF(NOT CYGWIN)
        IF(NOT DCMTK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
          ADD_DEFINITIONS(
            -D_CRT_FAR_MAPPINGS_NO_DEPRECATE
            -D_CRT_IS_WCTYPE_NO_DEPRECATE
            -D_CRT_MANAGED_FP_NO_DEPRECATE
            -D_CRT_NONSTDC_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE_GLOBALS
            -D_CRT_SETERRORMODE_BEEP_SLEEP_NO_DEPRECATE
            -D_CRT_TIME_FUNCTIONS_NO_DEPRECATE
            -D_CRT_VCCLRIT_NO_DEPRECATE
            -D_SCL_SECURE_NO_DEPRECATE
            )
        ENDIF(NOT DCMTK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
      ENDIF(NOT CYGWIN)
    ENDIF(NOT BORLAND)
  ENDIF(MINGW)

ELSE(WIN32)   # ... for non-Windows systems

  # Compiler flags for Mac OS X
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_XOPEN_SOURCE_EXTENDED -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_DARWIN_C_SOURCE")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_XOPEN_SOURCE_EXTENDED -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_DARWIN_C_SOURCE")
  ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # FreeBSD fails with these flags
    IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
      # Compiler flags for all other non-Windows systems
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_XOPEN_SOURCE_EXTENDED -D_XOPEN_SOURCE=500 -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_POSIX_C_SOURCE=199506L")
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_XOPEN_SOURCE_EXTENDED -D_XOPEN_SOURCE=500 -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_POSIX_C_SOURCE=199506L")
    ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

ENDIF(WIN32)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  LIST(APPEND dcmtk_SRCS   winconap.cc windebap.cc winsock.cc)
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  LIST(APPEND dcmtk_SRCS unixsock.cc)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

ADD_LIBRARY(D4CMTK ${dcmtk_SRCS})
TARGET_LINK_LIBRARIES(D4CMTK ${ZLIB_LIBRARIES})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  TARGET_LINK_LIBRARIES(D4CMTK wsock32 netapi32)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  
INSTALL(TARGETS D4CMTK
  RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT runtime
  LIBRARY DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT runtime
  ARCHIVE DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/dcmtk/*/*.h")
INSTALL(FILES ${files} DESTINATION ${CMTK_INSTALL_INCLUDE_DIR}/dcmtk COMPONENT headers)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/dicom.dic DESTINATION ${CMTK_INSTALL_DATA_DIR} COMPONENT runtime)
