##
##  Copyright 2012 SRI International
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

##
##
## THIS IS BASED ON THE KDE "FindGettext.cmake" FILE
##
##

# Try to find Gettext functionality
# Once done this will define
#
#  GETTEXT_FOUND - system has Gettext
#  GETTEXT_INCLUDE_DIR - Gettext include directory
#  GETTEXT_LIBRARIES - Libraries needed to use Gettext

# TODO: This will enable translations only if Gettext functionality is
# present in libc. Must have more robust system for release, where Gettext
# functionality can also reside in standalone Gettext library, or the one
# embedded within kdelibs (cf. gettext.m4 from Gettext source).
#
# Copyright (c) 2006, Chusslove Illich, <caslav.ilic@gmx.net>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


if (LIBC_HAS_GETTEXT OR LIBINTL_HAS_GETTEXT)

  # in cache already
  set(GETTEXT_FOUND TRUE)

else (LIBC_HAS_GETTEXT OR LIBINTL_HAS_GETTEXT)

  include(CheckLibraryExists)
  include(CheckFunctionExists)

  CHECK_INCLUDE_FILES(libintl.h HAVE_LIBINTL_H)

  set(GETTEXT_LIBRARIES)

  if (HAVE_LIBINTL_H)
    check_function_exists(_libintl_gettext LIBC_HAS_GETTEXT)
    if (LIBC_HAS_GETTEXT)
      set(GETTEXT_SOURCE "built in libc")
      set(GETTEXT_LIBRARIES "")
      set(GETTEXT_FOUND TRUE)
    else (LIBC_HAS_GETTEXT)
      find_library(LIBINTL_LIBRARY NAMES intl libintl )
      if(LIBINTL_LIBARY)
	check_library_exists(${LIBINTL_LIBRARY} "_libintl_gettext" "" LIBINTL_HAS_GETTEXT)
	if (LIBINTL_HAS_GETTEXT)
          set(GETTEXT_SOURCE "in ${LIBINTL_LIBRARY}")
          set(GETTEXT_LIBRARIES ${LIBINTL_LIBRARY})
          set(GETTEXT_FOUND TRUE)
	endif (LIBINTL_HAS_GETTEXT)
      endif(LIBINTL_LIBARY)
    endif (LIBC_HAS_GETTEXT)
  endif (HAVE_LIBINTL_H)
  
  if (GETTEXT_FOUND)
    if (NOT Gettext_FIND_QUIETLY)
      message(STATUS "Found Gettext: ${GETTEXT_SOURCE}")
    endif (NOT Gettext_FIND_QUIETLY)
  else (GETTEXT_FOUND)
    if (Gettext_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find Gettext")
    endif (Gettext_FIND_REQUIRED)
  endif (GETTEXT_FOUND)
  
  mark_as_advanced(GETTEXT_INCLUDE_DIR GETTEXT_LIBRARIES)

endif (LIBC_HAS_GETTEXT OR LIBINTL_HAS_GETTEXT)

