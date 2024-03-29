/*
 * "$Id: config.h.in 387 2009-04-18 17:05:52Z mike $"
 *
 * Configuration file for Mini-XML, a small XML-like file parsing library.
 *
 * Copyright 2003-2009 by Michael Sweet.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

/*
 * Include necessary headers...
 */

/*
 * add "cmtk_" prefix to exported symbols - we are building this library bundled with CMTK
 */
#include <cmtk_mxml_mangle.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>


/*
 * Version number...
 */

#define MXML_VERSION "@MXML_VERSION@"


/*
 * Inline function support...
 */

#define inline @C_INLINE@


/*
 * Long long support...
 */

#cmakedefine HAVE_LONG_LONG


/*
 * Do we have the snprintf() and vsnprintf() functions?
 */

#cmakedefine HAVE_SNPRINTF
#cmakedefine HAVE_VSNPRINTF


/*
 * Do we have the strXXX() functions?
 */

#cmakedefine HAVE_STRDUP


/*
 * Do we have threading support?
 */

#cmakedefine HAVE_PTHREAD_H


/*
 * Define prototypes for string functions as needed...
 */

#  ifndef HAVE_STRDUP
extern char	*_mxml_strdup(const char *);
#    define strdup _mxml_strdup
#  endif /* !HAVE_STRDUP */

extern char	*_mxml_strdupf(const char *, ...);
extern char	*_mxml_vstrdupf(const char *, va_list);

#  ifndef HAVE_SNPRINTF
extern int	_mxml_snprintf(char *, size_t, const char *, ...);
#    define snprintf _mxml_snprintf
#  endif /* !HAVE_SNPRINTF */

#  ifndef HAVE_VSNPRINTF
extern int	_mxml_vsnprintf(char *, size_t, const char *, va_list);
#    define vsnprintf _mxml_vsnprintf
#  endif /* !HAVE_VSNPRINTF */


#ifdef _MSC_VER
#if _MSC_VER >= 1400
/* disable warnings about "deprecated" C runtime functions  */
#pragma warning( disable : 4996 )
#endif
#endif

/*
 * End of "$Id: config.h.in 387 2009-04-18 17:05:52Z mike $".
 */
