/*
    Numdiff - compare putatively similar files, 
    ignoring small numeric differences
    Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010  Ivano Primi  <ivprimi@libero.it>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _NDSELECT_H_
#define _NDSELECT_H_

#include"config.h"

/* Error codes */
#define OPEN_ERROR     -1
#define READ_ERROR     -2
#define CLOSE_ERROR    -3

typedef struct {
  /* Mask of the options */
  unsigned long optmask;

  /* Begin, end, step */
  unsigned long begin_line, end_line, step;

  /* File to handle */
  const char *file;
} Argslist ; /* A structure of this type is used to store the options */
/* set by the user                                                     */ 

#define ___H_MASK   0x00000001 /* -h option, used to recall help */
#define ___B_MASK   0x00000002 /* -b option, used to set the start line */
#define ___E_MASK   0x00000004 /* -e option, used to set the end line */
#define ___S_MASK   0x00000008 /* -s option, used to explicitly set the step */
#define ___V_MASK   0x00000010 /* -v option, used to show version number,
				  Copyright and No-Warranty */

/* I18N and L10N support */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) gettext (String)
#define gettext_noop(String) String
#define N_(String) gettext_noop (String)
#else
#define _(String) (String)
#define N_(String) String
#define textdomain(Domain)
#define bindtextdomain(Package, Directory)
#endif

#ifndef PACKAGE2
#define PACKAGE2 "ndselect"
#endif

#ifndef LOCALEDIR
#define LOCALEDIR "/usr/local/share/locale/"
#endif

#endif /* _NDSELECT_H_ */
