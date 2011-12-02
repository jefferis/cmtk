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

#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<errno.h>
#ifdef ENABLE_NLS
#include<locale.h>
#endif
#include"getopt.h"
#include"numdiff.h" /* Only for LINE_INTERR definition */
#include"ndselect.h"

#include"read_line.c"

static int scan_file (const Argslist* data)
{
  FILE *fp;
  char *line_buffer;
  unsigned long lineno;
  int errcode;

  if (!data->file || !*data->file)
    fp = stdin;
  else
    {
      if ( !(fp = fopen (data->file, "r")) )
	return OPEN_ERROR;
    }

  if (!data->end_line)
    {
      for (lineno = 1; 
	   (line_buffer = read_line(fp, &errcode), errcode <= LINE_INTERR); 
	   lineno++)
	{
	  if (lineno >= data->begin_line && (lineno - data->begin_line) % data->step == 0)
	    fputs (line_buffer, stdout);
	  free ((void*)line_buffer);
	}
    }
  else
    {
      for (lineno = 1; 
	   lineno <= data->end_line && 
	     (line_buffer = read_line(fp, &errcode), errcode <= LINE_INTERR); 
	   lineno++)
	{
	  if (lineno >= data->begin_line && (lineno - data->begin_line) % data->step == 0)
	    fputs (line_buffer, stdout);
	  free ((void*)line_buffer);
	}
    }

  if (errcode <= EOF_REACHED)
    {
      if (fp == stdin)
	return OK;
      else
	return fclose (fp) == EOF ? CLOSE_ERROR : OK;
    }
  else
    {
      if((line_buffer))
	free ((void*)line_buffer);
      if (fp != stdin)
	fclose (fp);
      return READ_ERROR;
    }
}

static void print_selversion (const char* progname)
{
  printf ("%s %s\n", progname, VERSION);
  printf ("Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010  %s <ivprimi@libero.it>\n", 
	  /* TRANSLATORS: This is a proper name.  See the gettext
	     manual, section Names.
	     Pronounciation is like "evaa-no pree-me".  */
	  _("Ivano Primi"));
  printf (_("\
License GPLv3+: GNU GPL version 3 or later,\n\
see <http://gnu.org/licenses/gpl.html>.\n\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.\n"));
}

static void print_selhelp (const char* progname)
{
  puts (_("Usage:"));
  printf ("%s -h|--help|-v|--version   %s\n\n", progname, _("or"));
  printf ("%s %s\n", progname, "[-b FIRST][-e LAST][-s STEP][-l PATH][-o PATH] [FILE]");
  /* %%% */
  printf(_("\n\
For a given range of line numbers and a given step\n\
print on the standard output all lines of a file,\n\
starting with the first line of the range and ending within\n\
its last line, whose line number is such that the difference\n\
between it and the start point is an integer multiple\n\
of the given step\n\n"));
  printf ("-b, --beginning, --start=FIRST\n    %s\n    %s\n", 
	  _("Specify the first line to print"),
	  _("(The default behavior is to start with line number 1)"));
  printf ("-e, --end=LAST\n    %s\n    %s\n", 
	  _("Specify the last line that can be printed"),
	  _("(The default behavior is to arrive till to the end of the file)"));
  printf ("-s, --step=STEP\n    %s\n    %s\n", 
	  _("Specify the step to use when selecting the lines to print"),
	  _("(The default value for the step is 1)"));
  printf ("-l, --warnings-to=PATH\n    %s\n",
	  _("Redirect warning and error messages from stderr to the indicated file"));
  printf ("-o, --output=PATH\n    %s\n",
	  _("Redirect output from stdout to the indicated file"));
  printf ("-h, --help\n    %s\n", _("Show this help message"));
  printf ("-v, --version\n    %s\n", _("Show version number, Copyright, Distribution Terms and NO-Warranty"));
  printf ("\n%s\n%s\n%s\n%s\n\n",
	  _("The argument after the options is the name of the file to scan."),
	  _("The complete path of the file should be given,\na directory name is not accepted."),
	  _("If no input file is specified, the program reads from the standard input."),
	  _("Exit status: 0 in case of normal termination, -1 (255) in case of error"));
}

extern int errno;
extern char *optarg; 
extern int optind;

static int set_args (int argc, char* argv[], Argslist *list)
{
  const char *optstring = "hb:e:s:l:o:v";
  struct option long_options[] = {
    {"help",         0, NULL, 'h'},
    {"beginning",    1, NULL, 'b'},
    {"start",        1, NULL, 'b'},
    {"end",          1, NULL, 'e'},
    {"step",         1, NULL, 's'},
    {"warnings-to",  1, NULL, 'l'},
    {"output",       1, NULL, 'o'},
    {"version",      0, NULL, 'v'},
    {0, 0, 0, 0}
  };
  int option_index=0;
  char *endptr;
  int optch; 

  /*
    We start by loading the default values
    for the user settable options.
  */
  list->optmask = 0x0;
  list->begin_line=1;
  list->end_line=0;
  list->step=1;
  list->file=NULL;

  while ( (optch = getopt_long (argc, argv, optstring, long_options, &option_index)) != -1 )
    {
      switch (optch)
	{
	case 'h':
	  list->optmask |= ___H_MASK;
	  break;
	case 'b':
	  list->optmask |= ___B_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->begin_line = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) || list->begin_line == 0 ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 'e':
	  list->optmask |= ___E_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->end_line = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 's':
	  list->optmask |= ___S_MASK;
	  errno = 0;
	  for (endptr = optarg; is_space (*endptr) != 0; endptr++);
	  if (*endptr == '-' ||
	      (list->step = strtoul (optarg, &endptr, 10), 
	       errno == ERANGE) || list->step == 0 ||
	      *endptr != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE2, optch);
	      return -1;
	    }
	  break;
	case 'l':
	  if (!freopen (optarg, "w", stderr))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE2, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'o':
	  if (!freopen (optarg, "w", stdout))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE2, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'v':
	  list->optmask |= ___V_MASK;
	  break;
	default:
/*  	  fprintf (stderr,  */
/*  		   _("%s: unrecognized option `-%c\' \n"), PACKAGE2, optch);  */
	  return -1;
	}
    }
  if (!(list->optmask & (___H_MASK | ___V_MASK)) && argc - optind > 1)
    {
      print_selhelp (PACKAGE2);
      return -1;
    }
  else
    {
      if( !(list->optmask & (___H_MASK | ___V_MASK)) )
	list->file = (const char*) argv[optind];
      return 0;
    }
}

int main (int argc, char** argv)
{
  Argslist arg_list;

#ifdef ENABLE_NLS
  setlocale (LC_CTYPE, "");
  setlocale (LC_MESSAGES, "");
#endif
  bindtextdomain (PACKAGE2, LOCALEDIR);
  textdomain (PACKAGE2);
  if ( set_args (argc, argv, &arg_list) != 0 )
    return -1;
  else if ( (arg_list.optmask & (___H_MASK | ___V_MASK)) )
    {
      if ((arg_list.optmask & ___V_MASK))
      	print_selversion(PACKAGE2);
      if ((arg_list.optmask & ___H_MASK))
      	print_selhelp(PACKAGE2);
      if (argc > 2)
	return -1;
      else
	return 0;
    }
  else
    {
      switch (scan_file (&arg_list))
	{
	case OPEN_ERROR:
	  fprintf (stderr, _("%s: cannot open file \"%s\":\n"), PACKAGE2, arg_list.file);
	  perror(0);
	  return -1;
	case CLOSE_ERROR:
	  fprintf (stderr, _("%s: cannot close file \"%s\":\n"), PACKAGE2, arg_list.file);
	  perror(0);
	  return -1;
	case READ_ERROR:
	  fprintf (stderr, _("%s: Error occurred while reading from file \"%s\"\n\n"), PACKAGE2, arg_list.file);
	  return -1;
	default:
	  return 0;
	}
    }
}
