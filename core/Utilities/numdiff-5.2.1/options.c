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

#define GDIFF_OPTIONS 1

/* Leave this inclusion at the begin, otherwise problems */
/* with the symbol __USE_FILE_OFFSET64                   */
#include"numdiff.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<float.h>
#include"getopt.h"
#include"error.h"
#include"xalloc.h"

#ifdef _DMALLOC_
#include <dmalloc.h> /* Useful only for the debugging */
#endif

void print_version (const char* progname)
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
#ifdef USE_GMP
  printf ("\n%s %s.\n", _("The software has been linked against\n\
the GNU Multiple Precision Arithmetic Library,\n\
version number"), gmp_version);
#else  /* not USE_GMP */
  printf ("\n%s.\n", _("The software has been built with\n\
its own internal support for multiple precision arithmetic"));
#endif /* not USE_GMP */
}

void print_help (const char* progname)
{
  puts (_("Usage:"));
  printf ("%s -h|--help|-v|--version   %s\n\n", progname, _("or"));
  printf ("%s %s\n", progname, "[-s IFS][-a MAXERR][-r MAXERR][-2][-# PREC][-P][-N][-I][-d C1C2][-t C1C2][-g N1N2][-p C1C2][-n C1C2][-e C1C2][-i C1C2][-X 1:INT][-X 2:INT][-E][-D][-b][-V][-q][-S][-z 1:INT][-z 2:INT][-Z 1:INT][-Z 2:INT][-m][-H][-f[NUM]][-T][-l PATH][-o PATH] FILE1 FILE2");
  printf (_("\nINT stays for a positive integer value or for a range of integer values,\nlike 1-, 3-5 or -7\n"));
  /* %%% */
  printf (_("\nCompare putatively similar files line by line and field by field,\nignoring small numeric differences or/and different numeric formats\n\n"));
  printf ("-s, --separator=IFS\n    %s\n    %s\n",
	  _("Specify the set of characters to use\n    to split the input lines into fields"),
	  _("(The default set of characters is space, tab and newline)"));
  printf ("-a, --absolute-tolerance=MAXERR\n    %s\n    %s\n", 
	  _("Specify the maximum absolute difference permitted\n    before that two numeric fields are regarded as different"),
	  _("(The default value is zero)"));
  printf ("-r, --relative-tolerance=MAXERR\n    %s\n    %s\n", 
	  _("Specify the maximum relative difference permitted\n    before that two numeric fields are regarded as different"),
	  _("(The default value is zero)"));
  printf ("-2, --strict\n    %s\n",
	  _("Order that two numerical values are regarded as equal only if\n    both absolute and relative difference do not exceed\n    the corresponding tolerance threshold"));
  printf ("-#, --digits=PREC\n    %s\n",
	  _("Specify the number of digits in the significands\n    used in multiple precision arithmetic"));
  printf ("-P, --positive-differences\n    %s\n",
	  _("Ignore all differences due to numeric fields of the second file that\n    are less than the corresponding numeric fields in the first file"));
  printf ("-N, --negative-differences\n    %s\n",
	  _("Ignore all differences due to numeric fields of the second file that\n    are greater than the corresponding numeric fields in the first file"));
  printf ("-I, --ignore-case\n    %s\n",
	  _("Ignore changes in case while doing literal comparisons"));
  printf ("-d, --decimal-point=C1C2\n    %s\n",
	  _("Specify the characters representing the decimal point\n    in the two files to compare"));
  printf ("-t, --thousands-separator=C1C2\n    %s\n",
	  _("Specify the characters representing the thousands separator\n    in the two files to compare"));
  printf ("-g, --group-length=N1N2\n    %s\n",
	  _("Specify the number of digits forming each group of thousands\n    in the two files to compare"));
  printf ("-p, --plus-prefix=C1C2\n    %s\n",
	  _("Specify the (optional) prefixes for positive values\n    used in the two files to compare"));
  printf ("-n, --minus-prefix=C1C2\n    %s\n",
	  _("Specify the prefixes for negative values\n    used in the two files to compare"));
  printf ("-e, --exponent-letter=C1C2\n    %s\n",
	  _("Specify the exponent letters\n    used in the two files to compare"));
  printf ("-i, --imaginary-unit=C1C2\n    %s\n",
	  _("Specify the characters representing the imaginary unit\n    in the two files to compare"));
  printf ("-X, --exclude=1:INT\n    %s\n",
	  _("Select the fields of the first file that have to be ignored"));
  printf ("-X, --exclude=2:INT\n    %s\n",
	  _("Select the fields of the second file that have to be ignored"));
  printf ("-E, --essential\n    %s\n",
	  _("While printing the differences between the two compared files\n    show only the numerical ones"));
  printf ("-D, --dummy\n    %s\n",
	  _("While printing the differences between the two compared files\n    neglect all the numerical ones (dummy mode)"));
  printf ("-b, --brief\n    %s\n",
	  _("Suppress all messages concerning the differences discovered\n    in the structures of the two files"));
  printf ("-V, --verbose\n    %s\n",
	  _("For every couple of lines which differ in at least one field print\n    an header to show how these lines appear in the two compared files"));
  printf ("-q, --quiet, --silent\n    %s\n",
	  _("Suppress all the standard output"));
  printf ("-S, --statistics\n    %s\n",
	  _("Add some statistics to the standard output"));
  printf ("-z, --blur-if-numerical=1:INT\n    %s\n",
	  _("Select the fields of the first file that have to be\n    blurred during the synchronization procedure\n    only if they turn out to be numeric"));
  printf ("-z, --blur-if-numerical=2:INT\n    %s\n",
	  _("Select the fields of the second file that have to be\n    blurred during the synchronization procedure\n    only if they turn out to be numeric"));
  printf ("-Z, --blur-unconditionally=1:INT\n    %s\n",
	  _("Select the fields of the first file that have to be\n    unconditionally blurred during the synchronization procedure"));
  printf ("-Z, --blur-unconditionally=2:INT\n    %s\n",
	  _("Select the fields of the second file that have to be\n    unconditionally blurred during the synchronization procedure"));
  printf ("-m, --minimal\n    %s\n",
	  _("During synchronization try hard to find a smaller set of changes"));
  printf ("-H, --speed-large-files\n    %s\n",
	  _("During synchronization assume large files and\n    many scattered small changes"));
  printf ("-f, --test-filter[=NUM]\n    %s\n    %s\n    %s\n    %s\n",
	  _("Run only the filter and then show the results of its\n    attempt to synchronize the two files."),
	  _("If \'NUM\' is zero or is not specified, output at most 130 columns per line."),
	  _("If \'NUM\' is a positive number, output at most\n    \'NUM\' columns per line."),
	  _("If \'NUM\' is a negative number, do not output common lines\n    and display at most -\'NUM\' columns per line."));
  printf ("-T, --expand-tabs\n    %s\n",
	  _("Expand tabs to spaces in output while displaying the results of the\n    synchronization procedure (meaningful only together with option -f)"));
  printf ("-l, --warnings-to=PATH\n    %s\n",
	  _("Redirect warning and error messages from stderr to the indicated file"));
  printf ("-o, --output=PATH\n    %s\n",
	  _("Redirect output from stdout to the indicated file"));
  printf ("-h, --help\n    %s\n", _("Show help message and predefined settings"));
  printf ("-v, --version\n    %s\n", _("Show version number, Copyright, Distribution Terms and NO-Warranty"));
  printf ("\n%s\n%s\n%s\n%s\n",
	  _("The two arguments after the options are the names of the files to compare."),
	  _("The complete paths of the files should be given,\na directory name is not accepted."),
	  _("They cannot refer to the same file but one of them can be \"-\",\nwhich refers to stdin."),
	  _("Exit status: 1 if files differ, 0 if they are equal, -1 (255) in case of error"));
  /* %%% */
  puts (_("\n  Default numeric format (for both files to compare):\n"));
  printf (_("Decimal point = `%c\'\n"), DP);

  printf (_("Thousands separator = `%c\'\n"), THSEP);
  printf (_("Number of digits in each thousands group = %u\n"), GROUPING);

  printf (_("Leading positive sign = `%c\'\n"), POS_SIGN);
  printf (_("Leading negative sign = `%c\'\n"), NEG_SIGN);
  printf (_("Prefix for decimal exponent = `%c\'\n"), ECH);
  printf (_("Symbol used to denote the imaginary unit = `%c\'\n\n"), IU);
}

static
char* return_ifs (const char* optarg)
{
  char *s = (char*) malloc ((strlen(optarg) + 1) * sizeof(char));
  
  if(!s)
    return NULL;
  else
    {
      char *t, *u;

      strcpy (s, optarg);
      for (t = s; *t != '\0'; t++)
	{
	  if (*t == '\\')
	    {
	      switch (*(t+1))
		{
		case 'f':
		  *t = '\f';
		  break;
		case 'n':
		  *t = '\n';
		  break;
		case 'r':
		  *t = '\r';
		  break;
		case 't':
		  *t = '\t';
		  break;
		case 'v':
		  *t = '\v';
		  break;
		default:
		  *t = *(t+1);
		}
	      for (u = t+1; *u != '\0'; *u = *(u+1), u++);
	    }
	}
      return s;
    }
}

static
int nfset (int opt_ch, const char* opt_arg, argslist* arg_list)
{
  if (strlen(opt_arg) <= 2)
    {
      char _1st = *opt_arg, _2nd = *(opt_arg+1);

      switch (opt_ch)
	{
	case 'd':
	  if ( (is_punct(_1st)) && (_2nd == '\0' || is_punct(_2nd)) )
	    {
	      arg_list->optmask |= _D_MASK;
	      arg_list->nf1.dp = _1st;
	      arg_list->nf2.dp = (_2nd) ? _2nd : _1st; 
	      return 0;
	    }
	  break;
	case 't':
	  if ( (is_punct(_1st)) && (_2nd == '\0' || is_punct(_2nd)) )
	    {
	      arg_list->optmask |= _T_MASK;
	      arg_list->nf1.thsep = _1st;
	      arg_list->nf2.thsep = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'e':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _E_MASK;
	      arg_list->nf1.ech = _1st;
	      arg_list->nf2.ech = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'n':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _N_MASK;
	      arg_list->nf1.neg_sign = _1st;
	      arg_list->nf2.neg_sign = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'i':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _I_MASK;
	      arg_list->nf1.iu = _1st;
	      arg_list->nf2.iu = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'p':
	  if ( is_print(_1st) && (_2nd == '\0' || is_print(_2nd)) )
	    {
	      arg_list->optmask |= _P_MASK;
	      arg_list->nf1.pos_sign = _1st;
	      arg_list->nf2.pos_sign = (_2nd) ? _2nd : _1st;
	      return 0;
	    }
	  break;
	case 'g':
	  if ( (is_digit(_1st)) && (_2nd == '\0' || is_digit(_2nd)) )
	    {
	      arg_list->optmask |= _G_MASK;
	      arg_list->nf1.grouping = _1st - '0';
	      arg_list->nf2.grouping = (_2nd) ? _2nd - '0': _1st - '0';
	      return 0;
	    }
	  break;
	}
    }
  return -1;
}

static
int fselect (const char* str, unsigned char* mask, int mask_size)
{
  long beg, end;
  unsigned long n;
  char *ptr, *endptr;

  beg = end = -1;
  if (!str || !*str)
    return 0; /* no field selected */
  /* If we arrive here we are sure that *str != '\0' ! */
  if ( strcmp (str, "@") == 0 )
    {
      /* select all fields */
      for (mask_size /= 8; mask_size > 0; mask_size--, mask[mask_size] = 0xFF);       
      return 1;
    }
  if ((beg = strtol (str, &endptr, 10)) == 0
      || beg > mask_size || beg < -mask_size)
    return -1; /* illegal input */
  else if (beg < 0)
    {
      if (*endptr == '\0')
	{
	  end = -beg;
	  beg = 1;
	}
      else
	return -1;
    }
  else if (*endptr == '\0')
    end = beg;
  else if (*endptr == '-')
    {
      if (*(ptr = endptr + 1) == '\0')
	end = mask_size; 
      else
	{
	  if ((end = strtol (ptr, &endptr, 10)) <= 0
	      || *endptr != '\0' || end > mask_size)
	    return -1; /* illegal input */
	}
    }
  if (beg > end)
    return -1;
  else
    {
      for (n = beg - 1; n <= end - 1; n++)
	mask[n >> 3] |= 0x80 >> (n & 0x7);
      return 1;
    }
}

static
int valid_numfmt (const struct numfmt* pnf)
{
  char store[NUMFMT_CHARS];
  int i, j;

  store[0] = pnf->dp;
  store[1] = pnf->thsep;
  store[2] = pnf->pos_sign;
  store[3] = pnf->neg_sign;
  store[4] = pnf->ech;
  store[5] = pnf->iu;
  for (i=0; i < NUMFMT_CHARS; i++)
    {
      for (j = i+1; j < NUMFMT_CHARS; j++)
	if (store[i] == store[j])
	  return 0;
    }
  return 1;
}

extern int optind;

int setargs (int argc, char* argv[], argslist *list)
{
  const char *optstring = "h2bVqDESIPNz:Z:mHT#:s:a:r:d:t:g:p:n:e:i:f::X:l:o:v";
  struct option long_options[] = {
    {"help",                 0, NULL, 'h'},
    {"strict",               0, NULL, '2'},
    {"brief",                0, NULL, 'b'},
    {"verbose",              0, NULL, 'V'},
    {"quiet",                0, NULL, 'q'},
    {"silent",               0, NULL, 'q'},
    {"dummy",                0, NULL, 'D'},
    {"essential",            0, NULL, 'E'},
    {"statistics",           0, NULL, 'S'},
    {"ignore-case",          0, NULL, 'I'},
    {"positive-differences", 0, NULL, 'P'},
    {"negative-differences", 0, NULL, 'N'},
    {"blur-if-numerical",    1, NULL, 'z'},
    {"blur-unconditionally", 1, NULL, 'Z'},
    {"minimal",              0, NULL, 'm'},
    {"speed-large-files",    0, NULL, 'H'},
    {"expand-tabs",          0, NULL, 'T'},
    {"digits",               1, NULL, '#'},
    {"separator",            1, NULL, 's'},
    {"absolute-tolerance",   1, NULL, 'a'},
    {"relative-tolerance",   1, NULL, 'r'},
    {"decimal-point",        1, NULL, 'd'},
    {"thousands-separator",  1, NULL, 't'},
    {"group-length",         1, NULL, 'g'},
    {"plus-prefix",          1, NULL, 'p'},
    {"minus-prefix",         1, NULL, 'n'},
    {"exponent-letter",      1, NULL, 'e'},
    {"imaginary-unit",       1, NULL, 'i'},
    {"test-filter",          2, NULL, 'f'},
    {"exclude",              1, NULL, 'X'},
    {"warnings-to",          1, NULL, 'l'},
    {"output",               1, NULL, 'o'},
    {"version",              0, NULL, 'v'},
    {0, 0, 0, 0}
  };
  int option_index=0;
  char *tail;
  int i, optch, off; 
  unsigned int t;
  long w;
  unsigned char *bitmask;  
  struct numfmt defaults;

  /*
    We start by loading the default values
    for the user settable options.

    The initialization of these variables:

    list->maxrelerr, list->maxabserr,
    list->Labserr, list->Crelerr, list->Lrelerr, list->Cabserr,
    list->N1abserr, list->N1disperr, list->N2abserr, list->N2disperr 

    to Zero is done within main()
    through init_mpa_support() .
  */

  suppress_common_lines = 0;
  ignore_white_space = IGNORE_NO_WHITE_SPACE;
  expand_tabs = 0;
  w = DEF_ATMOST_NCOLS;
  speed_large_files = 0;
  program_name = PACKAGE;

  list->optmask = 0x0;
  list->output_mode = OUTMODE_NORMAL;
  for (i=0; i < FIELDMASK_SIZE; 
       list->ghostmask1[i] = list->ghostmask2[i] = list->tblurmask1[i] = list->tblurmask2[i] = list->pblurmask1[i] = list->pblurmask2[i] = 0x0, i++); 
  list->relerr_formula = CLASSIC_FORMULA;
  list->Nentries = list->Ndisperr = 0;
  list->flag = 0;
  list->ifs1 = list->ifs2 = NULL;
  list->iscale = ISCALE;
  list->nf1.dp = DP;
  list->nf1.thsep = THSEP;
  list->nf1.grouping = GROUPING;
  list->nf1.pos_sign = POS_SIGN;
  list->nf1.neg_sign = NEG_SIGN;
  list->nf1.ech = ECH;
  list->nf1.iu = IU;
  list->file1 = list->file2 = NULL;
  defaults = list->nf2 = list->nf1;

  /*
    defaults.dp == DP
    defaults.thsep == THSEP;
    defaults.grouping == GROUPING;
    defaults.pos_sign == POS_SIGN;
    defaults.neg_sign == NEG_SIGN;
    defaults.ech == ECH;
    defaults.iu == IU;

    Since 'defaults' is not modified in the following
    we are sure that it will always contain
    these default values.
  */

  while ( (optch = getopt_long (argc, argv, optstring, long_options, &option_index)) != -1 )
    {
      switch (optch)
	{
	case 'h':
	  list->optmask |= _H_MASK;
	  break;
	case '2':
	  list->optmask |= _2_MASK;
	  break;	  
	case 'b':
	  list->optmask |= _B_MASK;
	  break;
	case 'V':
	  list->optmask |= _SV_MASK;
	  break;
	case 'q':
	  list->optmask |= _Q_MASK;
	  break;
	case 'D':
	  list->optmask |= _SD_MASK;
	  break;
	case 'E':
	  list->optmask |= _SE_MASK;
	  break;
	case 'S':
	  list->optmask |= _SS_MASK;
	  break;
	case 'I':
	  list->optmask |= _SI_MASK;
	  break;
	case 'P':
	  list->optmask |= _SP_MASK;
	  list->flag = 1;
	  break;
	case 'N':
	  list->optmask |= _SN_MASK;
	  list->flag = -1;
	  break;
	case 'z':
	  if ( (i = strncmp (optarg, "1:", 2)) && (strncmp (optarg, "2:", 2)) )
	    {
	      /*
		None of the prefixes 1: and 2: has been used,
		then we have to select fields for both files
	      */
	      if (fselect (optarg, list->pblurmask1, FIELDMASK_SIZE*8) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		{
		  fselect (optarg, list->pblurmask2, FIELDMASK_SIZE*8);
		  list->optmask |= _Z_MASK;
		}
	    }
	  else 
	    {
	      bitmask = i == 0 ? list->pblurmask1 : list->pblurmask2;
	      if (fselect (optarg+2, bitmask, FIELDMASK_SIZE*8) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		list->optmask |= _Z_MASK;
	    }
	  break;
	case 'Z':
	  if ( (i = strncmp (optarg, "1:", 2)) && (strncmp (optarg, "2:", 2)) )
	    {
	      /*
		None of the prefixes 1: and 2: has been used,
		then we have to select fields for both files
	      */
	      if (fselect (optarg, list->tblurmask1, FIELDMASK_SIZE*8) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		{
		  fselect (optarg, list->tblurmask2, FIELDMASK_SIZE*8);
		  list->optmask |= _SZ_MASK;
		}
	    }
	  else 
	    {
	      bitmask = i == 0 ? list->tblurmask1 : list->tblurmask2;
	      if (fselect (optarg+2, bitmask, FIELDMASK_SIZE*8) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		list->optmask |= _SZ_MASK;
	    }
	  break;
	case 'm':
	  list->optmask |= _M_MASK;
	  break;
	case 'H':
	  list->optmask |= _SH_MASK;
	  speed_large_files = 1;
	  break;
	case 'T':
	  expand_tabs = 1;
	  break;
	case '#':
	  list->iscale = strtol (optarg, &tail, 10);
	  if (*tail != '\0' || list->iscale < 0 || list->iscale > MAX_ISCALE)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE, optch);
	      return -1;
	    }
	  else
	    list->optmask |= _X_MASK;	  
	  break;
	case 's':
	  if (*optarg == '1' && *(optarg+1) == ':')
	    {
	      if ((list->ifs1))
		{
		  free((void*)list->ifs1);
		  list->ifs1 = NULL;
		}
	      list->ifs1 = return_ifs (optarg+2);
	    }
	  else if (*optarg == '2' && *(optarg+1) == ':')
	    {
	      if ((list->ifs2))
		{
		  free((void*)list->ifs2);
		  list->ifs2 = NULL;
		}
	      list->ifs2 = return_ifs (optarg+2);
	    }
	  else
	    {
	      if ((list->ifs1))
		{
		  free((void*)list->ifs1);
		  list->ifs1 = NULL;
		}
	      if ((list->ifs2))
		{
		  free((void*)list->ifs2);
		  list->ifs2 = NULL;
		}
	      list->ifs1 = return_ifs (optarg);
	      list->ifs2 = return_ifs (optarg);
	    }
	  if ( ((list->ifs1) && !strchr(list->ifs1, NEWLINE)) ||
	       ((list->ifs2) && !strchr(list->ifs2, NEWLINE)) )
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option:\n"),
		       PACKAGE, optch);
	      fprintf (stderr, _("  A list of field delimiters can not be empty and\n  must always include the newline character (\'\\n\')\n"));
	      return -1;
	    }
	  else
	    list->optmask |= _S_MASK;
	  break;
	case 'a':
#ifndef USE_GMP
	  delR (&list->maxabserr); /* To avoid unpleasant memory leaks */
#endif /* not USE_GMP */
	  str2R (optarg, &tail, ISCALE, &defaults, &list->maxabserr);

	  if (*tail != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE, optch);
	      return -1;
	    }
	  else
	    list->optmask |= _A_MASK;
	  break;
	case 'r':
#ifndef USE_GMP
	  delR (&list->maxrelerr); /* To avoid unpleasant memory leaks */
#endif /* not USE_GMP */
	  if (*optarg == '1' && *(optarg+1) ==':')
	    {
	      list->relerr_formula = WR_TO_FIRST_FILE;
	      str2R (optarg+2, &tail, ISCALE, &defaults, &list->maxrelerr);
	    }
	  else if (*optarg == '2' && *(optarg+1) ==':')
	    {
	      list->relerr_formula = WR_TO_SECOND_FILE;
	      str2R (optarg+2, &tail, ISCALE, &defaults, &list->maxrelerr);
	    }
	  else
	    str2R (optarg, &tail, ISCALE, &defaults, &list->maxrelerr);

	  if (*tail != '\0')
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE, optch);
	      return -1;
	    }
	  else
	    list->optmask |= _R_MASK;
	  break;
	case 'd':
	case 't':
	case 'g':
	case 'p':
	case 'n':
	case 'e':
	case 'i':
	  if (nfset (optch, optarg, list) < 0)
	    {
	      fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
		       PACKAGE, optch);
	      return -1;
	    }
	  break;
	case 'f':
	  if(!optarg)
	    {
	      /* There is no optional argument, then set */
	      /* 'w' to 'DEF_ATMOST_NCOLS'.              */
	      list->optmask |= _F_MASK;
	      w = DEF_ATMOST_NCOLS;
	    }
	  else
	    {
	      /* An argument follows */
	      w = strtol (optarg, &tail, 10);
	      /* If the argument of the option is not a valid number, */
	      /* then exit after printing a suitable error message.   */
	      if (*tail != '\0')
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		list->optmask |= _F_MASK;
	      /* Otherwise you have to set 'w' appropriately. */
	      /* If the given argument is less than -MAX_ATMOST_NCOLS */
	      /* then set 'w' to 'DEF_ATMOST_NCOLS' and 'suppress_common_lines' to 'TRUE'. */
	      if (w < -MAX_ATMOST_NCOLS)
		{
		  w = DEF_ATMOST_NCOLS;
		  suppress_common_lines = 1;
		}
	      /* If the argument were negative, then remove the sign */
	      /* and set 'suppress_common_lines' to 'TRUE'.          */
	      if (w < 0)
		{
		  w *= -1;
		  suppress_common_lines = 1;
		}
	      /* If the given argument is too small or too big in absolute value, */
	      /* then set 'w' to 'DEF_ATMOST_NCOLS'.                              */
	      if (w < MIN_ATMOST_NCOLS || w > MAX_ATMOST_NCOLS)
		w = DEF_ATMOST_NCOLS;
	      /* Otherwise leave 'w' set to the value of the argument. */
	    } /* end optarg != 0 */
	  break;
	case 'X':
	  if ( (i = strncmp (optarg, "1:", 2)) && (strncmp (optarg, "2:", 2)) )
	    {
	      /*
		None of the prefixes 1: and 2: has been used,
		then we have to select fields for both files
	      */
	      if (fselect (optarg, list->ghostmask1, FIELDMASK_SIZE*8) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		{
		  fselect (optarg, list->ghostmask2, FIELDMASK_SIZE*8);
		  list->optmask |= _SX_MASK;
		}
	    }
	  else 
	    {
	      bitmask = i == 0 ? list->ghostmask1 : list->ghostmask2;
	      if (fselect (optarg+2, bitmask, FIELDMASK_SIZE*8) <= 0)
		{
		  fprintf (stderr, _("%s: invalid argument after `-%c\' option\n"),
			   PACKAGE, optch);
		  return -1;
		}
	      else
		list->optmask |= _SX_MASK; 
	    }
	  break;
	case 'l':
	  if (!freopen (optarg, "w", stderr))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'o':
	  if (!freopen (optarg, "w", stdout))
	    {
	      fprintf (stderr, _("%s: cannot open file \"%s\":\n"),
		       PACKAGE, optarg);
	      perror(0);
	      return -1;
	    }
	  break;
	case 'v':
	  list->optmask |= _V_MASK;
	  break;
	default:
	  /* 	  
		  fprintf (stderr, 
		  _("%s: unrecognized option `-%c\' \n"), PACKAGE, optch); 
	  */
	  return -1;
	}
    }

  t = expand_tabs ? 1 : TAB_WIDTH;
  off = (w + t + 3) / (2 * t)  *  t;
  sdiff_half_width = MAX (0, MIN (off - 3, w - off)),
    sdiff_column2_offset = sdiff_half_width ? off : w;

  if ( list->optmask & _SV_MASK )
    list->output_mode = OUTMODE_VERBOSE;
  if ( list->optmask & _B_MASK )
    list->output_mode = OUTMODE_BRIEF;
  if (list->optmask & _B_MASK && list->optmask & _SV_MASK)
    list->output_mode = OUTMODE_COINCISE;
  if ( list->optmask & _Q_MASK )
    list->output_mode = OUTMODE_QUIET;

  if (!(list->optmask & (_H_MASK | _V_MASK)) && argc - optind != 2)
    {
      print_help (PACKAGE);
      return -1;
    }
  else if ( !valid_numfmt(&list->nf1) )
    {
      fprintf (stderr, 
	       _("The numeric format specified for the first file is illegal,\n"));
      fprintf (stderr,
	       _("the following symbols should be all different\nwhile two or more of them are actually equal:\n"));
      fprintf (stderr, _("\nDecimal point = `%c\'\n"), list->nf1.dp);
      fprintf (stderr, _("Thousands separator = `%c\'\n"), list->nf1.thsep);
      fprintf (stderr, _("Leading positive sign = `%c\'\n"), list->nf1.pos_sign);
      fprintf (stderr, _("Leading negative sign = `%c\'\n"), list->nf1.neg_sign);
      fprintf (stderr, _("Prefix for decimal exponent = `%c\'\n"), 
	       list->nf1.ech);
      fprintf (stderr, 
	       _("Symbol used to denote the imaginary unit = `%c\'\n\n"), 
	       list->nf1.iu);
      return -1;
    }
  else if ( !valid_numfmt(&list->nf2) )
    {
      fprintf (stderr, 
	       _("The numeric format specified for the second file is illegal,\n"));
      fprintf (stderr,
	       _("the following symbols should be all different\nwhile two or more of them are actually equal:\n"));
      fprintf (stderr, _("\nDecimal point = `%c\'\n"), list->nf2.dp);
      fprintf (stderr, _("Thousands separator = `%c\'\n"), list->nf2.thsep);
      fprintf (stderr, _("Leading positive sign = `%c\'\n"), list->nf2.pos_sign);
      fprintf (stderr, _("Leading negative sign = `%c\'\n"), list->nf2.neg_sign);
      fprintf (stderr, _("Prefix for decimal exponent = `%c\'\n"), 
	       list->nf2.ech);
      fprintf (stderr, 
	       _("Symbol used to denote the imaginary unit = `%c\'\n\n"), 
	       list->nf2.iu);
      return -1;
    }
  else
    {
      if( !(list->optmask & (_H_MASK | _V_MASK)) )
	{
	  list->file1 = (const char*) argv[optind];
	  list->file2 = (const char*) argv[optind+1];
	}
      return 0;
    }
}  
