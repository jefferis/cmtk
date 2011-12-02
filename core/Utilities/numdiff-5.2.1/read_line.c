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

char* read_line (FILE* pf, int* errcode)
{
  char buffer[BUFF_SIZE];
  char *ptr, *line = NULL;
  size_t lline = 0;

  while ((ptr = fgets (buffer, BUFF_SIZE,pf)))
    {
      lline += strlen (buffer);
      if (!line)
	ptr = (char*) calloc (lline + 1, sizeof(char));
      else
	ptr = (char*) realloc ((void*)line, (lline + 1) * sizeof(char));
      if (!ptr)
	{
	  if ((line))
	    free ((void*)line);
	  *errcode = OUT_OF_MEMORY;
	  return NULL;
	}
      else
	{
	  line = ptr;
	  strcat (line, buffer);
	}
      if (lline > 0 && line[lline-1] == '\n')
	break;
    }
  if (!ptr)
    {
      if ( (ferror(pf)) )
	*errcode = READING_ERROR;
      else if (lline > 0)
	*errcode = LINE_INTERR;
      else
	*errcode = EOF_REACHED; 
    }
  else
    *errcode = OK;
  return line;
} 
