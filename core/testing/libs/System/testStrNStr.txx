/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <System/cmtkStrUtility.h>
#include <System/cmtkConsole.h>

int
testStrNStr()
{
  const char* needle = "needle";
  
  const char* haystack1 = "this is the needle";
  // make sure we find needle
  if ( cmtk::StrNStr( haystack1, strlen( haystack1 ), needle ) == NULL )
    {
    cmtk::StdErr << "StrNStr test #1 failed\n";
    return 1;
    }

  // make sure we stop looking after nBytes needle
  if ( cmtk::StrNStr( haystack1, strlen( haystack1 ) - 1, needle ) != NULL )
    {
    cmtk::StdErr << "StrNStr test #2 failed\n";
    return 1;
    }

  // make sure we find only needle
  const char* haystack2 = "this is not the Needle";
  if ( cmtk::StrNStr( haystack2, strlen( haystack2 ), needle ) != NULL )
    {
    cmtk::StdErr << "StrNStr test #3 failed\n";
    return 1;
    }

  // make sure we find only needle
  const char* haystack3 = "this is not the needl either";
  if ( cmtk::StrNStr( haystack3, strlen( haystack3 ), needle ) != NULL )
    {
    cmtk::StdErr << "StrNStr test #4 failed\n";
    return 1;
    }

  // make sure we find only needle
  const char haystack4[] = "first put \x00 then put needle";
  if ( cmtk::StrNStr( haystack4, sizeof( haystack4 ), needle ) == NULL )
    {
    cmtk::StdErr << "StrNStr test #5 failed\n";
    return 1;
    }

  // make sure we can find prefixes
  const char* haystack5 = "needle first";
  if ( cmtk::StrNStr( haystack5, strlen( haystack5 ), needle ) == NULL )
    {
    cmtk::StdErr << "StrNStr test #6 failed\n";
    return 1;
    }
  
  return 0;
}

