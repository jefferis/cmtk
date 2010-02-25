/*
//
//  Copyright 2004-2010 SRI International
//  Copyright 1997-2009 Torsten Rohlfing
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
//  $Revision: 1226 $
//
//  $LastChangedDate: 2010-02-08 11:28:13 -0800 (Mon, 08 Feb 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cstring>

#include "cmtkReadPGMTests.txx"
#include "cmtkSQLiteTests.txx"

/** Set up table of test names and function pointers */
typedef int (*testFuncPtr)();

typedef 
struct __testNameAndFunctionPointer
{
  const char* name;
  const testFuncPtr func;
} testNameAndFunctionPointer;

const testNameAndFunctionPointer testTable[] =
{
  { "ReadPGM",               &testReadPGM },
  { "SQLiteCreate",          &testSQLiteCreate },
  { "SQLiteOpen",            &testSQLiteOpen },
  { NULL, NULL }
};

int
main( const int argc, const char* argv[] )
{
  int testNumber = -1;
  // is test name given on command line?
  if ( argc < 2 )
    {
    // no: ask user in dialog mode
    for ( size_t i = 0; testTable[i].name; ++i )
      {
      std::cout << i << ". " << testTable[i].name << std::endl;
      }
    std::cout << "Run test number: ";
    std::cin >> testNumber;
    }
  else
    {
    // batch mode: find test by name given on command line
    for ( size_t i = 0; testTable[i].name; ++i )
      {
      if ( !std::strcmp( argv[1], testTable[i].name ) )
	testNumber = i;
      }
    }

  // run test, or return error if none found
  if ( testNumber < 0 )
    return 2;
  else
    return testTable[testNumber].func();
}
