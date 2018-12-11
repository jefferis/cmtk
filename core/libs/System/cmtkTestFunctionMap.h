/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkTestFunctionMap_h_included_
#define __cmtkTestFunctionMap_h_included_

#include <cmtkconfig.h>

#include <string>
#include <map>
#include <iostream>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Map from test name to test function.
class TestFunctionMap
{
public:
  /// Test function pointer.
  typedef int (*TestFunctionPtr)();

  /// Add test.
  void AddTest( const std::string& testName, TestFunctionPtr testFunction )
  {
    this->m_Map[testName] = testFunction;
  }

  /// Run a test by name.
  int RunTestByName( const std::string& testName )
  {
    MapType::iterator test = this->m_Map.find( testName );
    if ( test == this->m_Map.end() )
      {
      std::cerr << "Test '" << testName << "' not found.";
      return 2;
      }
    return test->second();
  }

private:
  /// Map from test name to test pointer.
  typedef std::map< std::string, TestFunctionPtr > MapType;

  /// Map from test name to test function.
  MapType m_Map;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTestFunctionMap_h_included_
