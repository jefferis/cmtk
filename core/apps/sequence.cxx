/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <Base/cmtkValueSequence.h>
#include <Base/cmtkTypes.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <cfloat>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <list>
#include <vector>
#include <limits>

double MaxThreshold = 0;
bool UseMaxThreshold = false;
bool AbsoluteValues = false;

const char* OutputFormat = "%.6f";

int
doMain( const int argc, const char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Value sequence" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Analyze sequence of numerical values, which is read from standard input" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( 't', "thresh" ), &MaxThreshold, "Maximum value threshold. All values above are ignored.", &UseMaxThreshold );
    cl.AddSwitch( Key( 'a', "abs" ), &AbsoluteValues, true, "Use absolute values." );
    cl.AddOption( Key( 'f', "format" ), &OutputFormat, "Output number format in printf() style." );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << "Command line parse error: " << e << "\n";
    return 1;
  }

  cmtk::ValueSequence<double> seq;
  std::list<double> list;

  unsigned int countOverThreshold = 0;
  double f;
  while ( ! std::cin.eof() ) 
    {
    std::cin >> f; 

    if ( ! finite( f ) ) break;

    if ( AbsoluteValues ) f = fabs( f );
    
    if ( UseMaxThreshold && (f > MaxThreshold) )
      ++countOverThreshold;
    else 
      {
      seq.Proceed( f );
      list.push_back( f );
      }
    
    f = std::numeric_limits<double>::signaling_NaN();
    }
  
  const size_t totalNumberOfValues = seq.GetNValues() + countOverThreshold;
  printf( "Number of Values:\t%d\n", (int)totalNumberOfValues );
  printf( "Values over Threshold:\t%d (%.2f%%)\n", countOverThreshold, 100.0 * countOverThreshold / totalNumberOfValues );

  char format[120];
  snprintf( format, sizeof( format ), "\nSTAT\tMin\tMax\tMean\tStdDev\nSTATval\t%s\t%s\t%s\t%s\n", OutputFormat, OutputFormat, OutputFormat, OutputFormat );
  printf( format, seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(), sqrt( seq.GetVariance() ) );

  list.sort();
  std::vector<double> sorted;
  for ( std::list<double>::const_iterator it = list.begin(); it != list.end(); ++it )
    {
    sorted.push_back( *it );
    }
  
  printf( "\nPERC" );
  const int percentiles[] = { 5, 10, 25, 50, 75, 90, 95, -1 };
  for ( size_t idx = 0; percentiles[idx] > 0; ++idx )
    {
    printf( "\t%d", percentiles[idx] );
    }

  snprintf( format, sizeof( format ), "\t%s", OutputFormat );
  printf( "\nPERCval" );
  for ( size_t idx = 0; percentiles[idx] > 0; ++idx )
    {
    if ( percentiles[idx] == 50 )
      {
      const size_t medianIdx = sorted.size() / 2;
      if ( sorted.size() & 1 )
	{
	printf( format, sorted[medianIdx] );
	}
      else
	{
	printf( format, 0.5*(sorted[medianIdx] + sorted[1+medianIdx]) );
	}
      }
    else
      {								       
      printf( format, sorted[ (size_t)( sorted.size() * 0.01 * percentiles[idx]) ] );
      }
    }
  printf( "\n" );
    
  return 0;
}

#include "cmtkSafeMain"
