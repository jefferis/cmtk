/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkconfig.h>

#include <stdio.h>
#include <string.h>

#include <cmtkTimers.h>
#include <cmtkElasticRegistrationCommandLine.h>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace warp
{
#endif
int
main ( int argc, char *argv[] )
{
  try 
    {
    cmtk::ElasticRegistrationCommandLine Registration( argc, argv );
    
    const double baselineTime = cmtk::Timers::GetTimeProcess();
    Registration.Register();
    const int elapsed = static_cast<int>( cmtk::Timers::GetTimeProcess() - baselineTime );
    
    fprintf( stderr, "Time: %d (%f : %f) sec.\n", elapsed, (float)Registration.GetTotalElapsedTime(), (float)Registration.GetThreadTotalElapsedTime() );
    }
  catch ( cmtk::VoxelRegistration::ConstructorFailed ) 
    {
    return 1;
    }
  catch (...) 
    {
    fputs( "Uncaught exception during registration.", stderr );
    return 1;
    }
  
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace warp
} // namespace apps
} // namespace cmtk
#endif

