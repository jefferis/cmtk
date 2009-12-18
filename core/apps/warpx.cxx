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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkconfig.h>

#include <stdio.h>
#include <string.h>

#include <cmtkImagePairNonrigidRegistrationCommandLine.h>
#include <cmtkProgressConsole.h>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace warpx
{
#endif
int
main ( int argc, char *argv[] )
{
  try 
    {
    cmtk::ImagePairNonrigidRegistrationCommandLine Registration( argc, argv );
    
    // set up console progress reporting
    cmtk::ProgressConsole progressInstance( "BSplineImageRegistration" );
    Registration.Register();
    }
  catch ( cmtk::ImagePairRegistration::ConstructorFailed ) 
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
} // namespace warpx
} // namespace apps
} // namespace cmtk
#endif

