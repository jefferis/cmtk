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

#ifndef __cmtkSafeCounter_h_included_
#define __cmtkSafeCounter_h_included_

#include <cmtkconfig.h>

#ifdef CMTK_USE_GCD
#  include "cmtkSafeCounterGCD.h"
namespace cmtk { typedef SafeCounterGCD SafeCounter; }
#else
#  include "cmtkSafeCounterMutex.h"
namespace cmtk { typedef SafeCounterMutex SafeCounter; }
#endif // #ifdef CMTK_USE_GCD


#endif // #ifndef __cmtkSafePtr_h_included_
