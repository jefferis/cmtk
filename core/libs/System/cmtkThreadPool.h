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
//  $Revision: -1 $
//
//  $LastChangedDate: $
//
//  $LastChangedBy: $
//
*/

#ifndef __cmtkThreadPool_h_included_
#define __cmtkThreadPool_h_included_

#include <cmtkconfig.h>

#ifdef CMTK_USE_GCD
#  include "cmtkThreadPoolGCD.h"
namespace cmtk { typedef ThreadPoolGCD ThreadPool; }
#else
#  include "cmtkThreadPoolThreads.h"
namespace cmtk { typedef ThreadPoolThreads ThreadPool; }
#endif

#endif // #ifndef __cmtkThreadPool_h_included_
