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

#ifndef __cmtkThreadParameters_h_included_
#define __cmtkThreadParameters_h_included_

#include "System/cmtkThreadSystemTypes.h"

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Base class for thread parameter blocks.
 * This doesn't hurt to have, even if we're building without thread support.
 */
template<class T> 
class ThreadParameters
{
public:
  /// Template pointer to the object starting this thread.
  T* thisObject;
  /// Unique index of this thread instance among all threads.
  unsigned int ThisThreadIndex;
  /// Total number of threads created.
  unsigned int NumberOfThreads;

  /// Thread ID of a started thread.
  ThreadIDType m_ThreadID;

#ifdef _MSC_VER
  void* m_Handle;
#endif

  /// Default constructor.
  ThreadParameters() 
  {
    this->m_ThreadID = 0;
#ifdef _MSC_VER
    m_Handle = NULL;
#endif
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkThreadParameters_h_included_
