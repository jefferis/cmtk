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

#ifndef __cmtkSafeCounterGCD_h_included_
#define __cmtkSafeCounterGCD_h_included_

#include <cmtkconfig.h>

#include <dispatch/dispatch.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Thread-safe counter using Grand Central Dispatch FIFO queue.
 */
class SafeCounterGCD
{
public:
  /// Constructor.
  SafeCounterGCD( const unsigned int counter = 0 ) : 
    m_Counter( counter ), 
    m_Queue( dispatch_queue_create( "SafeCounterGCD", NULL) ) 
  {}

  /// Destructor.
  ~SafeCounterGCD()
  {
    dispatch_release( this->m_Queue );
  }
  
  /// Retrieve counter value.
  unsigned int Get() const;

  /// Increment and return new counter value.
  unsigned int Increment();
  
  /// Decrement and return new counter value.
  unsigned int Decrement();
  
private:
  /// The actual counter.
  unsigned int m_Counter;
  
  /// GCD for thread-safe exclusive access to counter.
  mutable dispatch_queue_t m_Queue;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSafeCounterGCD_h_included_
