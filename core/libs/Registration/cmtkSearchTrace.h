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

#ifndef __cmtkSearchTrace_h_included_
#define __cmtkSearchTrace_h_included_

#include <cstdlib>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for traces in the optimization search space.
 * An object of this class runs a list of those locations in search space that
 * have already been visited during the optimization. This may save evaluations
 * of the target function when these locations are approached again later in
 * the optimum search. The type of relative vectors is defined by a template
 * argument R. For one step size only, the default setting of "short", even 
 * "signed char" may be enough. To keep track of multiresolution schemes,
 * however, we shall need to use "float" or "double".
 */
template<class R = short>
class SearchTrace 
{
private:
  typedef struct _TraceListEntry 
  {
    /// Vector pointing to a relative location that has already been visited.
    R *RelativePosition;

    /// Value of the target function that was found there.
    double FunctionValue;

    /// Pointer to the next entry in the linked list.
    struct _TraceListEntry *Next;
  } TraceListEntry;

  /** Dimension of the search space.
   */
  int DOF;

  /** Pointer to the first element of the list of visited locations.
   */
  TraceListEntry *List;

  /** Compare an item from the track list to a given location.
   * As with every move in the search space, relative locations of the
   * previously seen samples are updated, all components of the tested
   * location must be zero except for the one giving the current search
   * direction.
   *@param entry The entry in the track list to be tested.
   *@param dir Direction, i.e. index of the parameter, in which we are moving.
   *@param step Size of the intended step.
   *@return 1, if the list entry pointed to by "entry" is the location we
   * would be in when making the step defined by "dir" and "step", 0 otherwise.
   */
  int IsHit ( const TraceListEntry* entry, const int dir, const R step ) const 
  {
    for ( int idx=0; idx<DOF; ++idx )
      if ( entry->RelativePosition[idx] && ( (dir != idx) || (entry->RelativePosition[idx] != step) ) )
	return 0;
    
    return 1;
  }
  
public:
  /** Constructor.
   * Set dimension of search space and initialize trace list.
   */
  SearchTrace ( const int _DOF ) {
    DOF = _DOF;
    List = NULL;
  }

  /** Destructor.
   * Call Clear() to remove list from memory.
   */
  ~SearchTrace () { Clear(); }

  /** Add a location to the trace list.
   *@param value The value of the target function at the location to be
   * added to the list.
   *@param dir Direction of the location to add with respect to the current
   * position in search space. This is the index of the parameter we are 
   * modifying.
   *@step Size of the step, ie. distance of the new location from the current
   * position in search space.
   */
  void Add ( const double value, const int dir = 0, const R step = 0 ) 
  {
    TraceListEntry *add = new TraceListEntry;
    add->RelativePosition = Memory::AllocateArray<R>( DOF );
    memset( add->RelativePosition, 0, sizeof(R) );
    add->RelativePosition[dir] += step;
    add->FunctionValue = value;
    add->Next = List;
    List = add;
  }

  /** Get a previously visited location from the list.
   *@param value This reference is used to return the target function's value
   * at the location that was asked for. This value is only valid, if the
   * location was in the list. In case the function returns 0, value is
   * undefined.
   *@param dir Direction in search space towards the location we ask for.
   *@param step Size of the step to make in the given direction.
   *@return 1 if the desired location was in the list, 0 otherwise.
   */
  int Get ( double& value, const int dir = 0, const R step = 0 ) const 
  {
    TraceListEntry *cursor = List;
    while ( cursor ) 
      {
      if ( IsHit( cursor, dir, step ) ) 
	{
	value = cursor->FunctionValue;
	return 1;
	}
      cursor = cursor ->Next;
      }
    return 0;
  }

  /** Move current position in search space.
   * All entries in the table of visited positions are modified accordingly to
   * keep their relative positions up-to-date.
   *@param dir Parameter modified to do the move.
   *@param step Size of the update step in the direction defined by 'dir'.
   */
  void Move ( const int dir, const R step ) 
  {
    TraceListEntry *cursor = List;
    while ( cursor ) 
      {
      cursor->RelativePosition[dir] -= step;
      cursor = cursor ->Next;
      }    
  }

  /** Clear list from memory.
   * All list entries as well as the location vectors stored in them are 
   * deleted and the list pointer is reset to NULL.
   */
  void Clear () 
  {
    while ( List ) 
      {
      TraceListEntry *save = List->Next;
      delete[] List->RelativePosition;
      delete List;
      List = save;
      }
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSearchTrace_h_included_
