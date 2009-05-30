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

#ifndef __cmtkPipelineObject_h_included_
#define __cmtkPipelineObject_h_included_

#include <cmtkconfig.h>

#include <cmtkObject.h>

#include <stdlib.h>
#include <stdio.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Global flag for pipeline debugging.
 * If this flag is set, virtual functions of derived classes will generate
 * debugging output when they are invoked. This flag is only present and
 * effective if the library was compiled with the DEBUG symbol defined.
 */
static bool PipelineDebugMode = 0;

/** Base class for all reference counted and linked objects.
 */
class PipelineObject : 
  /// Inherit reference-counting object.
  public Object 
{
protected:
  /** This object's owner.
   * The owner is the object queried for updates first when this object is 
   * asked to update itself.
   */
  PipelineObject *Owner;
  
public:
  virtual const char *GetClassName() const { return "PipelineObject"; }

  //  void SetOwner( PipelineObject *const owner ) { Owner = owner; }
  const PipelineObject* GetOwner() const { return Owner; }

  /** Register another object as this objects owner.
   * The reference counter of this object is also incremented.
   *@param owner The object to be registered as the owner of this object. 
   * If this parameter is not given, the current owner is left untouched. In
   * this case, only the reference counter is modified.
   *@return The new value of the reference counter.
   *@see ReferenceCount
   */
  int Register( PipelineObject *const owner = NULL );

  /** Unregister one owner object.
   * This function decrements this object's reference counter. If the updated
   * reference count is zero, this object is destroyed.
   */
  void Unregister( PipelineObject *const owner = NULL );

  /** Check for update.
   * This function first checks whether since its last execution its owner
   * has been modified. In this case, the Execute() function is called to
   * update the current object with the new input. Derived classes may override
   * this function if they have more than one input object, for instance.
   *
   * Such derived implementations can then use the CheckForUpdate() and
   * ExecuteIfNecessary() member functions for convenient state checking and
   * execution.
   *@see Execute
   *@see CheckInputForUpdate
   *@see ExecuteIfNecessary
   */
  virtual long Update ();

  /** Execute the current object.
   * Derived classes need to override this function in order to make the
   * respective instance up-to-date.
   */
  virtual void Execute () { this->UpdateExecuteTime(); }

  /// Set global flag for pipeline debugging.
  static void SetDebugMode( const bool debugMode = true ) 
  {
    PipelineDebugMode = debugMode;
  }

protected:
  /** Default constructor.
   * Set the reference counter to zero initially.
   */
  PipelineObject();

  /** Destructor.
   * This is defined virtual so that derived classes are enabled to provide
   * their own virtual destructor functions.
   */
  virtual ~PipelineObject() {};

  /// Set time of last execution to current time.
  void UpdateExecuteTime() 
  { 
    ExecuteTime = this->GetCurrentTime(); 
    ExecutePending = 0; 
  }

  /** Compare input for update.
   * For the given input object (Owner for example), the Update() function
   * is called. Afterwards, the returned execution time of the input object
   * is compared to the current object's modification time. The later of both
   * times is then set as the current object's modification time.
   */
  virtual int CheckInputForUpdate( PipelineObject *const object );

  /** Execute an update if object was modified after last execution.
   *@return The new time of last execution.
   */
  virtual long ExecuteIfNecessary();

private:
  /** Last execution time.
   * This is the time of the latest execution of this objects Execute() 
   * function, ie. the time when this objects state or output was last
   * updated according to the input parameters.
   */
  long ExecuteTime;

  /** Flag for pending updates.
   * This field is set to 1 if Update() discovers a changed input object.
   * ExecuteIfNecessary() then evaluates this flag in addition to this objects
   * modification time and calls Exedcute() even if only this flag is set.
   * UpdateExecuteTime() then finally resets this field to 0.
   */
  int ExecutePending;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkPipelineObject_h_included_
