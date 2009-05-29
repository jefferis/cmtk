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

#ifndef __cmtkObject_h_included_
#define __cmtkObject_h_included_

#include <cmtkconfig.h>

#ifndef NULL
#define NULL 0
#endif

/** Macro to define a protected scalar class parameter and public read and 
 * write access member functions.
 */
#define igsClassParameter(type,name) \
  public: type name; \
             void Set##name( const type v ) { if ( name != v ) this->UpdateModifiedTime(); name = v; } \
             type Get##name() const { return name; } \
             void Get##name( type& to ) { to = name; }

/** Macro to define a protected scalar class parameter and public read and 
 * write access member functions.
 */
#define cmtkClassParameter(type,name) \
  public: type m_##name; \
             void Set##name( const type v ) { if ( this->m_##name != v ) this->UpdateModifiedTime(); this->m_##name = v; } \
             type Get##name() const { return this->m_##name; } \
             void Get##name( type& to ) { to = this->m_##name; }

/** Macro to define a protected string class parameter and public read and 
 * write access member functions.
 */
#define igsClassParameterString(name) \
  public: char *name; \
             void Set##name( const char* v ) { \
                this->SetParameterString( name, v ); } \
             char* Get##name() { return name; } \
             const char* Get##name() const { return name; } \
             void Get##name( char*& to ) { to = name; }

/** Macro to define a protected pointer class parameter and public read and 
 * write access member functions.
 */
#define igsClassParameterPtr(t,name) \
  protected: t *name; \
  public:    void Set##name( t *const v ) { if ( name != v ) this->UpdateModifiedTime(); name = v; } \
             t* Get##name() { return name; } \
             const t* Get##name() const { return name; } \
             void Get##name( t*& to ) { to = name; }

/** Macro to define a protected pointer class parameter and public read and 
 * write access member functions.
 */
#define igsClassParameterObject(t,name) \
  protected: t *name; \
  public:    void Set##name( t *const v ) { if ( name != v ) this->UpdateModifiedTime(); if (v) v->Reference(); if (name) name->Delete(); name = v; } \
             t* Get##name() { return name; } \
             const t* Get##name() const { return name; } \
             void Get##name( t*& to ) { to = name; }

/** Macro to define a protected 2D array class parameter and public read and 
 * write access member functions.
 */
#define igsClassParameter2Array(type,name) \
  protected: type name[2]; \
  public:    void Set##name( const type v0, const type v1 ) \
                  { if ( (name[0] != v0) || (name[1] != v1) ) this->UpdateModifiedTime(); name[0] = v0; name[1] = v1; } \
             void SetByIndex##name( const int axis, const type v ) \
                  { if (name[axis] != v) this->UpdateModifiedTime(); name[axis] = v; } \
             void Set##name( const type* v ) \
                  { if ( (name[0] != v[0]) || (name[1] != v[1]) ) this->UpdateModifiedTime(); name[0] = v[0]; name[1] = v[1]; } \
             void Get##name( type &v0, type &v1 ) const \
                  { v0 = name[0]; v1 = name[1]; } \
             void Get##name( type *const v ) const \
                  { v[0] = name[0]; v[1] = name[1]; } \
             const type* Get##name() const { return name; } \
             type* Get##name() { return name; } \
             type Get##name( const int i ) const { return name[i]; }

/** Macro to define a protected 3D array class parameter and public read and 
 * write access member functions.
 */
#define igsClassParameter3Array(type,name) \
  protected: type name[3]; \
  public:    void Set##name( const type v0, const type v1, const type v2 ) \
                  { if ( (name[0] != v0) || (name[1] != v1) || (name[2] != v2) ) this->UpdateModifiedTime(); name[0] = v0; name[1] = v1; name[2] = v2; } \
             void SetByIndex##name( const int axis, const type v ) \
                  { if (name[axis] != v) this->UpdateModifiedTime(); name[axis] = v; } \
             void Set##name( const type* v ) \
                  { if ( (name[0] != v[0]) || (name[1] != v[1]) || (name[2] != v[2]) ) this->UpdateModifiedTime(); name[0] = v[0]; name[1] = v[1]; name[2] = v[2]; } \
             void Get##name( type &v0, type &v1, type &v2 ) const \
                  { v0 = name[0]; v1 = name[1]; v2 = name[2]; } \
             void Get##name( type *const v ) const \
                  { v[0] = name[0]; v[1] = name[1]; v[2] = name[2]; } \
             const type* Get##name() const { return name; } \
             type* Get##name() { return name; } \
             type Get##name( const int i ) const { return name[i]; }

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Static counter representing the current (discrete) processing time.
 * This variable is incremented with every query using 
 * Object::GetCurrentTime(). Direct access MUST NOT occur in order 
 * to prevent inconsistent object times.
 *@name igsCurrentTime
 *@see PipelineObject#GetCurrentTime
 */
extern long CurrentTime;

/** Base class for all reference counted objects.
 * This class to some extent shadows VTK's respective class. However, our
 * class should make be a little more runtime- and memory-efficient.
 */
class Object 
{
public:
  /** Register another object as this objects owner.
   * The reference counter of this object is also incremented.
   *@param owner The object to be registered as the owner of this object. 
   * If this parameter is not given, the current owner is left untouched. In
   * this case, only the reference counter is modified.
   *@return The new value of the reference counter.
   *@see ReferenceCount
   */
  virtual int Reference() const 
  {
    return ++ReferenceCount;
  }

  /** Destroy this object.
   * This function decrements this object's reference counter. If the updated
   * reference count is zero, this object is deleted.
   */
  virtual void Delete() 
  {
    if ( (--ReferenceCount) <= 0 ) delete this;
  }
  
  /** Directly set the reference counter.
   * This function should be used CAREFULLY.
   */
  void SetReferenceCount ( const int referenceCount ) 
  {
    ReferenceCount = referenceCount;
  }
  
  /** Get the reference counter.
   */
  int GetReferenceCount () const 
  { 
    return ReferenceCount;
  }
  
  /// Return this objects last modification time.
  long GetModifiedTime() const
  { 
    return ModifiedTime;
  }
  
  /// Set time of last modification to current time.
  void UpdateModifiedTime() 
  {
    ModifiedTime = this->GetCurrentTime(); 
  }
  
  /** Explicitly set time of last modification.
   * Only monotone updates, ie. setting the update time to a more recent time,
   * are allowed.
   */
  void UpdateModifiedTime( long modifiedTime ) 
  { 
    if ( modifiedTime > ModifiedTime ) 
      ModifiedTime = modifiedTime;
  }
  
  /** Utility function: Replace one reference counted object by another.
   * References are updated accordingly. Substitution of a pointer by itself
   * does not have any effect. The function is NULL safe.
   *@param C Template parameter: The type of the objects pointed to by both
   * function parameters.
   *@param to Reference to a pointer to be replaced by the pointer given as
   * "from" parameter. If both pointers are different, the reference counter
   * of this object is decremented before overwriting the pointer. It is safe
   * to pass references to NULL pointers here.
   *@param from The reference counter of the object pointed to by this pointer
   * is increased if a pointer substitution takes place. It is safe to pass a
   * NULL pointer here.
   *@return "True" is returned if there was an actual replacement. If an object
   * was basically replaced by itself, "false" is returned instead. Note that
   * only POINTERS are compared to find out which of both has happened.
   */
  template<class C> bool ReplaceObject( C*& to, C *const from ) 
  {
    if ( from == to ) return false;
    if ( from ) from->Reference();    
    if ( to ) to->Delete();
    to = from;
    this->UpdateModifiedTime();
    return true;
  }

  /** Replace string parameter.
   * This function replaces a string parameter of an arbitrary derived class.
   * It is safe in that it will accept any combination of NULL and non-NULL
   * pointers for previous and new value of the respective parameter.
   * If a change actually occurred, the derived object's UpdateModifiedTime()
   * function is called.
   *@param to Reference to the char pointer to be set to the new value.
   *@param from Pointer to the new value to be assigned to "to". Assignment
   * will be done by a call to "strdup", so the caller is free to use and
   * deallocate the source string after return from this function.
   *@see PipelineObject#UpdateModifiedTime
   */
  void SetParameterString( char*& to, const char *from );

  /** Default constructor.
   * Set the reference counter to zero initially.
   */
  Object() 
  { 
    ReferenceCount = 1; ModifiedTime = 0;
  }

  /** Virtual constructor.
   * This function is only present to ensure that derived classes can have a
   * virtual destructor function.
   */
  virtual ~Object() {};

  /** Query the time processing counter.
   * Immediately after returning the present time, the time counter is
   * advanced.
   *@return The present time.
   *@see igsCurrentTime
   */
  static long GetCurrentTime () 
  { 
    return CurrentTime++; 
  }
  
private:
  /// The actual reference counter.
  mutable int ReferenceCount;

  /** Last modification time.
   * This is the time of the last modification made to this object or one of
   * its input parameters.
   */
  long ModifiedTime;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkObject_h_included_
