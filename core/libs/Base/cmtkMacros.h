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

#ifndef __cmtkMacros_h_included_
#define __cmtkMacros_h_included_

#include <cmtkconfig.h>

#include <stdlib.h>
#include <string.h>
#include <vector>

#include <cmtkTypes.h>

/** Trap for abstract class member functions.
 * In DEBUG mode (that is, when compiled with DEBUG defined, a call to this
 * macro will terminate the application. Without DEBUG mode, it has no effect.
 * The purpose of this is to detect calls to member functions of otherwise 
 * purely virtual classes. That is, we are trying to detect calls to member of
 * a class that should never have any instances anyway.
 */
inline void abstract() {
#ifdef DEBUG
  abort();
#endif
}

/** Macro to define a protected scalar class parameter and public read and 
 * write access member functions.
 */
#define igsGetSetMacro(type,name) \
  public: type name; \
             void Set##name( const type& v ) { name = v; } \
             const type& Get##name() const { return name; } \
             type& Get##name() { return name; } \
             void Get##name( type& to ) { to = name; }

/** Macro to define a protected scalar class parameter and public read and 
 * write access member functions.
 */
#define cmtkGetSetMacro(type,name) \
  public: type m_##name; \
             void Set##name( const type& v ) { this->m_##name = v; } \
             const type& Get##name() const { return this->m_##name; } \
             type& Get##name() { return this->m_##name; } \
             void Get##name( type& to ) { to = this->m_##name; }

/** Macro to define a protected vector parameter and public read and 
 * write access member functions.
 */
#define igsGetSetMacroVector(type,name) \
  public: std::vector<type> name; type Primary##name; \
             void Set##name( const type& v, const unsigned index = 0 ) \
               { name[index] = v; if ( !index ) Primary##name = v; } \
             const unsigned int GetCount##name() const { return name.size(); }\
             const type& Get##name( const unsigned index = 0 ) const \
               { return name[index]; } \
             type& Get##name( const unsigned index = 0 ) \
               { return name[index]; } \
             void Get##name( type& to, const unsigned index = 0 ) \
               { to = name[index]; }

/** Macro to define a protected scalar class parameter and public read and 
 * write access member functions.
 */
#define igsGetSetMacroDefault(type,name,default) \
  public: type name; \
             void Set##name( const type v = default ) { name = v; }    \
             const type& Get##name() const { return name; } \
             type& Get##name() { return name; } \
             void Get##name( type& to ) { to = name; }

/** Macro to define a protected string class parameter and public read and 
 * write access member functions.
 */
#define igsGetSetMacroString(name) \
  protected: char *name;		       \
  public:    void Set##name( const char* v ) {		   \
                  if ( (name != NULL) ) { if ( v != NULL ) \
                                          if ( !strcmp( name, v ) ) return; \
                                          free(name); name = NULL; \
                  } else { if ( v == NULL ) return; } \
                  if ( v != NULL ) { name = strdup(v); } } \
             char* Get##name() { return name; } \
             const char* Get##name() const { return name; } \
             void Get##name( char*& to ) { to = name; }

/** Macro to define a protected pointer class parameter and public read and 
 * write access member functions.
 */
#define igsGetSetMacroPtr(t,name) \
  public: t *name; \
             void Set##name( t *const v ) { name = v; } \
             t* Get##name() { return name; } \
             const t* Get##name() const { return name; } \
             void Get##name( t*& to ) { to = name; }

/** Macro to define a protected pointer class parameter and public read and 
 * write access member functions.
 */
#define igsGetSetMacroObject(t,name) \
  protected: t *name; \
  public:    void Set##name( t *const v ) { if (v) v->Reference(); if (name) name->Delete(); name = v; } \
             t* Get##name() { return name; } \
             const t* Get##name() const { return name; } \
             void Get##name( t*& to ) { to = name; }

/** Macro to define a protected 2D array class parameter and public read and 
 * write access member functions.
 */
#define igsGetSetMacro2Array(type,name) \
  public: type name[2];				    \
             void Set##name( const type v0, const type v1 ) \
                  { name[0] = v0; name[1] = v1; } \
             void SetByIndex##name( const int axis, const type v )  \
                  { name[axis] = v; } \
             void Set##name( const type* v ) \
                  { name[0] = v[0]; name[1] = v[1]; } \
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
#define igsGetSetMacro3Array(type,name) \
  public: type name[3]; \
             void Set##name( const type v0, const type v1, const type v2 ) \
                  { name[0] = v0; name[1] = v1; name[2] = v2; } \
             void SetByIndex##name( const int axis, const type v ) \
                  { name[axis] = v; } \
             void Set##name( const type* v ) \
                  { name[0] = v[0]; name[1] = v[1]; name[2] = v[2]; } \
             void Get##name( type &v0, type &v1, type &v2 ) const \
                  { v0 = name[0]; v1 = name[1]; v2 = name[2]; } \
             void Get##name( type *const v ) const \
                  { v[0] = name[0]; v[1] = name[1]; v[2] = name[2]; } \
             const type* Get##name() const { return name; } \
             type* Get##name() { return name; } \
             type Get##name( const int i ) const { return name[i]; }

/** Macro to define a protected 3D array class parameter and public read and 
 * write access member functions.
 */
#define igsGetMacro3Array(type,name) \
  public: type name[3]; \
             void Get##name( type &v0, type &v1, type &v2 ) const \
                  { v0 = name[0]; v1 = name[1]; v2 = name[2]; } \
             void Get##name( type *const v ) const \
                  { v[0] = name[0]; v[1] = name[1]; v[2] = name[2]; } \
             const type* Get##name() const { return name; } \
             type* Get##name() { return name; } \
             type Get##name( const int i ) const { return name[i]; }
//@}

#endif // #ifndef __cmtkMacros_h_included_
