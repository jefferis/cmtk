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

#ifndef __cmtkRendererCollection_h_included_
#define __cmtkRendererCollection_h_included_

#include <cmtkconfig.h>

#include <cmtkObject.h>
#include <cmtkRenderer.h>

#include <list>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class to collect all active render views in the current application.
 * All renderers should register themselves with the global instance of this
 * class upon initialization.
 *@see RendererCollectionInstance
 */
class RendererCollection : 
  public Object
{
public:
  /// Create new renderer collection object.
  static RendererCollection *New() { return new RendererCollection; }

  /// Return virtual class name.
  virtual const char* GetClassName() { return "igsRenderCollection"; }

  /** Add a new renderer to the renderer collection.
   * The reference counter of the added renderer object is incremented
   * and the object pointer is added to the end of the renderer list.
   *@param renderer Pointer to the renderer to add.
   */
  void AddRenderer( Renderer *const renderer );

  /** Remove a renderer from the renderer collection.
   *@param renderer Pointer to the renderer to be removed from the collection.
   * It is safe to pass a pointer to an object that is not actually IN the
   * collection.
   */
  void RemoveRenderer( Renderer *const renderer );

  /** Render all active views.
   * This function walks through all renderers in the collection from front
   * to back of the renderer list. For each renderer, its "Render()" function
   * is called so that the respective object updates its view.
   *@see Renderer#Render
   */
  virtual void Render();

protected:
  /// Default constructor.
  RendererCollection();

  /** Destructor.
   * Clear the renderer list and unregister from all stored objects.
   */
  virtual ~RendererCollection();

private:
  /// Define a list as the container type for the renderers.
  typedef std::list<Renderer*> RendererList;

  /// The list of renderers.
  RendererList m_RendererList;
};

/// The global instance of the application-wide renderer collection.
extern RendererCollection *RendererCollectionInstance;

/** Call the global collections "Render()" function.
 * This function may be called by the application in order to update all render
 * views at virtually the same time.
 *@see RendererCollection#Render
 */
int igsRenderAll();

//@}

} // namespace cmtk

#endif
