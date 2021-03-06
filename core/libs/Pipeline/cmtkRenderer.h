/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkRenderer_h_included_
#define __cmtkRenderer_h_included_

#include <Pipeline/cmtkPipelineObject.h>
#include <Pipeline/cmtkImageRGB.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** General renderer template class.
 * This class provides a virtual "Render()" function as a common access path
 * for client code. It also provides an ImageRGB input object representing
 * the actual image to be rendered.
 */
class Renderer 
  : public PipelineObject 
{
public:
  /// Calls to this function should make derived objects update their display.
  virtual void Render();

  /// Set image to display.
  void SetInput( ImageRGB *const input );

  /// The actual Update() function.
  virtual long Update();

  /** Make this renderer active.
   *\see Active
   */
  virtual void SetActive() {Active = true; }

  /** Make this renderer inactive.
   *\see Active
   */
  virtual void SetInactive() { Active = false; }

  /** Query activity state of this renderer.
   *\see Active
   */
  virtual int IsActive() const { return Active; }

protected:
  /// Default constructor.
  Renderer();

  /// Destructor.
  ~Renderer();

  /// Image to be displayed.
  ImageRGB* Input;

private:
  /** Active flag.
   * If this flag is set, the renderer is active. It will display its 
   * associated image if Render() is called. If this flag is not
   * set, calls to Render() will result in this object setting its
   * output to size 0.
   */
  bool Active;

  /** Recursion flag.
   * This flag is used by the Render() method to check whether this viewer is
   * being rendered already. This is necessary, as checking ancestors in the
   * pipeline can lead to calls back to the environment, Tcl for example, for
   * status output. That, in turn can cause calls to the renderer again.
   *
   * So while this object is in its Render() method, this flag is set to 
   * 'true'.
   * If Render() is then called again, the function simply returns without
   * actually updating anything, thus avoiding recursion.
   */
  bool RenderPending;

  /// Convenience declaration for calls to parent class' functions.
  typedef PipelineObject Superclass;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRenderer_h_included_
