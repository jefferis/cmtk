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

#ifndef __cmtkRenderer_h_included_
#define __cmtkRenderer_h_included_

#include "Pipeline/cmtkPipelineObject.h"
#include "Pipeline/cmtkImageRGB.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Export color mode.
 */
typedef enum {
  /// Automatically determine if image is color or greyscale.
  EXPORTCOLOR_AUTO,
  /// Always export color image.
  EXPORTCOLOR_RGB,
  /// Always export greyscale image.
  EXPORTCOLOR_GREY
} ExportColorMode;

/** General renderer template class.
 * This class provides a virtual "Render()" function as a common access path
 * for client code. It also provides an ImageRGB input object representing
 * the actual image to be rendered.
 */
class Renderer 
  : public PipelineObject 
{
public:
  /// Export color mode.
  cmtkClassParameter(ExportColorMode,ExportColorMode);

  /// Image compression mode.
  igsClassParameter(int,CompressionMode);

  /// Calls to this function should make derived objects update their display.
  virtual void Render();

  /// Set image to display.
  void SetInput( ImageRGB *const input );

  /// The actual Update() function.
  virtual long Update();

  /** Make this renderer active.
   *@see Active
   */
  virtual void SetActive() {Active = true; }

  /** Make this renderer inactive.
   *@see Active
   */
  virtual void SetInactive() { Active = false; }

  /** Query activity state of this renderer.
   *@see Active
   */
  virtual int IsActive() const { return Active; }

  /** Write the input image to an RGB PPM file.
   *@param filename Output filename.
   *@param cmtc Number of descriptive comment strings to write into the output
   * file. These comments may be used to add additional information such as
   * pixel size, image history etc. to the image file.
   *@param cmtv Array of pointers to the comment strings to be written to the
   * input file. This array has to contain at least "cmtc" valid entries.
   */
  virtual void WritePPM( const char* filename, const int cmtc = 0, const char** cmtv = NULL );

protected:
  /// Default constructor.
  Renderer();

  /// Destructor.
  ~Renderer();

  /// Register to global list of all active renderers.
  virtual void RegisterToCollection();

  /// Image to be displayed.
  ImageRGB* Input;

  /** Capture displayed RGB image.
   * By default, the captured image is only the input image. However, derived
   * classes can capture actual displays which may contain annotations etc.
   * If the input image is returned, its reference counter is incremented.
   * Callers to this funciton can and should therefore call "Delete()" once
   * they are done with the object returned by this function.
   */
  virtual ImageRGB* CaptureDisplay() 
  {
    if ( Input ) Input->Reference();
    return Input;
  }

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

  /// Determine effective color mode.
  ExportColorMode GetEffectiveExportColorMode( const ImageRGB* capture ) const;

  /// Convenience declaration for calls to parent class' functions.
  typedef PipelineObject Superclass;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRenderer_h_included_
