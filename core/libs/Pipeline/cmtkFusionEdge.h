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

#ifndef __cmtkFusionEdge_h_included_
#define __cmtkFusionEdge_h_included_

#include <cmtkconfig.h>

#include <cmtkArrayFilter.h>

#include <cmtkImageRGB.h>
#include <cmtkMacros.h>

#include <cmtkImageEdgeOperator.h>
#include <cmtkColormap.h>
#include <cmtkImageToImageRGB.h>
#include <cmtkFusionAlpha.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class for edge image overlay fusion.
 * As large portions of the computations to be done by this class are already
 * available from other classes (such as alpha-blending image fusion), this new
 * class simply builds a private "background-pipeline" of filters. These are
 * then used for the actual processing, FusionEdge serving only as an
 * interface to their inputs and outputs.
 */
class FusionEdge : 
  /// This is a filter from two scalar images to an RGB image.
  public ArrayFilter<Image,ImageRGB,2> 
{
public:
  /// Create new class instance.
  static FusionEdge* New() { return new FusionEdge; }

  /// Return virtual class name.
  virtual const char* GetClassName() const { return "FusionEdge"; }

  /// Perform fusion.
  virtual void Execute();

  /// Constructor.
  FusionEdge();

  /// Virtual destructor.
  virtual ~FusionEdge();

  /// Filter object that performs the actual edge computation.
  ImageEdgeOperator *EdgeOperator;

  /** Colormaps for the two original images.
   * Depending on which image the edge detection is applied to, one of these
   * colormaps is connected to the image lookup module to generate the RGB
   * representation of the other image.
   */
  Colormap *ColormapImage;

  /// Colormap for the edge image.
  Colormap *ColormapEdge;

  /// Colormap lookup filter for the original bottom image.
  ImageToImageRGB *LookupImage;

  /// Colormap lookup filter for the edge image.
  ImageToImageRGB *LookupEdge;

  /** Alpha-blending fusion filter.
   * The actual edge image overlay is identical to the image fusion technique
   * implemented by the "FusionAlpha" class. Therefore, an instance of that
   * class is used as a "slave filter" to produce the output RGB image of the
   * edge image overlay. The output object of the "FusionAlpha" instance is 
   * simply forwarded to clients of the edge overlay class. This is achieved by
   * simply setting this class's "Output" field to the pointer returned by
   * FusionAlpha->GetOutput().
   */
  FusionAlpha *m_FusionAlpha;

protected:
  /// Update internal links after input links changed.
  virtual void InputLinkChanged( const int inputIndex );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFusionEdge_h_included_
