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

#ifndef __cmtkFusionIsolines_h_included_
#define __cmtkFusionIsolines_h_included_

#include <cmtkconfig.h>

#include <cmtkArrayFilter.h>

#include <cmtkImageRGB.h>
#include <cmtkMacros.h>

#include <cmtkIsolineFilter.h>
#include <cmtkColormap.h>
#include <cmtkImageToImageRGB.h>
#include <cmtkFusionOverlay.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class for edge image overlay fusion.
 * As large portions of the computations to be done by this class are already
 * available from other classes (such as alpha-blending image fusion), this new
 * class simply builds a private "background-pipeline" of filters. These are
 * then used for the actual processing, FusionIsolines serving only as an
 * interface to their inputs and outputs.
 */
class FusionIsolines : 
  /// Inherit filter functions.
  public ArrayFilter<Image,ImageRGB,2> 
{
public:
  /// Create new class instance.
  static FusionIsolines* New() { return new FusionIsolines; }

  /// Parameter selecting which modality to display on top.
  int IsolineImageIndex;

  /// Select which image to apply the edge operator to.
  void SetIsolineImageIndex( const int isolineImageIndex );

  /// Return the index of the image that the edge operator is applied to.
  int GetIsolineImageIndex() const { return IsolineImageIndex; }

  /// Perform fusion.
  virtual void Execute();

protected:
  /// Constructor.
  FusionIsolines();

  /// Virtual destructor.
  virtual ~FusionIsolines();

  /// Filter object that performs the actual edge computation.
  IsolineFilter *m_IsolineFilter;

  /** Colormaps for the two original images.
   * Depending on which image the edge detection is applied to, one of these
   * colormaps is connected to the image lookup module to generate the RGB
   * representation of the other image.
   */
  Colormap *ColormapImage[2];

  /// Colormap lookup filter for the original bottom image.
  ImageToImageRGB *LookupImage;

  /// Colormap lookup filter for the edge image.
  ImageToImageRGB *LookupIsolines;

  /** Overlay fusion filter.
   * The actual isoline image overlay is identical to the image fusion 
   * technique implemented by the "cmtkFusionOverlay" class. Therefore, an 
   * instance of that class is used as a "slave filter" to produce the output 
   * RGB image of the isoline overlay. The output object of the "FusionOverlay"
   * instance is simply forwarded to clients of the edge overlay class. This
   * is achieved by simply setting this class's "Output" field to the pointer
   * returned by FusionOverlay->GetOutput().
   */
  FusionOverlay *m_FusionOverlay;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFusionIsolines_h_included_
