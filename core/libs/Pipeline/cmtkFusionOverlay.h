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

#ifndef __cmtkFusionOverlay_h_included_
#define __cmtkFusionOverlay_h_included_

#include "cmtkArrayFilter.h"

#include "cmtkImageRGB.h"
#include "cmtkMacros.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/// Class for overlay image fusion.
class FusionOverlay : 
    public ArrayFilter<ImageRGB,ImageRGB,2> 
{
public:
  /// Create new class instance.
  static FusionOverlay* New() { return new FusionOverlay; }

  /// Return virtual class name.
  virtual const char* GetClassName() const { return "cmtkFusionOverlay"; }

  /// Parameter selecting which modality to display on top.
  igsClassParameter(int,TopImageIndex);

  /** Lower overlay threshold.
   * All pixels below this threshold in the top image are removed and reveal
   * the corresponding pixel of the bottom image.
   */
  igsClassParameter(double,LowerThreshold);

  /** Upper overlay threshold.
   * All pixels above this threshold in the top image are removed and reveal
   * the corresponding pixel of the bottom image.
   */
  igsClassParameter(double,UpperThreshold);

  /** Set pointer to original image data.
   */
  void SetOriginalImage( const int index, Image *const originalImage );

  /// Perform fusion.
  virtual void Execute();

protected:
  /// Constructor.
  FusionOverlay();

  /// Virtual destructor.
  virtual ~FusionOverlay();

private:
  /// Pointers to the original image data before RGB lookup.
  Image *OriginalImage[2];
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFusionOverlay_h_included_
