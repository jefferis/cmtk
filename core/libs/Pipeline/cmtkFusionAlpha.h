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

#ifndef __cmtkFusionAlpha_h_included_
#define __cmtkFusionAlpha_h_included_

#include <cmtkconfig.h>

#include <cmtkArrayFilter.h>

#include <cmtkImageRGB.h>
#include <cmtkMacros.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/// Class for alpha-blending image fusion.
class FusionAlpha : 
  /// This is a filter from two RGB images to an RGB image.
  public ArrayFilter<ImageRGB,ImageRGB,2> 
{
private:
  /// Direct base class.
  typedef ArrayFilter<ImageRGB,ImageRGB,2> Superclass;

public:
  /// Transperency study (optional).
  igsClassParameterObject(ImageRGB,TransparencyImage);

  /// Create new class instance.
  static FusionAlpha* New() { return new FusionAlpha; }

  /// Return virtual class name.
  virtual const char* GetClassName() const { return "FusionAlpha"; }

  /// Parameter selecting which modality to display on top.
  igsClassParameter(int,TopImageIndex);

  /// Flag whether to generate RGB or RGBA output.
  igsClassParameter(bool,OutputHasAlpha);

  /// Perform fusion.
  virtual void Execute();

  /** Update this object.
   * Check for changes in the TransparencyImage object first, then call
   * inherited Update() function.
   *@see Object#Update
   */
  virtual long Update () {
    this->CheckInputForUpdate( TransparencyImage );
    return this->Superclass::Update();
  }

protected:
  /// Constructor.
  FusionAlpha() : TransparencyImage( NULL )
  { TopImageIndex = 0; OutputHasAlpha = false; };

  /// Virtual destructor.
  virtual ~FusionAlpha() {};
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFusionAlpha_h_included_
