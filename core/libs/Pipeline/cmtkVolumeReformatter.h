/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkVolumeReformatter_h_included_
#define __cmtkVolumeReformatter_h_included_

#include <cmtkconfig.h>

#include <cmtkArrayFilter.h>
#include <cmtkVolumeWrapper.h>
#include <cmtkImage.h>

#include <cmtkMacros.h>
#include <cmtkReformatVolume.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class for a filter to reformat registered volume data.
 */
class VolumeReformatter : 
  public ArrayFilter<VolumeWrapper,Image,2> 
{
public:
  /// Constructor function.
  static VolumeReformatter* New() { return new VolumeReformatter; }

  /// Return virtual class name.
  virtual const char* GetClassName() const { return "VolumeReformatter"; }

  /// Index of the reference dataset.
  igsClassParameter(int,ReferenceIndex);

  /// Index of the "top" dataset for overlay and subtraction.
  igsClassParameter(int,TopIndex);

  /// Lower threshold for overlay.
  igsClassParameter(double,LowerThreshold);

  /// Upper threshold for overlay.
  igsClassParameter(double,UpperThreshold);

  /// Output image format.
  igsClassParameter(int,OutputFormat);

  /// Checkerboard filling of missing data regions.
  igsClassParameter(bool,CheckerboardMode);

  /// Path to write image files.
  igsClassParameterString(ImagePath);

  /// Pattern for naming of image files.
  igsClassParameterString(FilenamePattern);

  /** Execute function.
   * This function DOES NOT perform the actual writing of reformatted slices.
   * It rather computes the preview image that constitutes the filter output
   * of this class.
   */
  virtual void Execute();

  /** Write reformatted slices images.
   */
  void ExecuteReformat();

protected:
  /// Default constructor.
  VolumeReformatter();

  /// Virtual destructor.
  virtual ~VolumeReformatter();

private:
  /// Actual reformatter object from the Registration library.
  cmtk::ReformatVolume m_ReformatVolume;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVolumeReformatter_h_included_
