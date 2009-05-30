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

#ifndef __cmtkVolumeToVtkStructuredPoints_h_included__
#define __cmtkVolumeToVtkStructuredPoints_h_included__

#include <cmtkconfig.h>

#include <vtkStructuredPointsSource.h>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>
#include <cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

/** Class to encapsulate IGS' "Volume" data into vtk's structured points.
 * Only instances of UniformVolume can be encapsulated, as vtkStructuredPoints
 * is limited to representing uniformly spaced data.
 */
class VolumeToVtkStructuredPoints : 
  /// Inherit from VTK's structured points source.
  public vtkStructuredPointsSource 
{
public:
  /// Create new instance.
  static VolumeToVtkStructuredPoints *New();

  vtkTypeMacro(VolumeToVtkStructuredPoints,vtkStructuredPointsSource);

  /// Set associated volume object.
  void SetVolume( UniformVolume::SmartPtr& );

  /// Set affine transformation associated with volume object.
  void SetXform( AffineXform::SmartPtr& );

protected:
  /// Default constructor.
  VolumeToVtkStructuredPoints();

  /// Virtual destructor.
  virtual ~VolumeToVtkStructuredPoints() {};

  VolumeToVtkStructuredPoints(const VolumeToVtkStructuredPoints&);
  void operator=(const VolumeToVtkStructuredPoints&);

  /** Perform actual encapsulation.
   * The volume data is not copied during encapsulation. Instead, vtk is set
   * to share the pointers to the existing data. Consequently, vtk is asked
   * NOT to free the data when its objects are deleted.
   */
  void Execute();

private:
  /// Pointer to the input volume object.
  UniformVolume::SmartPtr m_Volume;

  /// Pointer to the associated transformation object.
  AffineXform::SmartPtr m_AffineXform;
};

//@}

} // namespace cmtk

#endif
