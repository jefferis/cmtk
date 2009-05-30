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

#ifndef __cmtkVtkIsoSurface_h_included_
#define __cmtkVtkIsoSurface_h_included_

#include <cmtkconfig.h>

#include <cmtkVtkObject.h>

#include <vtkDataSet.h>
#include <vtkStructuredPoints.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkWindowedSincPolyDataFilter.h>

#include <vtkActor.h>

#include <cmtkStudy.h>

#include <cmtkVolumeToVtkStructuredPoints.h>
#include <cmtkWarpVtkPolyData.h>
#include <cmtkVtkLookupTable.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

/// This class handles isosurface objects within a VTK scene.
class VtkIsoSurface :
  /// This is a VTK graphics object.
  public VtkObject
{
public:
  /// Constructor.
  VtkIsoSurface();

  /// Destructor.
  virtual ~VtkIsoSurface();

  /// Set volume.
  void SetStudy( Study::SmartPtr& study );

  /// Set isolevel.
  void SetIsoLevel( const Types::DataItem isoLevel );

  /// Set mapping of scalars.
  void SetMapScalars( const bool mapScalars ) 
  { ToggleMapScalars = mapScalars; }

  /// Set decimation.
  void SetDecimate( float targetReduction, float maximumError );

  /// Disable decimation.
  void DisableDecimate() { ToggleDecimate = false; }

  /// Set smoothing.
  void SetSmooth();

  /// Disable smoothing.
  void DisableSmooth() { ToggleSmooth = false; }

  /// Return VTK actor for this object.
  vtkActor* GetVtkActor() {
    if ( !Actor ) Actor = vtkActor::New();
    this->UpdateConnections(); 
    return Actor; 
  }

private:
  /// The study that this object was created from.
  Study::SmartPtr m_Study;

  /// Wrapper for Volume into VTK representation.
  VolumeToVtkStructuredPoints* Wrapper;

  /// Iso-Contour filter.
  vtkContourFilter* ContourFilter;

  /// Grid smoothing filter.
  vtkSmoothPolyDataFilter* Smooth;

  /// Toggle for decimation.
  bool ToggleSmooth;

  /// Windowed sinc smoothing filter.
  vtkWindowedSincPolyDataFilter* WindowedSincFilter; 

  /// Grid decimation filter.
  vtkDecimatePro* Decimate;

  /// Toggle for decimation.
  bool ToggleDecimate;

  /// VTK mapper.
  vtkPolyDataMapper* Mapper;

  /// Color lookup table.
  VtkLookupTable* LookupTable;

  /// Scalar mapping toggle.
  bool ToggleMapScalars;

  /// Actor to be added to VTK renderer.
  vtkActor* Actor;

  /// Update all internal connections.
  void UpdateConnections();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVtkIsoSurface_h_included_
