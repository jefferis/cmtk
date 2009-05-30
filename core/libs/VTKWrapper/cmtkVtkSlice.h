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

#ifndef __cmtkVtkSlice_h_included_
#define __cmtkVtkSlice_h_included_

#include <cmtkconfig.h>

#include <cmtkVtkObject.h>

#include <vtkDataSet.h>
#include <vtkStructuredPoints.h>
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
class VtkSlice :
  /// This is a VTK graphics object.
  public VtkObject
{
public:
  /// Constructor.
  VtkSlice();

  /// Destructor.
  virtual ~VtkSlice();

  /// Set volume.
  void SetStudy( Study_P& study );

  /// Return VTK actor for this object.
  vtkActor* GetVtkActor() 
  {
    if ( !Actor ) Actor = vtkActor::New();
    this->UpdateConnections(); 
    return Actor; 
  }

private:
  /// The study that this object was created from.
  Study_P Study;

  /// Actor to be added to VTK renderer.
  vtkActor* Actor;

  /// Update all internal connections.
  void UpdateConnections();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVtkSlice_h_included_
