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

#ifndef __cmtkVtkADM_h_included_
#define __cmtkVtkADM_h_included_

#include <cmtkconfig.h>

#include <cmtkMacros.h>
#include <cmtkVolume.h>
#include <cmtkActiveDeformationModel.h>

#include <cmtkWarpVtkPolyData.h>

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

/// VTK representation of an Active Deformation Model.
class VtkADM
{
  /// Resolution before iso-surface computation.
  igsGetSetMacro(Types::Coordinate,Resolution);

  /// Iso-level.
  igsGetSetMacro(Types::DataItem,Level);

public:
  /// Constructor.
  VtkADM( SplineActiveDeformationModel::SmartPtr& adm, const UniformVolume::SmartPtr& volume ) 
  {
    this->m_Volume = volume; 
    ADM = adm; 
    VolMapper = NULL;
  }

  /// Return VTK Actor.
  vtkActor* GetActor() 
  { 
    if ( ! VolMapper ) 
      this->CreatePipeline();
    Actor = vtkActor::New();
    Actor->SetMapper( VolMapper );
    Actor->GetProperty()->SetInterpolationToGouraud();
    Actor->GetProperty()->SetDiffuse( 0.3 );
    Actor->GetProperty()->SetSpecular( 0.3 );
    return Actor; 
  }

  /// Set ADM mode weights.
  void SetModeWeights( const Types::Coordinate* weights );

private:
  /// The image to visualize.
  UniformVolume::SmartPtr m_Volume;

  /// The active deformation model.
  SplineActiveDeformationModel::SmartPtr ADM;

  /// VTK actor to render.
  vtkActor* Actor;

  /// VTK volume mapper.
  vtkPolyDataMapper* VolMapper;

  /// The poly data warp filter.
  WarpVtkPolyData* WarpPolyData;

  /// Create VTK pipeline.
  void CreatePipeline();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVtkADM_h_included_
