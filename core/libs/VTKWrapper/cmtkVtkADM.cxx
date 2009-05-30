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

#include <cmtkVtkADM.h>

#include <cmtkVolumeToVtkStructuredPoints.h>

#include <vtkDataSet.h>
#include <vtkStructuredPoints.h>
#include <vtkDiscreteMarchingCubes.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkPolyDataNormals.h>

#include <vtkLookupTable.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

void
VtkADM::CreatePipeline()
{
  this->m_Volume->GetData()->SetDataClass( DATACLASS_LABEL );
  UniformVolume::SmartPtr uniformVolume( new UniformVolume( *this->m_Volume, this->m_Volume->GetMinDelta() * 4 ) );
  if ( ! uniformVolume )
    return;

  VolumeToVtkStructuredPoints *wrapper = VolumeToVtkStructuredPoints::New();
  wrapper->SetVolume( uniformVolume );

  Types::DataItem min, max;
  uniformVolume->GetData()->GetRange( min, max );
  int maxContour = static_cast<int>( max );

  vtkDiscreteMarchingCubes *marchingCubes = vtkDiscreteMarchingCubes::New();
  marchingCubes->SetInput( wrapper->GetOutput() );
  marchingCubes->GenerateValues( maxContour, 1, maxContour );
  marchingCubes->ComputeScalarsOn();
  marchingCubes->ComputeNormalsOff();
  marchingCubes->ComputeGradientsOff();

  WarpPolyData = WarpVtkPolyData::New();

  WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( ADM );
  ADM->Compose( NULL );
  WarpPolyData->SetWarpXform( warpXform );

  vtkSmoothPolyDataFilter* Smooth = vtkSmoothPolyDataFilter::New();
  Smooth->SetNumberOfIterations( 100 );
  Smooth->SetInput( marchingCubes->GetOutput() );

  vtkDecimatePro* Decimate = vtkDecimatePro::New();
  Decimate->SetInput( Smooth->GetOutput() );
  Decimate->SetTargetReduction( 0.9 );
  Decimate->SetPreserveTopology( true );
  Decimate->SetMaximumError( uniformVolume->GetMinDelta() / 2 );

  vtkPolyDataNormals* polyDataNormals = vtkPolyDataNormals::New();
  polyDataNormals->SetInput( Decimate->GetOutput() );
  polyDataNormals->SplittingOff();
  polyDataNormals->SetFeatureAngle( 15 ); // according to Bill Lorensen post

  WarpPolyData->SetInput( polyDataNormals->GetOutput() );

  vtkLookupTable* LookupTable = vtkLookupTable::New();
  LookupTable->SetNumberOfTableValues( maxContour + 1 );
  LookupTable->SetTableRange( 0, max );
  LookupTable->SetHueRange( 0.66, 0 );
  LookupTable->SetSaturationRange( 1, 1 );
  LookupTable->SetValueRange( 1, 1 );
  LookupTable->Build();

  VolMapper = vtkPolyDataMapper::New();
  VolMapper->SetInput( WarpPolyData->GetOutput() );
  VolMapper->ScalarVisibilityOn();
  VolMapper->SetLookupTable( LookupTable );
  VolMapper->SetColorModeToMapScalars();
  VolMapper->SetScalarRange( 0, max );
}

void
VtkADM::SetModeWeights( const Types::Coordinate* weights )
{
  ADM->Compose( weights );
  if ( WarpPolyData ) 
    WarpPolyData->Execute();
}

} // namespace cmtk
