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

#include <cmtkVtkIsoSurface.h>

#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

VtkIsoSurface::VtkIsoSurface() :
  m_Study( NULL )
{
  Wrapper = NULL;
  ContourFilter = NULL;
  Smooth = NULL;
  WindowedSincFilter = NULL;
  Decimate = NULL;
  Mapper = NULL;
  LookupTable = NULL;
  Actor = NULL;
  ToggleDecimate = false;
  ToggleSmooth = false;
}

VtkIsoSurface::~VtkIsoSurface()
{
  if ( Actor ) Actor->Delete();
  if ( Mapper ) Mapper->Delete();
  if ( ContourFilter ) ContourFilter->Delete();
  if ( LookupTable ) LookupTable->Delete();
  if ( Smooth ) Smooth->Delete();
  if ( Decimate ) Decimate->Delete();
  if ( Wrapper ) Wrapper->Delete();
}

void
VtkIsoSurface::SetStudy( Study::SmartPtr& study )
{
  m_Study = study;
  UniformVolume::SmartPtr uniformVolume = m_Study->GetVolume();
  if ( ! uniformVolume )
    return;

  if ( ! Wrapper ) Wrapper = VolumeToVtkStructuredPoints::New();
  Wrapper->SetVolume( uniformVolume );
}

void
VtkIsoSurface::SetIsoLevel( const Types::DataItem isoLevel )
{
  if ( ! ContourFilter ) ContourFilter = vtkContourFilter::New();
  ContourFilter->SetNumberOfContours( 1 );
  ContourFilter->SetValue( 0, isoLevel );
  ContourFilter->ComputeNormalsOff();
  ContourFilter->ComputeGradientsOff();
  ContourFilter->ComputeScalarsOff();
}

void
VtkIsoSurface::SetDecimate( float targetReduction, float maximumError )
{
  if ( ! Decimate ) Decimate = vtkDecimatePro::New();
  Decimate->SetTargetReduction( targetReduction );
  Decimate->SetMaximumError( maximumError );

  ToggleDecimate = true;
}

void
VtkIsoSurface::SetSmooth()
{
  if ( ! Smooth ) Smooth = vtkSmoothPolyDataFilter::New();

  ToggleSmooth = true;
}

void
VtkIsoSurface::UpdateConnections()
{
  if ( ! Wrapper ) return;
  if ( ! ContourFilter ) return;
  if ( ! Actor ) return;
  if ( ToggleDecimate && ! Decimate ) return;
  if ( ToggleSmooth && ! Smooth ) return;

  ContourFilter->SetInput( dynamic_cast<vtkDataSet*>( Wrapper->GetOutput() ) );

  if ( ! Mapper ) Mapper = vtkPolyDataMapper::New();
  if ( ToggleMapScalars ) 
    {
    ContourFilter->ComputeScalarsOn();
    Mapper->ScalarVisibilityOn();
    Mapper->SetColorModeToMapScalars();
    
    if ( ! LookupTable ) 
      {
      LookupTable = new VtkLookupTable();
      }
    LookupTable->SetFromStudy( m_Study );
    Mapper->SetLookupTable( LookupTable );
    } 
  else
    {
    // do not reset contour filter's setting, since that would waste time!
    Mapper->ScalarVisibilityOff();
    Mapper->SetColorModeToDefault();
    }
  
  if ( ToggleDecimate ) 
    {
    if ( ToggleSmooth ) 
      {
      Smooth->SetInput( ContourFilter->GetOutput() );
      Decimate->SetInput( Smooth->GetOutput() );
      } 
    else
      {
      Decimate->SetInput( ContourFilter->GetOutput() );
      }
    Mapper->SetInput( Decimate->GetOutput() );
    } 
  else
    {
    if ( ToggleSmooth ) 
      {
      Smooth->SetInput( ContourFilter->GetOutput() );
      Mapper->SetInput( Smooth->GetOutput() );
      } 
    else
      {
      Mapper->SetInput( ContourFilter->GetOutput() );
      }
    }
  Actor->SetMapper( Mapper );
}

} // namespace cmtk
