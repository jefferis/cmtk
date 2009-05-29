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

#include <cmtkFusionSlicers.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

SlicerPipeline::SlicerPipeline
( Plane *const slicePlane, Study::SmartPtr& study, 
  AffineXform::SmartPtr& affineXform, WarpXform::SmartPtr& warpXform )
{
  this->m_VolumeWrapper = VolumeWrapper::New();
  UniformVolume::SmartPtr volume = study->GetVolume();
  this->m_VolumeWrapper->SetVolume( volume );

  this->m_VolumeWrapper->SetAffineXform( affineXform );
  this->m_VolumeWrapper->SetWarpXform( warpXform );

  this->m_Slicer = Slicer::New();
  this->m_Slicer->SetInput( this->m_VolumeWrapper );
  this->m_Slicer->SetPlane( slicePlane );
}

SlicerPipeline::~SlicerPipeline()
{
  if ( this->m_Slicer ) 
    this->m_Slicer->Delete();
  if ( this->m_VolumeWrapper ) 
    this->m_VolumeWrapper->Delete();
}

Image*
SlicerPipeline::GetOutput()
{
  return this->m_Slicer->GetOutput();
}

FusionSlicers::FusionSlicers() :
  ApplyWarp( false ),
  ReferenceStudy( NULL )
{
  InterpolationMode = cmtk::Interpolators::LINEAR;
  SliceNormal = AXIS_Z;
  ReferenceSlice = Image::New();
}

FusionSlicers::~FusionSlicers()
{
  if ( ReferenceSlice ) 
    ReferenceSlice->Delete();
}

void
FusionSlicers::SetReferenceStudy( Study::SmartPtr& referenceStudy )
{
  if ( ReferenceStudy != referenceStudy ) 
    {
    ReferenceStudy = referenceStudy;
    
    if ( ReferenceStudy ) 
      {
      UniformVolume::SmartPtr volume = ReferenceStudy->GetVolume();
      if ( volume ) 
	{
	ReferenceSlicePosition = volume->Size[SliceNormal] / 2;
	}
      this->UpdateSlicePlane();
      this->ClearAndDelete();
      }
    }
}

void
FusionSlicers::SetStudyList( StudyList::SmartPtr& studyList )
{
  this->m_StudyList = studyList;
  this->ClearAndDelete();
}

Image* 
FusionSlicers::GetOutput( Study::SmartPtr& study )
{
  if ( study == ReferenceStudy ) return ReferenceSlice;

  study->ReadVolume();

  if ( this->find( study ) == this->end() ) 
    {
    AffineXform::SmartPtr affineXform;
    WarpXform::SmartPtr warpXform;
    
    StudyToXform::iterator it = (*this->m_StudyList)[ReferenceStudy].find( study );
    while ( it != (*this->m_StudyList)[ReferenceStudy].end() ) 
      {
      if ( AffineXform::SmartPtr::DynamicCastFrom(it->second) )
	affineXform = AffineXform::SmartPtr::DynamicCastFrom(it->second);
      else
	if ( WarpXform::SmartPtr::DynamicCastFrom(it->second) )
	  warpXform = WarpXform::SmartPtr::DynamicCastFrom(it->second);
      ++it;
    }
    
    SlicerPipeline *pipeline = new SlicerPipeline( ReferenceSlice, study, affineXform, warpXform );
    Slicer *slicer = pipeline->Getm_Slicer();
    Image *slice = slicer->GetOutput();
    (*this)[ study ] = pipeline;
    return slice;
    } 
  else 
    {
    return (*this)[study]->Getm_Slicer()->GetOutput();
    }
}

void 
FusionSlicers::SetReslicePlane( const Plane* plane )
{
  ReferenceSlice->CopyStructure( plane );
}

void
FusionSlicers::SetInterpolationMode
( const cmtk::Interpolators::InterpolationEnum interpolationMode )
{
  InterpolationMode = interpolationMode;
  iterator it = this->begin();
  while ( it != this->end() ) 
    {
    if ( it->second ) 
      it->second->Getm_Slicer()->SetInterpolationMode( InterpolationMode );
    ++it;
    }  
}

void
FusionSlicers::SetApplyWarp
( const bool applyWarp )
{
  ApplyWarp = applyWarp;
  iterator it = this->begin();
  while ( it != this->end() ) 
    {
    if ( it->second ) 
      it->second->Getm_Slicer()->SetApplyWarp( ApplyWarp );
    ++it;
    }  
}

void 
FusionSlicers::ClearAndDelete()
{
  iterator it = this->begin();
  while ( it != this->end() ) 
    {
    if ( it->second ) delete it->second;
    this->erase( it );
    it = this->begin();
    }
}

void
FusionSlicers::UpdateSlicePlane()
{
  if ( ! ReferenceStudy ) return;
  UniformVolume::SmartPtr volume = ReferenceStudy->GetVolume();
  if ( volume ) 
    {
    ScalarImage::SmartPtr orthoSlice( volume->GetOrthoSliceInterp( SliceNormal, ReferenceSlicePosition ) );
    
    if ( SliceNormal != AXIS_Z )
      orthoSlice->Mirror( false /* horizontal */, true /* vertical */ );
    orthoSlice->AdjustAspect();

    ReferenceSlice->SetFromScalarImage( orthoSlice.GetPtr(), false /* copyData */ );
  }
}

} // namespace cmtk
