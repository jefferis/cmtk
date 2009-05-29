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

#include <cmtkAffineRegistration.h>

#include <cmtkVector.h>

#include <cmtkXform.h>
#include <cmtkAffineXform.h>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>
#include <cmtkFunctional.h>

#include <cmtkVoxelMatchingAffineFunctional.h>

#include <cmtkOptimizer.h>
#include <cmtkBestNeighbourOptimizer.h>

#include <cmtkTimers.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

AffineRegistration::AffineRegistration () 
{ 
  InitialAlignCenters = false;
  HistogramEqualization1 = HistogramEqualization2 = 0;
  NoSwitch = false;
}

AffineRegistration::~AffineRegistration () 
{
}

CallbackResult 
AffineRegistration::InitRegistration ()
{
  CallbackResult result = this->Superclass::InitRegistration();
  if ( result != CALLBACK_OK ) return result;

  UniformVolume::SmartPtr refVolume;
  UniformVolume::SmartPtr fltVolume;

  if ( NoSwitch || (Volume_1->AverageVoxelVolume() >= Volume_2->AverageVoxelVolume()) ) 
    {
    refVolume = Volume_1;
    fltVolume = Volume_2;
    SwitchVolumes = false;
    } 
  else
    {
    refVolume = Volume_2;
    fltVolume = Volume_1;
    SwitchVolumes = true;
    }
  
  AffineXform::SmartPtr affineXform;

  if ( InitialXform )
    {
    if ( SwitchVolumes ^ InitialXformIsInverse )
      {
      affineXform = AffineXform::SmartPtr( InitialXform->MakeInverse() );
      } 
    else
      {
      affineXform = AffineXform::SmartPtr( InitialXform );
      }
  } 
  else 
    {
    affineXform = AffineXform::SmartPtr( new AffineXform );
    }

  if ( InitialAlignCenters ) 
    {
    Vector3D deltaCenter = ( refVolume->GetCenterCropRegion() - fltVolume->GetCenterCropRegion() );
    affineXform->SetXlate( deltaCenter.XYZ );
    }
  
  // explicit cast needed for MIPSpro comp.
  Xform = Xform::SmartPtr::DynamicCastFrom( affineXform );
  
  Vector3D center = fltVolume->GetCenterCropRegion();
  affineXform->ChangeCenter( center.XYZ );

  if ( UseOriginalData ) 
    {  
    VoxelMatchingAffineFunctional *newFunctional = CreateAffineFunctional( Metric, refVolume, fltVolume, affineXform );
    FunctionalStack.push( Functional::SmartPtr( newFunctional ) );
    }
  
  Types::Coordinate currSampling = std::max( Sampling, 2 * std::min( refVolume->GetMinDelta(), fltVolume->GetMinDelta()));
  CallbackResult irq = CALLBACK_OK;
  
  double coarsest = CoarsestResolution;
  if ( coarsest <= 0 ) coarsest = Exploration;
  
  for ( ; ( irq == CALLBACK_OK ) && (currSampling<=coarsest); currSampling *= 2 ) 
    {
    UniformVolume::SmartPtr nextRef( NULL );
    UniformVolume::SmartPtr nextMod( NULL );
    try 
      {
      this->ReportProgress( "resampling", 0 );
      nextRef = UniformVolume::SmartPtr( new UniformVolume( *refVolume, currSampling ) );
      irq = this->ReportProgress(  "resampling", 50 );
      nextMod = UniformVolume::SmartPtr( new UniformVolume( *fltVolume, currSampling ) );
      }
    catch (...) 
      {
      }
    
    UniformVolume::SmartPtr useRef( nextRef );
    UniformVolume::SmartPtr useFlt( nextMod );
    if ( HistogramEqualization1 ) 
      {
      useRef = UniformVolume::SmartPtr( useRef->Clone() );
      useRef->GetData()->HistogramEqualization();
      }
    if ( HistogramEqualization2 )
      {
      useFlt = UniformVolume::SmartPtr( useFlt->Clone() );
      useFlt->GetData()->HistogramEqualization();
      }

    VoxelMatchingAffineFunctional *newFunctional = CreateAffineFunctional( Metric, useRef, useFlt, affineXform );
    FunctionalStack.push( Functional::SmartPtr( newFunctional ) );
    
    refVolume = nextRef;
    fltVolume = nextMod;
    }

  Optimizer = Optimizer::SmartPtr( new BestNeighbourOptimizer( OptimizerStepFactor ) );   
  Optimizer->SetCallback( Callback );
  
  // default to rigid transformation
  if ( NumberDOFs.empty() )
    NumberDOFs.push_back( 6 );
  
  // push guard elements
  NumberDOFs.push_back( -1 );
  NumberDOFsFinal.push_back( -1 );
  // intialize iterator.
  NumberDOFsIterator = NumberDOFs.begin();

  return irq;
}

void
AffineRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, const int level, const int total ) 
{
  if ( *NumberDOFsIterator < 0 )
    {
    if ( (level == total) && (NumberDOFsFinal.size()>1) )
      NumberDOFsIterator = NumberDOFsFinal.begin();
    else
      NumberDOFsIterator = NumberDOFs.begin();
    }

  AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( Xform );
  if ( affineXform ) 
    {
    int numberDOFs = *NumberDOFsIterator;
    affineXform->SetNumberDOFs( numberDOFs );
    if ( Callback ) 
      {
      char buffer[64];
      snprintf( buffer, sizeof( buffer ), "Setting Number DOFs to %d.", numberDOFs );
      Callback->Comment( buffer );
      }
    }
  this->Superclass::EnterResolution( v, f, level, total );
}

int 
AffineRegistration::DoneResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f,
  const int level, const int total )
{
  this->Superclass::DoneResolution( v, f, level, total );
  
  NumberDOFsIterator++;
  return (*NumberDOFsIterator < 0);
}

AffineXform::SmartPtr
AffineRegistration::GetTransformation() const
{
  AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( Xform );
  if ( affineXform && SwitchVolumes ) 
    {
    return affineXform->GetInverse();
    } 
  else 
    {
    return affineXform;
    }
}

} // namespace cmtk
