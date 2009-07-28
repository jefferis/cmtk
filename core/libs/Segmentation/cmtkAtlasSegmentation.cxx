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

#include <cmtkAtlasSegmentation.h>

#include <cmtkRegistrationCallback.h>
#include <cmtkAffineRegistration.h>
#include <cmtkElasticRegistration.h>

cmtk::AtlasSegmentation::AtlasSegmentation
( UniformVolume::SmartPtr& targetImage, UniformVolume::SmartPtr& atlasImg, UniformVolume::SmartPtr& atlasLbl, const bool verbose )
  : m_LabelMap( NULL )
{
  AffineXform::SmartPtr affine;
  {
  AffineRegistration ar;
  ar.SetVolume_1( targetImage );
  ar.SetVolume_2( atlasImg );
  
  // run 6 DOFs first, then 9 at each level.
  ar.AddNumberDOFs( 6 );
  ar.AddNumberDOFs( 9 );
  
  ar.SetInitialAlignCenters( true );
  ar.SetExploration( 4 * targetImage->GetMaxDelta() );
  ar.SetAccuracy( .1 * targetImage->GetMinDelta() );
  ar.SetSampling( 2 * targetImage->GetMaxDelta() );
  
  if ( verbose ) 
    {
    StdErr << "Affine registration...";
    StdErr.flush();
    }
  ar.Register();
  if ( verbose ) 
    {
    StdErr << " done.\n";
    }
  
  affine = ar.GetTransformation();
  }
  
  WarpXform::SmartPtr xform;
  {
  ElasticRegistration er;
  er.SetVolume_1( targetImage );
  er.SetVolume_2( atlasImg );
  er.SetInitialXform( affine );
  
  er.SetAlgorithm( 3 );
  er.SetExploration( 4 * targetImage->GetMaxDelta() );
  er.SetAccuracy( .1 * targetImage->GetMinDelta() );
  er.SetSampling( 2 * targetImage->GetMaxDelta() );
  
  er.SetGridSpacing( targetImage->Size[0] / 4 );
  er.SetRefineGrid( 3 );
  //  er.SetUseOriginalData( false );
  er.SetGridEnergyWeight( 1e-1 );
  
  er.SetFastMode( true );
  er.SetAdaptiveFixParameters( true );
  
  if ( verbose ) 
    {
    StdErr << "Nonrigid registration...";
    StdErr.flush();
    }
  er.Register();
  if ( verbose ) 
    {
    StdErr << " done.\n";
    }
  xform = WarpXform::SmartPtr::DynamicCastFrom( er.m_Xform );
  }
  
}

