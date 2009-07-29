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
#include <cmtkReformatVolume.cxx>

cmtk::AtlasSegmentation::AtlasSegmentation
( UniformVolume::SmartPtr& targetImage, UniformVolume::SmartPtr& atlasImage, UniformVolume::SmartPtr& atlasLabels )  
  : m_Verbose( false ),
    m_Fast( false ),
    m_TargetImage( targetImage ),
    m_AtlasImage( atlasImage ),
    m_AtlasLabels( atlasLabels ),
    m_LabelMap( NULL )
{
}

void
cmtk::AtlasSegmentation
::RegisterAffine()
{
  AffineRegistration ar;
  ar.SetVolume_1( this->m_TargetImage );
  ar.SetVolume_2( this->m_AtlasImage );
  
  // run 6 DOFs first, then 9 at each level.
  ar.AddNumberDOFs( 6 );
  ar.AddNumberDOFs( 9 );
  
  ar.SetInitialAlignCenters( true );
  ar.SetExploration( 4 * this->m_TargetImage->GetMaxDelta() );
  ar.SetAccuracy( .1 * this->m_TargetImage->GetMinDelta() );
  ar.SetSampling( 2 * this->m_TargetImage->GetMaxDelta() );

  ar.SetUseOriginalData( !this->m_Fast );
  
  if ( this->m_Verbose ) 
    {
    StdErr << "Affine registration...";
    StdErr.flush();
    }
  ar.Register();
  if ( this->m_Verbose ) 
    {
    StdErr << " done.\n";
    }
  
  this->m_AffineXform = ar.GetTransformation();
}

void
cmtk::AtlasSegmentation
::RegisterSpline()
{
  ElasticRegistration er;
  er.SetVolume_1( this->m_TargetImage );
  er.SetVolume_2( this->m_AtlasImage );
  er.SetInitialXform( this->GetAffineXform() );
  
  er.SetAlgorithm( 3 );
  er.SetExploration( 8 * this->m_TargetImage->GetMaxDelta() );
  er.SetAccuracy( .1 * this->m_TargetImage->GetMinDelta() );
  er.SetSampling( 2 * this->m_TargetImage->GetMaxDelta() );
  
  er.SetUseOriginalData( !this->m_Fast );
  er.SetFastMode( this->m_Fast );

  const Types::Coordinate minSize = std::min( std::min( this->m_TargetImage->Size[0], this->m_TargetImage->Size[1] ), this->m_TargetImage->Size[2] );
  er.SetGridSpacing( minSize / 2 );
  er.SetRefineGrid( std::max( 0, static_cast<int>( (log( minSize / this->m_TargetImage->GetMinDelta() ) / log(2)) - 2 ) ) );
  er.SetDelayRefineGrid( !this->m_Fast );
  
  er.SetGridEnergyWeight( 1e-1 );
  
  er.SetAdaptiveFixParameters( true );
  
  if ( this->m_Verbose ) 
    {
    StdErr << "Nonrigid registration...";
    StdErr.flush();
    }
  er.Register();
  if ( this->m_Verbose ) 
    {
    StdErr << " done.\n";
    }
  this->m_WarpXform = er.GetTransformation();
}

void
cmtk::AtlasSegmentation
::ReformatLabels()
{
  ReformatVolume reformat;
  reformat.SetInterpolation( Interpolators::PARTIALVOLUME );
  reformat.SetPaddingValue( 0 );
  reformat.SetReferenceVolume( this->m_TargetImage );
  reformat.SetFloatingVolume( this->m_AtlasLabels );

  WarpXform::SmartPtr warpXform = this->GetWarpXform();
  reformat.SetWarpXform( warpXform );

  this->m_LabelMap = UniformVolume::SmartPtr( reformat.PlainReformat() );
}
