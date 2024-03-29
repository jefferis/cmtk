/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkAtlasSegmentation.h"

#include <Registration/cmtkRegistrationCallback.h>
#include <Registration/cmtkAffineRegistration.h>
#include <Registration/cmtkElasticRegistration.h>
#include <Registration/cmtkReformatVolume.h>

#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>

cmtk::AtlasSegmentation::AtlasSegmentation
( UniformVolume::SmartPtr& targetImage, UniformVolume::SmartPtr& atlasImage, UniformVolume::SmartPtr& atlasLabels ) :
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
  ar.SetAccuracy( .1 * this->m_TargetImage->GetMaxDelta() );
  ar.SetSampling( 2 * this->m_TargetImage->GetMaxDelta() );

  ar.SetUseOriginalData( !this->m_Fast );
  
  (DebugOutput( 1 ) << "Affine registration...").flush();
  ar.Register();
  DebugOutput( 1 ) << " done.\n";

  try
    {
    this->m_AffineXform = ar.GetTransformation();
    }
  catch ( const AffineXform::MatrixType::SingularMatrixException& )
    {
    StdErr << "ERROR: singular matrix in AffineRegistration::GetTransformation()\n";
    throw ExitException( 1 );
    }
}

void
cmtk::AtlasSegmentation
::RegisterSpline()
{
  ElasticRegistration er;
  er.SetVolume_1( this->m_TargetImage );
  er.SetVolume_2( this->m_AtlasImage );
  er.SetInitialTransformation( this->GetAffineXform() );
  
  er.SetUseOriginalData( !this->m_Fast );
  er.SetFastMode( this->m_Fast );

  const Types::Coordinate minSize = std::min( std::min( this->m_TargetImage->m_Size[0], this->m_TargetImage->m_Size[1] ), this->m_TargetImage->m_Size[2] );
  er.SetGridSpacing( minSize / 2 );
  er.SetRefineGrid( std::max<int>( 0, static_cast<int>( (log( minSize / this->m_TargetImage->GetMaxDelta() ) / log(2.0)) - 3 ) ) );
  er.SetDelayRefineGrid( !this->m_Fast );
  
  er.SetGridEnergyWeight( 1e-1f );
  er.SetAdaptiveFixParameters( true );

  er.SetAlgorithm( 3 );
  er.SetExploration( minSize / 8 );
  er.SetAccuracy( .1 * this->m_TargetImage->GetMinDelta() );
  er.SetSampling( 2 * this->m_TargetImage->GetMaxDelta() );  
  
  (DebugOutput( 1 ) << "Nonrigid registration...").flush();
  er.Register();
  DebugOutput( 1 ) << " done.\n";

  this->m_WarpXform = er.GetTransformation();
}

void
cmtk::AtlasSegmentation
::ReformatLabels()
{
  ReformatVolume reformat;
  reformat.SetInterpolation( Interpolators::PARTIALVOLUME );
  reformat.SetPaddingValue( 0 );
  reformat.SetUsePaddingValue( true );
  reformat.SetReferenceVolume( this->m_TargetImage );
  reformat.SetFloatingVolume( this->m_AtlasLabels );

  WarpXform::SmartPtr warpXform = this->GetWarpXform();
  reformat.SetWarpXform( warpXform );

  this->m_LabelMap = UniformVolume::SmartPtr( reformat.PlainReformat() );
}
