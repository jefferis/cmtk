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

#include <cmtkFusionIsolines.h>
#include <cmtkRGB.h>

#include <limits.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

FusionIsolines::FusionIsolines() 
{ 
  IsolineImageIndex = 0; 

  this->m_IsolineFilter = IsolineFilter::New();

  ColormapImage[0] = Colormap::New();
  ColormapImage[1] = Colormap::New();

  LookupImage = ImageToImageRGB::New();
  LookupImage->SetColormap( ColormapImage[1] );

  LookupIsolines = ImageToImageRGB::New();
  LookupIsolines->SetColormap( ColormapImage[1] );
  LookupIsolines->SetInput( this->m_IsolineFilter->GetOutput() );

  this->m_FusionOverlay = FusionOverlay::New();
  this->m_FusionOverlay->SetInput( 0, LookupImage->GetOutput() );
  this->m_FusionOverlay->SetInput( 1, LookupIsolines->GetOutput() );
  this->m_FusionOverlay->SetTopImageIndex( 1 );
  this->m_FusionOverlay->SetLowerThreshold( 1 );
  this->m_FusionOverlay->SetUpperThreshold( FLT_MAX );

  Output = this->m_FusionOverlay->GetOutput();
  Output->Reference();
};

FusionIsolines::~FusionIsolines()
{
  this->m_FusionOverlay->Delete();
  LookupIsolines->Delete();
  LookupImage->Delete();
  ColormapImage[1]->Delete();
  ColormapImage[0]->Delete();
  this->m_IsolineFilter->Delete();
  Output->Delete();
}

void
FusionIsolines::SetIsolineImageIndex( const int isolineImageIndex )
{
  IsolineImageIndex = isolineImageIndex;
  
  this->m_IsolineFilter->SetInput( Input[IsolineImageIndex] );
  this->m_FusionOverlay->SetOriginalImage( 0, Input[1-IsolineImageIndex] );
  this->m_FusionOverlay->SetOriginalImage( 1, this->m_IsolineFilter->GetOutput() );
  LookupImage->SetInput( Input[1-IsolineImageIndex] );
  LookupImage->SetColormap( ColormapImage[1-IsolineImageIndex] );
  LookupIsolines->SetColormap( ColormapImage[IsolineImageIndex] );
}

void
FusionIsolines::Execute()
{
  this->m_FusionOverlay->Update();
  this->UpdateExecuteTime();
}

} // namespace cmtk
