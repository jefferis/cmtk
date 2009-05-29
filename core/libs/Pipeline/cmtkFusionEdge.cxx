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

#include <cmtkFusionEdge.h>

#include <cmtkRGB.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

FusionEdge::FusionEdge() 
{ 
  EdgeOperator = ImageEdgeOperator::New();

  ColormapImage = Colormap::New();
  ColormapEdge = Colormap::New();

  LookupImage = ImageToImageRGB::New();
  LookupImage->SetColormap( ColormapImage );
  LookupImage->SetAlphaMode( ImageToImageRGB::AlphaModeConst );

  LookupEdge = ImageToImageRGB::New();
  LookupEdge->SetColormap( ColormapEdge );
  LookupEdge->SetInput( EdgeOperator->GetOutput() );
  LookupEdge->SetAlphaMode( ImageToImageRGB::AlphaModeRamp );

  this->m_FusionAlpha = FusionAlpha::New();
  this->m_FusionAlpha->SetInput( 0, LookupImage->GetOutput() );
  this->m_FusionAlpha->SetInput( 1, LookupEdge->GetOutput() );
  this->m_FusionAlpha->SetTopImageIndex( 1 );
  this->m_FusionAlpha->SetOutputHasAlpha( true );

  Output = this->m_FusionAlpha->GetOutput();
  Output->Reference();
};

FusionEdge::~FusionEdge()
{
  this->m_FusionAlpha->Delete();
  LookupEdge->Delete();
  LookupImage->Delete();
  ColormapEdge->Delete();
  ColormapImage->Delete();
  EdgeOperator->Delete();
  Output->Delete();
}

void
FusionEdge::Execute()
{
  this->m_FusionAlpha->Update();
  this->UpdateExecuteTime();
}

void
FusionEdge::InputLinkChanged( const int inputIndex )
{
  if ( inputIndex == 0 )
    EdgeOperator->SetInput( Input[0] );
  else
    LookupImage->SetInput( Input[1] );
}

} // namespace cmtk
