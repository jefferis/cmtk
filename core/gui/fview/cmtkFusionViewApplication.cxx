/*
//
//  Copyright 2010 SRI International
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

#include "cmtkFusionViewApplication.h"

#include <System/cmtkCommandLine.h>

#include <IO/cmtkVolumeIO.h>

#include <QtGui/QActionGroup>

cmtk::FusionViewApplication
::FusionViewApplication( int argc, char* argv[] ) 
  : QApplication( argc, argv ),
    m_MainWindow( new QMainWindow ),
    m_SliceAxis( 2 )
{
  CommandLine cl;
  cl.SetProgramInfo( CommandLine::PRG_TITLE, "Fusion viewer." );

  const char* imagePathFix;
  cl.AddParameter( &imagePathFix, "FixedImage", "Fixed image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

  const char* imagePathMov;
  cl.AddParameter( &imagePathMov, "MovingImage", "Moving image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

  try
    {
    cl.Parse( argc, const_cast<const char**>( argv ) );
    }
  catch ( const CommandLine::Exception& ex )
    {
    throw(ex);
    }

  this->m_FixedVolume = VolumeIO::Read( imagePathFix );
  if ( ! this->m_FixedVolume )
    {
    exit( 1 );
    }

  this->m_MovingVolume = VolumeIO::Read( imagePathMov );
  if ( ! this->m_MovingVolume )
    {
    exit( 1 );
    }

  this->m_MainWindowUI.setupUi( this->m_MainWindow );

  this->m_MainWindowUI.alphaSlider->setRange( 0, 1000 );
  this->m_MainWindowUI.sliceSlider->setRange( 0, this->m_FixedVolume->GetDims()[this->m_SliceAxis] );

  const Types::DataItemRange rangeFix = this->m_FixedVolume->GetData()->GetRange();
  this->m_MainWindowUI.blackSliderFix->setRange( rangeFix.m_LowerBound, rangeFix.m_UpperBound );
  this->m_MainWindowUI.whiteSliderFix->setRange( rangeFix.m_LowerBound, rangeFix.m_UpperBound );
  
  const Types::DataItemRange rangeMov = this->m_MovingVolume->GetData()->GetRange();
  this->m_MainWindowUI.blackSliderMov->setRange( rangeMov.m_LowerBound, rangeMov.m_UpperBound );
  this->m_MainWindowUI.whiteSliderMov->setRange( rangeMov.m_LowerBound, rangeMov.m_UpperBound );
  
  QActionGroup* zoomGroup = new QActionGroup( this->m_MainWindow );
  zoomGroup->setExclusive( true );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom25 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom50 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom100 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom200 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom300 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom400 );
  
  QActionGroup* sliceGroup = new QActionGroup( this->m_MainWindow );
  sliceGroup->setExclusive( true );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceAxial_XY );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceCoronal_XZ );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceSagittal_YZ );
  
  this->m_MainWindow->show();
}

