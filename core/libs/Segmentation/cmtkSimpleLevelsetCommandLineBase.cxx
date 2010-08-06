/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkSimpleLevelsetCommandLineBase.h"

#include "System/cmtkConsole.h"

#include "IO/cmtkVolumeIO.h"

cmtk::SimpleLevelsetCommandLineBase::SimpleLevelsetCommandLineBase()
  : m_Verbose( false ),
    m_FilterSigma( 2.0 ),    
    m_TimeDelta( 0.1 ),
    m_LevelsetThreshold( 1.0 ),
    m_NumberOfIterations( 100 ),
    m_ForceIterations( false ),
    m_Binarize( false ),    
    m_InFile( NULL ),
    m_OutFile( NULL ),
    m_CommandLine( cmtk::CommandLine::PROPS_XML )
{
#ifdef CMTK_USE_SQLITE
  this->m_UpdateDB = NULL;
#endif

  this->m_CommandLine.SetProgramInfo( CommandLine::PRG_TITLE, "Levelset segmentation" );
  this->m_CommandLine.SetProgramInfo( CommandLine::PRG_DESCR, "Levelset-type segmentation of foreground/background using minimum regional variance energy" );
  this->m_CommandLine.SetProgramInfo( CommandLine::PRG_CATEG, "CMTK.Segmentation" );
  
  typedef CommandLine::Key Key;
  this->m_CommandLine.AddSwitch( Key( 'v', "verbose" ), &this->m_Verbose, true, "Verbose mode" )->SetProperties( CommandLine::PROPS_NOXML );
  
  this->m_CommandLine.AddSwitch( Key( 'b', "binarize" ), &this->m_Binarize, true, "Binarize levelset and write as byte mask, rather than write floating-point levelset function itself." );
  
  this->m_CommandLine.BeginGroup( "Levelset Evolution Parameters", "These parameters of control the evolution of the levelset function" )->SetProperties( CommandLine::PROPS_ADVANCED );
  this->m_CommandLine.AddOption( Key( 'n', "iterations" ), &this->m_NumberOfIterations, "Maximum number of iterations" );
  this->m_CommandLine.AddSwitch( Key( 'f', "force-iterations" ), &this->m_ForceIterations, true, "Force given number of iterations, even when convergence has been detected" );
  
  this->m_CommandLine.AddOption( Key( 's', "filter-sigma" ), &this->m_FilterSigma, "Gaussian filter sigma in world coordinate units (e.g., mm)" );
  this->m_CommandLine.AddOption( Key( 'd', "delta" ), &this->m_TimeDelta, "Time constant for levelset evolution; must be > 0; larger is faster" );
  this->m_CommandLine.AddOption( Key( 't', "levelset-threshold" ), &this->m_LevelsetThreshold, "Levelset threshold: levelset function is truncated at +/- this value" );
  this->m_CommandLine.EndGroup();
  
#ifdef CMTK_USE_SQLITE
  this->m_CommandLine.BeginGroup( "Database", "Image/Transformation Database" );
  this->m_CommandLine.AddOption( Key( "db" ), &this->m_UpdateDB, "Path to image/transformation database that should be updated with the newly created image." );
  this->m_CommandLine.EndGroup();
#endif
  
  this->m_CommandLine.AddParameter( &this->m_InFile, "InputImage", "Input image path" )->SetProperties( CommandLine::PROPS_IMAGE );
  this->m_CommandLine.AddParameter( &this->m_OutFile, "OutputImage", "Output image path" )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_LABELS | CommandLine::PROPS_OUTPUT );
}

int
cmtk::SimpleLevelsetCommandLineBase::Init( const int argc, const char* argv[] )
{
  try
    {
    this->m_CommandLine.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex;
    return 1;
    }
  
  this->m_Volume = cmtk::VolumeIO::ReadOriented( this->m_InFile, this->m_Verbose );

  if ( !this->m_Volume )
    return 1;

  return 0;
}
