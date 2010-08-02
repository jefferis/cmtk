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
#include "System/cmtkCommandLine.h"

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
    m_OutFile( NULL )
{
#ifdef CMTK_USE_SQLITE
  this->m_UpdateDB = NULL;
#endif
}

int
cmtk::SimpleLevelsetCommandLineBase::Init( const int argc, const char* argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv, cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Levelset segmentation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Levelset-type segmentation of foreground/background using minimum regional variance energy" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Segmentation" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &this->m_Verbose, true, "Verbose mode" )->SetProperties( cmtk::CommandLine::PROPS_NOXML );

    cl.AddSwitch( Key( 'b', "binarize" ), &this->m_Binarize, true, "Binarize levelset and write as byte mask, rather than write floating-point levelset function itself." );

    cl.BeginGroup( "Levelset Evolution Parameters", "These parameters of control the evolution of the levelset function" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( 'n', "iterations" ), &this->m_NumberOfIterations, "Maximum number of iterations" );
    cl.AddSwitch( Key( 'f', "force-iterations" ), &this->m_ForceIterations, true, "Force given number of iterations, even when convergence has been detected" );

    cl.AddOption( Key( 's', "filter-sigma" ), &this->m_FilterSigma, "Gaussian filter sigma in world coordinate units (e.g., mm)" );
    cl.AddOption( Key( 'd', "delta" ), &this->m_TimeDelta, "Time constant for levelset evolution; must be > 0; larger is faster" );
    cl.AddOption( Key( 't', "levelset-threshold" ), &this->m_LevelsetThreshold, "Levelset threshold: levelset function is truncated at +/- this value" );
    cl.EndGroup();

#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &this->m_UpdateDB, "Path to image/transformation database that should be updated with the newly created image." );
    cl.EndGroup();
#endif

    cl.AddParameter( &this->m_InFile, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &this->m_OutFile, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_LABELS | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse();
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
