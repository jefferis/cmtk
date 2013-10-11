/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#ifndef __cmtkSimpleLevelsetCommandLine_h_included_
#define __cmtkSimpleLevelsetCommandLine_h_included_

#include <cmtkconfig.h>

#include <System/cmtkProgressConsole.h>
#include <System/cmtkExitException.h>

#include <Segmentation/cmtkSimpleLevelsetCommandLineBase.h>

#include <IO/cmtkVolumeIO.h>

#ifdef CMTK_USE_SQLITE
#  include <Registration/cmtkImageXformDB.h>
#endif


namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Command line interface class template for simple levelset segmentation with a particular implementation (CPU or GPU).
 */
template<class TImpl>
class SimpleLevelsetCommandLine
  : public SimpleLevelsetCommandLineBase
{
public:
  /// This class.
  typedef SimpleLevelsetCommandLine<TImpl> Self;

  /// The actual levelset implementation.
  typedef TImpl SimpleLevelsetImplementation;

  /// Execute levelset segmentation.
  void Execute()
  {
    // Instantiate programm progress indicator.
    cmtk::ProgressConsole progressIndicator( "LevelsetSegmentation" );
    
    SimpleLevelsetImplementation levelset( this->m_Volume );
    levelset.SetScaleInitialSphere( this->m_ScaleInitialSphere );
    levelset.SetFilterSigma( cmtk::Units::GaussianSigma( this->m_FilterSigma ) );
    levelset.SetTimeDelta( this->m_TimeDelta );
    levelset.SetLevelsetThreshold( this->m_LevelsetThreshold );
    
    levelset.InitializeCenteredSphere();

    try
      {
      levelset.Evolve( this->m_NumberOfIterations, this->m_ForceIterations );
      }
    catch ( const SimpleLevelset::DegenerateLevelsetException& ex )
      {
      StdErr << "ERROR: degenerate levelset (all foreground or all background).\n";
      throw ExitException( 1 );
      }
    
    cmtk::VolumeIO::Write( *levelset.GetLevelset( this->m_Binarize ), this->m_OutFile );
    
#ifdef CMTK_USE_SQLITE
    if ( this->m_UpdateDB )
      {
      try
	{
	cmtk::ImageXformDB db( this->m_UpdateDB );
	db.AddImage( this->m_OutFile, this->m_InFile );
	}
      catch ( const cmtk::SQLite::Exception& ex )
	{
	StdErr << "ERROR: cmtk::SQLite threw exception - " << ex.what() << "\n";	
	}
      }
#endif
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkSimpleLevelsetCommandLine_h_included_

