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

#ifndef __cmtkSimpleLevelsetCommandLineBase_h_included_
#define __cmtkSimpleLevelsetCommandLineBase_h_included_

#include <cmtkconfig.h>

#include "Segmentation/cmtkSimpleLevelset.h"

#include "System/cmtkCommandLine.h"

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Command line interface base class for simple levelset segmentation.
 */
class SimpleLevelsetCommandLineBase
{
public:
  /// This class.
  typedef SimpleLevelsetCommandLineBase Self;

  /// Default constructor.
  SimpleLevelsetCommandLineBase();

  /// Initialize from command line arguments.
  int Init( const int argc, const char* argv[] );

  /// Reference to command line object.
  CommandLine& GetCommandLine()
  {
    return this->m_CommandLine;
  }

protected:
  ///Verbose mode.
  bool m_Verbose;

  /// Initial sphere scale factor.
  Types::Coordinate m_ScaleInitialSphere;

  /// Gaussian smoothing kernel sigma in mm.
  Types::Coordinate m_FilterSigma;

  /// Levelset evolution time constant.
  Types::Coordinate m_TimeDelta;

  /// Levelset threshold: the levelset function is truncated to plus/minus this value at each iteration.
  Types::Coordinate m_LevelsetThreshold;
  
  /// Number of levelset evolution iterations.
  int m_NumberOfIterations;
  
  /// Flag to force given number of iterations even when premature (discrete) convergence is detected.
  bool m_ForceIterations;

  /// Binarize levelset before output.
  bool m_Binarize;
  
  /// Input image path.
  const char* m_InFile;

  /// Output image path.
  const char* m_OutFile;

  /// The input image volume.
  UniformVolume::SmartConstPtr m_Volume;
  
#ifdef CMTK_USE_SQLITE
  /// Update this image/transformation database with the newly created levelset image.
  const char* m_UpdateDB;
#endif  

private:
  /// The command line parser object.
  cmtk::CommandLine m_CommandLine;
};

} // namespace cmtk

#endif // #ifndef __cmtkSimpleLevelsetCommandLineBase_h_included_
