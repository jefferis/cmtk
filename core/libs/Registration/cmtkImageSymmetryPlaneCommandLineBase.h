/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkImageSymmetryPlaneCommandLineBase_h_included_
#define __cmtkImageSymmetryPlaneCommandLineBase_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkImageSymmetryPlaneFunctionalBase.h>

#include <System/cmtkCannotBeCopied.h>
#include <System/cmtkCommandLine.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for symmetry plane computation command line tools.
 */
class ImageSymmetryPlaneCommandLineBase :
    private CannotBeCopied
{
public:
  /// Default constructor.
  ImageSymmetryPlaneCommandLineBase();

  /** Run the symmetry plane computation based on provided command line arguments.
   *\return The return code for the symmetry plane program. This should be returned
   * from main() via either "return" or "exit()".
   */
  int Run( const int argc, const char* argv[] );

  /// Get reference to command line parser object.
  CommandLine& GetCommandLine()
  {
    return this->m_CommandLine;
  }

protected:
  /// Create functional for volume.
  virtual ImageSymmetryPlaneFunctionalBase::SmartPtr CreateFunctional( UniformVolume::SmartPtr& volume ) = 0;

  /// Create functional for volume and value range.
  virtual ImageSymmetryPlaneFunctionalBase::SmartPtr CreateFunctional( UniformVolume::SmartPtr& volume, const Types::DataItemRange& range ) = 0;

private:
  /// Parse the given command line.
  bool ParseCommandLine ( const int argc, const char* argv[] );

  /// The symmetry plane object.
  ParametricPlane m_SymmetryPlane;

  /// Minimum data value (lower threshold).
  float m_MinValue;

  /// Flag for valid (user-set) minimum data value.
  bool m_MinValueSet;

  /// Minimum data value (lower threshold).
  float m_MaxValue;

  /// Flag for valid (user-set) maximum data value.
  bool m_MaxValueSet;

  /// Image sampling (highest resampled resolution).
  Types::Coordinate m_Sampling;

  /// Optimization "accuracy" (really, precision)
  Types::Coordinate m_Accuracy;
  
  /// Interpolation method for reformatted image generation.
  Interpolators::InterpolationEnum m_Interpolation;

  /// Number of multi-resolution levels.
  int m_Levels;
  
  /// Flag to disable optimization.
  bool m_DisableOptimization;

  /// Symmetry plane "Rho" parameter (offset).
  Types::Coordinate m_Rho;

  /// Symmetry plane "Theta" angle parameter.
  Units::Degrees m_Theta;

  /// Symmetry plane "Phi" angle parameter.
  Units::Degrees m_Phi;
  
  /// Flag to fix symmetry plane offset parameter (i.e., do not optimize "Rho")
  bool m_FixOffset;

  /// Optional output path for mirrored file.
  const char* m_MirrorOutFile;

  /// Optional output path for aligned file.
  const char* m_AlignedOutFile;

  /// Flag for marking the symmetry plane in the aligned file.
  bool m_MarkPlaneAligned;

  /// Optional output path for file with marked symmetry plane.
  const char* m_MarkedOutFile;
  
  /// Optional output path for file with subtraction between input and mirrored input.
  const char* m_DifferenceOutFile;

  /// Optional output path for the alignment transformation.
  const char* m_WriteXformPath;

  /// Data value used for marking the symmetry plane in output images.
  Types::DataItem m_MarkPlaneValue;

  /// Flag for user-provided padding value.
  bool m_PadOutValueSet;
  
  /// User-provided padding value.
  Types::DataItem m_PadOutValue;

  /// Optional output path for symmetry plane.
  const char* m_SymmetryOutFileName;
  
  /// Optional output path for symmetry plane parameters.
  const char* m_SymmetryParameters;

  /// Optional input path for previously computed symmetry plane parameters.
  const char* m_SymmetryParametersFile;

  /// Input image file path.
  const char* m_InFileName;
  
/// Constants for initial plane orientation.
  typedef enum
  {
    /// XY plane (axial)
    SYMPL_INIT_XY,
    /// XZ plane (coronal)
    SYMPL_INIT_XZ,
    /// YZ plane (sagittal)
    SYMPL_INIT_YZ
  } InitialPlaneEnum;
  
/// Initial plane orientation: default to sagittal for human images.
  InitialPlaneEnum m_InitialPlane;

  /// Write difference image between original and mirrored image.
  void WriteDifference( UniformVolume::SmartConstPtr& originalVolume ) const;

  /// Write mirrored image.
  void WriteMirror( UniformVolume::SmartConstPtr& originalVolume ) const;

  /// Write original image with marked symmetry plane.
  void WriteMarkPlane( UniformVolume::SmartConstPtr& originalVolume ) const;

  /// Write image aligned w.r.t. the symmetry plane.
  void WriteAligned( UniformVolume::SmartConstPtr& originalVolume ) const;

private:
  CommandLine m_CommandLine;
};

} // namespace cmtk

#endif // #ifndef  __cmtkImageSymmetryPlaneCommandLineBase_h_included_
