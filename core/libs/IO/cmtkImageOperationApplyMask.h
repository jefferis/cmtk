/*
//
//  Copyright 2009-2010 SRI International
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

#ifndef __cmtkImageOperationApplyMask_h_included_
#define __cmtkImageOperationApplyMask_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkImageOperation.h"

namespace
cmtk
{

/// Apply mask image.
class ImageOperationApplyMask
/// Inherit generic image operation.
  : public ImageOperation
{
public:
  /// This class.
  typedef ImageOperationApplyMask Self;

  /// Constructor.
  ImageOperationApplyMask( const cmtk::UniformVolume::SmartPtr& maskVolume ) : m_MaskVolume( maskVolume ) {}

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );

  /// Create new mask operation.
  static void New( const char* maskFileName )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationApplyMask( Self::ReadMaskFile( maskFileName ) ) ) );
  }

  /// Create new inverse mask operation.
  static void NewInverse( const char* maskFileName )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationApplyMask( Self::ReadMaskFile( maskFileName, true /*inverse*/ ) ) ) );
  }

private:
  /// The mask volume.
  cmtk::UniformVolume::SmartPtr m_MaskVolume;

  /// Read the actual mask file.
  static cmtk::UniformVolume::SmartPtr ReadMaskFile( const char* maskFileName, const bool inverse = false );
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationApplyMask_h_included_
