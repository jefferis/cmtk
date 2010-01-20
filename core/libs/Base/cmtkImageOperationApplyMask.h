/*
//
//  Copyright 2009 SRI International
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

#include <cmtkImageOperation.h>

namespace
cmtk
{

/// Apply mask image.
class ImageOperationApplyMask
/// Inherit generic image operation.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationApplyMask( const cmtk::UniformVolume::SmartPtr& maskVolume ) : m_MaskVolume( maskVolume ) {}

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    const std::string maskOrientation = this->m_MaskVolume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
    const std::string workingOrientation = volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
    if ( maskOrientation != workingOrientation )
      {
      this->m_MaskVolume = cmtk::UniformVolume::SmartPtr( this->m_MaskVolume->GetReoriented( workingOrientation.c_str() ) );
      }
    
    for ( int dim = 0; dim < 3; ++dim )
      {
      if ( this->m_MaskVolume->m_Dims[dim] != volume->m_Dims[dim] )
	{
	cmtk::StdErr << "ERROR: mask volume dimensions do not match working volume dimensions.\n";
	exit( 1 );
	}
      }
    
    const cmtk::TypedArray* maskData = this->m_MaskVolume->GetData();
    cmtk::TypedArray::SmartPtr& volumeData = volume->GetData();

    const size_t nPixels = volume->GetNumberOfPixels();
    for ( size_t i = 0; i < nPixels; ++i )
      if ( maskData->IsPaddingOrZeroAt( i ) ) 
	volumeData->SetPaddingAt( i );
    return volume;
  }

  /// Create new mask operation.
  static void New( const char* maskFileName )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationApplyMask( ReadMaskFile( maskFileName ) ) ) );
  }

  /// Create new inverse mask operation.
  static void NewInverse( const char* maskFileName )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationApplyMask( ReadMaskFile( maskFileName, true /*inverse*/ ) ) ) );
  }

private:
  /// The mask volume.
  cmtk::UniformVolume::SmartPtr m_MaskVolume;

  /// Read the actual mask file.
  static cmtk::UniformVolume::SmartPtr ReadMaskFile( const char* maskFileName, const bool inverse = false )
  {
    cmtk::UniformVolume::SmartPtr maskVolume( cmtk::VolumeIO::ReadOriented( maskFileName ) );
    if ( !maskVolume || !maskVolume->GetData() ) 
      {
      cmtk::StdErr << "ERROR: could not read mask from file " << maskFileName << "\nProgram will terminate now, just to be safe.\n";
      exit( 1 );
      }
    
    // binarize mask to 1/0, convert to char, and also consider "inverse" flag in the process.
    cmtk::TypedArray::SmartPtr& maskData = maskVolume->GetData();
    const size_t nPixels = maskData->GetDataSize();
    for ( size_t n = 0; n < nPixels; ++n )
      {
      if ( maskData->IsPaddingOrZeroAt( n ) != inverse ) 
	maskData->Set( n, 1 );
      else
	maskData->Set( n, 0 );
      }
    maskVolume->SetData( cmtk::TypedArray::SmartPtr( maskData->Convert( cmtk::TYPE_BYTE ) ) );
    
    return maskVolume;
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationApplyMask_h_included_
