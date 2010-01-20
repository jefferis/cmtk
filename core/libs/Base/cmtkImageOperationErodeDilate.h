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

#ifndef __cmtkImageOperationErodeDilate_h_included_
#define __cmtkImageOperationErodaDilate_h_included_

#include <cmtkconfig.h>

#include <cmtkImageOperation.h>

namespace
cmtk
{

/// Image operation: erode or dilate.
class ImageOperationErodeDilate
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationErodeDilate( const int iterations ) : m_Iterations( iterations ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    if ( this->m_Iterations < 0 )
      volume->ApplyErode( -this->m_Iterations );
    else
      if ( this->m_Iterations > 0 )
	volume->ApplyDilate( this->m_Iterations );
    return volume;
  }

  /// Create new dilation operation.
  static void NewDilate( const long int iterations )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationErodeDilate( iterations ) ) );
  }

  /// Create new erosion operation.
  static void NewErode( const long int iterations )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationErodeDilate( -iterations ) ) );
  }
  
private:
  /// Number of iterations of erosion (if negative) or dilation (if positive).
  int m_Iterations;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationErodeDilate_h_included_
