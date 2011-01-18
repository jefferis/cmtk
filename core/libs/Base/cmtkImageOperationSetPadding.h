/*
//
//  Copyright 2010 Torsten Rohlfing
//
//  Copyright 2011 SRI International
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

#ifndef __cmtkImageOperationSetPadding_h_included_
#define __cmtkImageOperationSetPadding_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkDataGridMorphologicalOperators.h>

namespace
cmtk
{

/// Image operation: set padding flag and value.
class ImageOperationSetPadding
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationSetPadding( const bool flag, const double value = 0 ) : m_PaddingFlag( flag ), m_PaddingValue( value ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    if ( this->m_PaddingFlag )
      {
      volume->GetData()->SetPaddingValue( this->m_PaddingValue );
      }
    else
      {
      volume->GetData()->ClearPaddingFlag();
      }

    return volume;
  }
  
  /// Create new connected components operation.
  static void New( const double value )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationSetPadding( true, value ) ) );
  }

  /// Create new connected components operation.
  static void NewUnset()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationSetPadding( false ) ) );
  }

private:
  /// Set flag: if this is set, padding is activated, otherwise it is deactivated (and m_PaddingValue is ignored).
  bool m_PaddingFlag;

  /// Padding value.
  Types::DataItem m_PaddingValue;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationSetPadding_h_included_
