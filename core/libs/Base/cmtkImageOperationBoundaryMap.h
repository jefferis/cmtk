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
//  $Revision: 1150 $
//
//  $LastChangedDate: 2010-01-18 12:41:38 -0800 (Mon, 18 Jan 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkImageOperationBoundaryMap_h_included_
#define __cmtkImageOperationBoundaryMap_h_included_

#include <cmtkconfig.h>

#include <cmtkImageOperation.h>

namespace
cmtk
{

/// Image operation: create binary or multi-valued boundary map.
class ImageOperationBoundaryMap
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationBoundaryMap( const bool multiValued ) : m_MultiValued( multiValued ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->SetData( cmtk::TypedArray::SmartPtr( volume->GetBoundaryMap( this->m_MultiValued ) ) );
    return volume;
  }
  
  /// Create new binary boundary map operation.
  static void New()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationBoundaryMap( false ) ) );
  }

  /// Create new multi-valued boundary map operation.
  static void NewMulti()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationBoundaryMap( true ) ) );
  }

private:
  /// Multi-valued flag: if this is set, a multi-valued boundary map will be created, otherwise a binary map.
  bool m_MultiValued;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationBoundaryMap_h_included_
