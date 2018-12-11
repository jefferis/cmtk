/*
//
//  Copyright 2010 Torsten Rohlfing
//
//  Copyright 2011, 2014 SRI International
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
//  $Revision: 2902 $
//
//  $LastChangedDate: 2011-02-24 12:12:46 -0800 (Thu, 24 Feb 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkImageOperationSetDataClass_h_included_
#define __cmtkImageOperationSetDataClass_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkImageOperation.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: set padding flag and value.
class ImageOperationSetDataClass
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationSetDataClass( const DataClass dataClass ) : m_DataClass( dataClass ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->GetData()->SetDataClass( this->m_DataClass );
    return volume;
  }
  
  /// Create new operation to set class to labels.
  static void NewLabels()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationSetDataClass( DATACLASS_LABEL ) ) );
  }

  /// Create new operation to set class to grey values.
  static void NewGrey()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationSetDataClass( DATACLASS_GREY ) ) );
  }

private:
  /// Set flag: if this is set, padding is activated, otherwise it is deactivated (and m_DataClassValue is ignored).
  DataClass m_DataClass;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationSetDataClass_h_included_
