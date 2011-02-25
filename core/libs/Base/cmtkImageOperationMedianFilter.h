/*
//
//  Copyright 2009-2011 SRI International
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

#ifndef __cmtkImageOperationMedianFilter_h_included_
#define __cmtkImageOperationMedianFilter_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkDataGridFilter.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: median filtering.
class ImageOperationMedianFilter
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationMedianFilter( const int radiusX, const int radiusY, const int radiusZ ) : m_RadiusX( radiusX ), m_RadiusY( radiusY ), m_RadiusZ( radiusZ ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->SetData( DataGridFilter( volume ).GetDataMedianFiltered( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );    
    return volume;
  }
  
  /// Create a new median filter operation.
  static void New( const char* arg );
  
private:
  /// Downsampling radius in X direction.
  int m_RadiusX;

  /// Downsampling radius in Y direction.
  int m_RadiusY;

  /// Downsampling radius in Z direction.
  int m_RadiusZ;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationMedianFilter_h_included_
