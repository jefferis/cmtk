/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkImageTemplate_h_included_
#define __cmtkImageTemplate_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Uniform volume template.
 * This class is a templated generalization of the UniformVolume class. Here, the type of pixel data is determined by template instantiation.
 */
template<class TPixelType>
class ImageTemplate : 
  /// Inherit from generic Volume class.
  public UniformVolume 
{
public:
  /// Pixel data type.
  typedef TPixelType PixelType;

  /// This class.
  typedef ImageTemplate<TPixelType> Self;

  /// Superclass.
  typedef UniformVolume Superclass;

  /// Smart pointer to ImageTemplate.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const ImageTemplate.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Region type.
  typedef Superclass::CoordinateRegionType CoordinateRegionType;

  /// Index type.
  typedef Superclass::CoordinateVectorType CoordinateVectorType;

  /// Destructor.
  virtual ~ImageTemplate() {}

  /** Create volume "from scratch".
   *\param dims Number of grid elements for the three spatial dimensions.
   *\param size Size of the volume in real-world coordinates.
   */
  ImageTemplate( const DataGrid::IndexType& dims, const Self::CoordinateVectorType& size ) : Superclass( dims, size ) 
  {
    this->m_DataArray.resize( this->GetNumberOfPixels );
  }

  /** Create volume from base class instance.
   *\param dims Number of grid elements for the three spatial dimensions.
   *\param size Size of the volume in real-world coordinates.
   */
  ImageTemplate( const Superclass& from ) : Superclass( from ) 
  {
  }

  /// Access operator.
  PixelType& operator[]( const size_t idx )
  {
    return this->m_DataArray[idx];
  }

  /// Const access operator.
  const PixelType& operator[]( const size_t idx ) const
  {
    return this->m_DataArray[idx];
  }
  
private:
  /// The actual data array.
  std::vector<PixelType> m_DataArray;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageTemplate_h_included_
