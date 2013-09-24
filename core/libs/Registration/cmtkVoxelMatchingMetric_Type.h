/*
//
//  Copyright 1997-2012 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#ifndef __cmtkVoxelMatchingMetric_Type_h_included_
#define __cmtkVoxelMatchingMetric_Type_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkDataTypeTraits.h>
#include <Base/cmtkInterpolator.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypedArray.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for voxel metrics with pre-converted image data.
 */
template<class T,ScalarDataType DT>
class VoxelMatchingMetric_Type
{
public:
  /// This class.
  typedef VoxelMatchingMetric_Type<T,DT> Self;

  /// This is the type used internally for storing pre-processed image data.
  typedef T Exchange;
  
  /// Structure for handling the two images compared.
  class ImageData 
  {
  private:
    /// Padding data used in this instance.
    T Padding;

  public:
    /// This function returns the constant actually used for padding data.
    T padding() const { return this->Padding; }
    
    /// Precomputed reference voxel data bin indices.
    T* Data;
    
    /// Reference-counting wrapper for Data.
    TypedArray::SmartPtr DataArray;

    /// Bin offset.
    Types::DataItem BinOffset;
    
    /// Bin width.
    Types::DataItem BinWidth;
    
    /// Value range.
    Types::DataItemRange m_ValueRange;
    
    /// Number of pixels per dimension in the original image.
    DataGrid::IndexType ImageDims;

    /// Total number of pixels.
    size_t NumberOfSamples;

    /// Get value range of distribution as stored herein.
    const Types::DataItemRange GetValueRange() const
    {
      return this->m_ValueRange;
    }

    /// Convert value to bin.
    byte ValueToIndex( const Types::DataItem value )
    {
      return static_cast<byte>( (std::min( std::max( value, this->m_ValueRange.m_LowerBound ), this->m_ValueRange.m_UpperBound )- this->BinOffset) / BinWidth );
    }
    
    /// Default constructor.
    ImageData() : 
      Padding( DataTypeTraits<T>::ChoosePaddingValue() ), 
      Data( NULL ), 
      DataArray( NULL ), 
      m_ValueRange( 0, 0 ),
      NumberOfSamples( 0 )
    {
      nextJ = nextK = nextIJ = nextJK = nextIK = nextIJK = 0;
    }

    /// Allocate internal data array and create wrapper for reference counting.
    void AllocDataArray( const TypedArray* templateArray );

    /** Initialize internal storage for one volume.
     * The volume's data is converted into an array of byte values that
     * directly represent the bin to sort the respective sample into. In 
     * addition, the number of bins is determined and the bins array allocated.
     *\param volume The original volume data.
     *\param defNumBins The desired number of bins. If this parameter is
     * zero, a suitable number is automatically determined.
     *\param bounds User-specified bounds for data values. All values
     * outside this range will be set to the upper or lower limit and sorted into the
     * first or last histogram bin, respectively.
     *\return The number of bins required to hold the data. Note that there 
     * will be an extra bin allocated to hold non-existing data values.
     */
    size_t Init( const UniformVolume* volume, const size_t defNumBins, const Types::DataItemRange& bounds = Types::DataItemRange( -HUGE_VAL, HUGE_VAL ) );
    
    /** Initialize internal storage for one (reference of model) volume.
     * The volume's data is converted into an array of byte values that 
     * directly represent the bin to sort the respective sample into. In 
     * addition, the number of bins is determined and the bins array allocated.
     * This function can distinguish between different kinds of data 
     * (grey-level, binary, and labels) and handle these accordingly.
     *\param volume The original volume data.
     */
    void Init( const UniformVolume* volume );
    
    /// Precompute grid increments.
    void PrecomputeIncrements( const UniformVolume* volume );

    /// Offset of next voxel row.
    unsigned int nextJ;
    
    /// Offset for next row and column.
    unsigned int nextIJ;
    
    /// Offset for next plane.
    unsigned int nextK;
    
    /// Offset for next plane and column.
    unsigned int nextIK;
    
    /// Offset for next plane and row.
    unsigned int nextJK;
    
    /// Offset for next plane, row, and column.
    unsigned int nextIJK;
  };

  /// Data of the X distribution.
  ImageData DataX;

  /// Data of the Y distribution.
  ImageData DataY;

  /// Set data for the X distribution (reference image).
  void SetDataX( const UniformVolume* volume )
  { 
    DataX.PrecomputeIncrements( volume ); 
  }

  /// Set data for the Y distribution (floating image).
  void SetDataY( const UniformVolume* volume ) 
  { 
    DataY.PrecomputeIncrements( volume );
  }

  /// Set data for the X and Y distribution (both images).
  void SetDataXY( const UniformVolume* volumeX, const UniformVolume* volumeY ) 
  {
    DataX.PrecomputeIncrements( volumeX );
    DataY.PrecomputeIncrements( volumeY );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVoxelMatchingMetric_Type_h_included_
