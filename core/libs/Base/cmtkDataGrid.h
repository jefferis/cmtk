/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#ifndef __cmtkDataGrid_h_included_
#define __cmtkDataGrid_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkMacros.h"
#include "Base/cmtkTypes.h"
#include "Base/cmtkTypedArray.h"
#include "Base/cmtkFixedVector.h"
#include "Base/cmtkScalarImage.h"
#include "Base/cmtkRegion.h"
#include "Base/cmtkMetaInformationObject.h"
#include "Base/cmtkAnatomicalOrientation.h"

#include "System/cmtkSmartPtr.h"
#include "System/cmtkThreads.h"

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Grid topology of data arranged in a 3D lattice.
 * This class extends the plain data handling functions of TypedArray
 * with a 3D topology. Real world coordinates, however, are not considered and
 * need to be handled by derived classes. Thus, this class provides the coordinate
 * independent services such as median filtering and, to a certain extent,
 * interpolation.
 */
class DataGrid :
  /// Inherit class that handles meta information.
  public MetaInformationObject
{
public:
  /// This class.
  typedef DataGrid Self;

  /// Smart pointer to DataGrid.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const to DataGrid.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Region type.
  typedef Region<3,int> RegionType;

  /// Index type.
  typedef RegionType::IndexType IndexType;

  /// Space vector type.
  typedef FixedVector<3,Types::Coordinate> SpaceVectorType;

  /// Number of grid samples in the three spatial dimensions
  Self::IndexType m_Dims;

  /// Data array (element type is variable)
  cmtkGetSetMacro(TypedArray::SmartPtr,Data);

public:
  /// Default constructor.
  DataGrid() : 
    m_Data( NULL )
  {}
  
  /// Constructor.
  DataGrid( const Self::IndexType& dims ) 
    : m_Dims( dims ), 
      m_Data( NULL )
  {
    this->m_CropRegion = this->GetWholeImageRegion();
  }
  
  /// Virtual destructor.
  virtual ~DataGrid() {}

  /// Test whether this grid matches another one, i.e., has the same dimensions.
  bool GridMatches( const Self& other ) const
  {
    return (this->m_Dims == other.m_Dims);
  }

  /// Downsampling and pixel-averaging constructor function.
  virtual DataGrid* GetDownsampledAndAveraged( const int (&downsample)[3] ) const;

  /** Reorientation constructor function.
   *@param newOrientation Three letter orientation code that specifies the anatomically-based
   * orientation of the reoriented volume. Each letter can be one of the following: R, L, A, 
   * P, I, S. These stand for Right, Left, Anterior, Posterior, Inferior, Superior. 
   *
   * The three letters in the orientation string define the directions of the positive x, y, 
   * and z axes, in this order. For example, "RAS", the standard orientation for this software, 
   * means that the pixels along the x axis are arranged from the subject's Left to the Right 
   * side, along the y axis from the subject's Posterior (back) to Anterior (front), and along
   * the z axis from Inferior (feet) to Superior (head).
   *
   * The current orientation of this volume is to be taken from its meta information,
   * as this->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION]. This is also a three-letter string of the
   * same form as the one given to this function as a parameter.
   *
   * If the current orientation is not set, a warning message should be printed to StdErr, and
   * a NULL pointer returned.
   *
   *\return Reoriented data grid with permuted pixels in this->Data and permuted grid dimensions
   * in this->Dims. The returned pointers points to a newly allocated object, which can be
   * wrapped in an SmartPointer.
   */
  const DataGrid::SmartPtr GetReoriented( const char* newOrientation = AnatomicalOrientation::ORIENTATION_STANDARD ) const;
  
  /** Set dimensions array.
   * This function updates the internal offsets for fast access to adjacent
   * voxel rows, columns, planes etc.
   */
  void SetDims( const Self::IndexType& dims );

  /// Get dimensions array.
  const Self::IndexType GetDims() const
  {
    return this->m_Dims;
  }

  /** Create data array.
   *\param dataType ID of scalar data type for the array. This is the image pixel type.
   *\param setToZero If this flag is set, all values in the newly created array will be initialized to zero.
   */
  virtual TypedArray::SmartPtr CreateDataArray( const ScalarDataType dataType, const bool setToZero = false );

  /// Get number of data items in the volume.
  size_t GetNumberOfPixels () const { return this->m_Dims[0]*this->m_Dims[1]*this->m_Dims[2]; }

  /// Check whether given pixel index is inside grid.
  bool IndexIsInRange( const int x, const int y, const int z ) const
  {
    return (x>=0) && (x<this->m_Dims[0]) && (y>=0) && (y<this->m_Dims[1]) && (z>=0) && (z<this->m_Dims[2]);
  }

  /// Get offset of a pixel.
  size_t GetOffsetFromIndex( const int x, const int y, const int z ) const 
  {
    return x + nextJ * y + nextK * z;
  }

  /// Get index of a pixel identified by its offset.
  void GetIndexFromOffset( const size_t offset, int& x, int& y, int& z ) const 
  {
    z = offset / nextK;
    y = (offset % nextK) / nextJ;
    x = (offset % nextK) % nextJ;
  }

  /// Return data at specified offset
  bool GetDataAt ( Types::DataItem& data, const size_t offset ) const 
  {
    return this->m_Data->Get( data, offset );
  }

  /// Set data at specified offset
  void SetDataAt ( const Types::DataItem data, const size_t offset )
  {
    this->m_Data->Set( data, offset );
  }
  
  /// Return data at specified grid point.
  bool GetDataAt ( Types::DataItem& data, const int x, const int y, const int z ) const
  {
    return this->GetDataAt( data, this->GetOffsetFromIndex( x, y, z ) );
  }

  /// Set data at specified grid point.
  void SetDataAt ( const Types::DataItem data, const int x, const int y, const int z )
  {
    this->SetDataAt( data, this->GetOffsetFromIndex( x, y, z ) );
  }

  /// Return data at specified grid point, or a given default value if no data exists there.
  Types::DataItem GetDataAt ( const int x, const int y, const int z, const Types::DataItem defaultValue = 0.0 ) const
  {
    Types::DataItem value;
    if ( this->GetDataAt( value, this->GetOffsetFromIndex( x, y, z ) ) )
      return value;
    else
      return defaultValue;
  }

  /// Return data at specified grid offset, or a given default value if no data exists there.
  Types::DataItem GetDataAt ( const size_t offset, const Types::DataItem defaultValue = 0.0 ) const
  {
    Types::DataItem value;
    if ( this->GetDataAt( value, offset ) )
      return value;
    else
      return defaultValue;
  }

  /** Return data after mirroring.
   *@param axis Coordinate axis normal to mirror plane. Default is AXIS_X
   * (mirror about mid-sagittal plane).
   */
  TypedArray::SmartPtr GetDataMirrorPlane( const int axis = AXIS_X ) const;

  /// Replace data with mirrored version.
  void ApplyMirrorPlane( const int axis = AXIS_X );

private:
  /// Mirror about plane without allocating additional memory.
  static void MirrorPlaneInPlace( TypedArray& data, const Self::IndexType& dims, const int axis = AXIS_X );

public:
  /** Get cropped region reference.
   */
  Self::RegionType& CropRegion()
  {
    return this->m_CropRegion;
  }

  /** Set cropped region.
   * This function can deal with negative crop region values, which refer to the upper grid
   * boundary and are automatically converted.
   */
  virtual void SetCropRegion( const Self::RegionType& region );

  /** Get cropped region reference.
   */
  const Self::RegionType& CropRegion() const
  {
    return this->m_CropRegion;
  }

  /** Return number of voxels in the cropped remaining image.
   */
  int GetCropRegionNumVoxels() const;

  /// Get whole image region.
  const RegionType GetWholeImageRegion() const;

  /// Get index increments for crop region.
  const Self::IndexType GetCropRegionIncrements() const;

  /** Automatically crop to minimal region above given threshold.
   *@param threshold The cropping threshold.
   *@param recrop If this flag is true, then the cropping will be performed
   * inside an already existing cropping region. If this flag is false 
   * (default), then any pre-set crop region is ignored.
   */
  void AutoCrop( const Types::DataItem threshold, const bool recrop = false, const int margin = 0 );

  /// Fill volume outside current crop region with constant value.
  void FillCropBackground( const Types::DataItem value );

  /// Return data for cropped volume.
  TypedArray::SmartPtr GetCroppedData() const;

  /// Accessor functions for protected member variables
  int GetNextI() const { return nextI; }
  int GetNextJ() const { return nextJ; }
  int GetNextK() const { return nextK; }
  int GetNextIJ() const { return nextIJ; }
  int GetNextIK() const { return nextIK; }
  int GetNextJK() const { return nextJK; }
  int GetNextIJK() const { return nextIJK; }
  
  /// Get center of mass of pixel data.
  virtual FixedVector<3,Types::Coordinate> GetCenterOfMassGrid() const;
  
  /// Get center of mass and first-order moments of pixel data.
  virtual FixedVector<3,Types::Coordinate> GetCenterOfMassGrid( FixedVector<3,Types::Coordinate>& firstOrderMoment ) const;
  
  /// Return maximum-intensity projection MIP along one axis.
  template<class TAccumulator>
  ScalarImage* ComputeProjection( const int axis ) const;
  
  /** Return orthogonal slice as a 2D image.
   */
  virtual ScalarImage* GetOrthoSlice( const int axis, const unsigned int plane ) const;
  
  /** Set orthogonal slice from a 2D image.
   */
  virtual void SetOrthoSlice( const int axis, const unsigned int idx, const ScalarImage* slice );

  /// Return data after mirror operator.
  TypedArray::SmartPtr GetDataMirrored( const int axis = AXIS_X ) const;

  /// Print object.
  void Print() const;

protected:
  /** Utility function for trilinear interpolation.
   *@param data This reference is set to the interpolated data value. It is 
   * valid if and only if this function returns 1.
   *@param location 3D coordinate to interpolate data at.
   *@param gridPosition (x,y,z) indices of the voxel containing the given
   * location.
   *@param cellFrom 3D coordinate of the lower-left-front voxel of the cell
   * enclosing the given location.
   *@param cellTo 3D coordinate of the upper-right-rear voxel of the cell
   * enclosing the given location.
   *@return True if there is valid data for all eight voxels enclosing the 
   * given location, so that the interpolation could be completed successfully,
   * False otherwise.
   */
  bool TrilinearInterpolation( Types::DataItem&, const int, const int, const int, const Self::SpaceVectorType&, const Types::Coordinate*, const Types::Coordinate* ) const;

  /** Utility function for trilinear interpolation from a primitive data array.
   * This function is provided for computational efficiency when a large number 
   * of interpolations from a given data volume of known pixel type are required.
   *
   *\param dataPtr Pointer to the primitive data array.
   *@param location 3D coordinate to interpolate data at.
   *@param gridPosition (x,y,z) indices of the voxel containing the given
   * location.
   *@param cellFrom 3D coordinate of the lower-left-front voxel of the cell
   * enclosing the given location.
   *@param cellTo 3D coordinate of the upper-right-rear voxel of the cell
   * enclosing the given location.
   *@return The interpolated data value..
   */
  template<class TData>
  inline TData TrilinearInterpolation
  ( const TData* dataPtr, const int x, const int y, const int z, const Self::SpaceVectorType& gridPosition, const Types::Coordinate* cellFrom, 
    const Types::Coordinate* cellTo ) const;

  /** Utility function for trilinear interpolation from multiple primitive data arrays of identical grid structure.
   * This function is provided for computational efficiency when a large number 
   * of interpolations from a given data volume of known pixel type are required.
   */
  template<class TData,class TOutputIterator>
  inline void TrilinearInterpolation
  ( TOutputIterator result, const std::vector<TData*>& dataPtr, const int x, const int y, const int z,
    const Types::Coordinate fracX, const Types::Coordinate fracY, const Types::Coordinate fracZ ) const;

  /// Offset to next voxel column.
  int nextI;

  /// Offset to next voxel row.
  int nextJ;

  /// Offset to next voxel plane.
  int nextK;

  /// Offset to next column and row.
  int nextIJ;
  
  /// Offset to next column and plane.
  int nextIK;

  /// Offset to next row and plane.
  int nextJK;

  /// Offset to next column, row, and plane.
  int nextIJK;

private:
  /** Crop region.
   */
  Self::RegionType m_CropRegion;
};

//@}

} // namespace cmtk

#include "cmtkDataGrid.txx"

#endif // #ifndef __cmtkDataGrid_h_included_
