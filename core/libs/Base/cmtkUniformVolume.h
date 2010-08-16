/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkUniformVolume_h_included_
#define __cmtkUniformVolume_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkVolume.h"

#include "System/cmtkSmartPtr.h"
#include "System/cmtkSmartConstPtr.h"

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Forward declaration.
class VolumeGridToGridLookup;

/** Uniform volume.
 * This class handles 3D isotropic data of any scalar type.
 */
class UniformVolume : 
  /// Inherit from generic Volume class.
  public Volume 
{
public:
  /// This class.
  typedef UniformVolume Self;

  /// Superclass.
  typedef Volume Superclass;

  /// Smart pointer to UniformVolume.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const UniformVolume.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Region type.
  typedef Superclass:: CoordinateRegionType CoordinateRegionType;

  /// Index type.
  typedef Superclass::CoordinateVectorType CoordinateVectorType;

  /// Grid spacing for three dimensions
  CoordinateVectorType m_Delta; 

  /// Get pixel size vector.
  const CoordinateVectorType& Deltas() const
  {
    return this->m_Delta;
  }

  /// Get pixel size vector.
  CoordinateVectorType& Deltas()
  {
    return this->m_Delta;
  }

  /// Destructor.
  virtual ~UniformVolume() {}

  /** Resample other volume to given resolution.
   *\param other Original volume that will be resampled.
   *\param resolution Resolution of the newly created volume in world units.
   *\param allowUpsampling If true, then local upsampling is allowed in regions
   * where the original image resolution (non-uniform) was coarser than the
   * given resolution of the resampled volume.
   */
  UniformVolume( const UniformVolume& other, const Types::Coordinate resolution = 0, const bool allowUpsampling = false );

  /** Create uniform volume "from scratch".
   *@param dims Number of grid elements for the three spatial dimensions.
   *@param size Size of the volume in real-world coordinates.
   *@param data An existing TypedArray containing the scalar voxel data.
   */
  UniformVolume( const DataGrid::IndexType& dims, const Self::CoordinateVectorType& size, TypedArray::SmartPtr& data = TypedArray::SmartPtr::Null );

  /** Create uniform volume "from scratch".
   *@param dims Number of grid elements for the three spatial dimensions.
   *@param deltaX Pixel size in x direction.
   *@param deltaY Pixel size in y direction.
   *@param deltaZ Pixel size in z direction.
   *@param data An existing TypedArray containing the scalar voxel data.
   */
  UniformVolume( const DataGrid::IndexType& dims, const Types::Coordinate deltaX, const Types::Coordinate deltaY, const Types::Coordinate deltaZ, TypedArray::SmartPtr& data = TypedArray::SmartPtr::Null );

  /// Test whether this grid matches another one, i.e., has the same pixel sizes.
  bool GridMatches( const Self& other ) const
  {
    return Superclass::GridMatches( other ) && ((this->m_Delta-other.m_Delta).MaxAbsValue() < 1e-5) && ((this->m_Offset-other.m_Offset).MaxAbsValue() < 1e-5);
  }

  /** Coordinate transformation from index to physical position.
   * This incorporates image axis directions and first pixel offset.
   *\note Strictly, this is not a transformation from the pixel index
   *  to physical coordinates, but from pixel index times pixel size to
   *  physical coordinates. This way, the transformation maps from physical
   *  space back into physical space.
   */
  AffineXform::MatrixType m_IndexToPhysicalMatrix;

  /** Change volume coordinate space.
   * Re-arrange volume's direction vectors to refer to a different coordinate space.
   *\invariant This should keep the result of GetOrientationFromDirections invariant.
   */
  virtual void ChangeCoordinateSpace( const std::string& newSpace );
  
  /** Get anatomical orientation from space directions.
   */
  std::string GetOrientationFromDirections() const;

  /** Create a default index-to-physical transformation matrix.
   * The new matrix is the identity matrix with diagonal elements scaled according to pixel size.
   */
  virtual void CreateDefaultIndexToPhysicalMatrix();

  /** Get physical location of indexed pixel.
   * This function applies the index-to-physical transformation matrix to the given pixel index. 
   * The index itself can be fractional.
   */
  virtual Self::CoordinateVectorType IndexToPhysical( const Types::Coordinate i, const Types::Coordinate j, const Types::Coordinate k ) const;

  /** Get matrix that maps form image coordinates to physical space.
   * The returned matrix is computed by removing pixel size from the one stored
   * in m_IndexToPhysicalMatrix.
   */
  virtual AffineXform::MatrixType GetImageToPhysicalMatrix() const;

  /** Create a physical copy of this object.
   *@param copyData If true, the associated data array is also copied.
   */
  Self::SmartPtr Clone( const bool copyData )
  {
    return Self::SmartPtr( this->CloneVirtual( copyData ) );
  }

  Self::SmartPtr Clone() const
  {
    return Self::SmartPtr( this->CloneVirtual() );
  }

  /// Create igsUniformObject with identical geometry but no data.
  Self::SmartPtr CloneGrid() const
  {
    return Self::SmartPtr( this->CloneGridVirtual() );
  }

  virtual Types::Coordinate GetDelta( const int axis, const int = 0 ) const
  {
    return this->m_Delta[axis];
  }

  /** Get minimum grid spacing for all dimensions.
   * For each dimension, the minimum delta is retrieved by a call to the
   * derived class's corresponding method.
   */
  virtual Types::Coordinate GetMinDelta () const 
  {
    return this->m_Delta.MinValue();
  }
  
  /** Get maximum grid spacing for all dimensions.
   * For each dimension, the maximum delta is retrieved by a call to the
   * derived class's corresponding method.
   */
  virtual Types::Coordinate GetMaxDelta () const 
  {
    return this->m_Delta.MaxValue();
  }

  /** Resample other volume.
   * Resampling is done by a sliding window algorithm combining both
   * interpolation and averaging, even within one volume.
   *@return Pointer to a newly created TypedArray object holding the
   * resampled volume data. NULL is returned if sufficient memory is not
   * available.
   */
  virtual TypedArray::SmartPtr Resample ( const UniformVolume& ) const;

  /// Get volume reoriented to a different anatomical axis alignment.
  const UniformVolume::SmartPtr GetReoriented ( const char* newOrientation = AnatomicalOrientation::ORIENTATION_STANDARD ) const;

  /** Downsampling constructor function.
   *\param approxIsotropic If this is set (default: off), then the downsample
   * factors per dimension are adjusted so that the resulting output volume
   * is as close to isotropic as possible without interpolation.
   */
  virtual UniformVolume* GetDownsampledAndAveraged( const int downsample, const bool approxIsotropic = false ) const;
 
  /// Downsampling constructor function.
  virtual UniformVolume* GetDownsampledAndAveraged( const int (&downsample)[3] ) const;

  /// Get interleaved sub-volume along given axis and with given interleave offset.
  UniformVolume* GetInterleavedSubVolume( const int axis, const int factor, const int idx ) const;

  /// Get interleaved sub-volume along given axis and with given interleave offset, padded with empty image planes.
  UniformVolume* GetInterleavedPaddedSubVolume( const int axis, const int factor, const int idx ) const;

  /// Mirror volume and associated data.
  virtual void Mirror ( const int axis = AXIS_X );

  /** Return orthogonal slice.
   * This function calls its equivalent in DataGrid and adds calibration
   * info (i.e., correct pixel sizes) to the resulting image.
   *\NOTE The pixel size if taken from the size of the first grid element along
   * each axis -- non-uniform spacings will lead to incorrect results.
   */
  virtual ScalarImage* GetOrthoSlice( const int axis, const unsigned int plane ) const;

  /** Extract orthogonal slice as a new volume.
   */
  Self::SmartPtr ExtractSlice( const int axis, const unsigned int plane ) const
  {
    return Self::SmartPtr( this->ExtractSliceVirtual( axis, plane ) );
  }

  /** Return interpolated orthogonal slice.
   * This function calls its non-interpolating counterpart twice and performs
   * a 1-D interpolation on the results to form the interpolated slice.
   *\NOTE The pixel size if taken from the size of the first grid element along
   * each axis -- non-uniform spacings will lead to incorrect results.
   */
  virtual ScalarImage* GetOrthoSliceInterp( const int axis, const Types::Coordinate location ) const;

  /** Return orthogonal slice by location.
   * This function looks up a given orthogonal slice location and returns the 
   * nearest slice from this volume.
   *\NOTE The pixel size if taken from the size of the first grid element along
   * each axis -- non-uniform spacings will lead to incorrect results.
   */
  virtual ScalarImage* GetNearestOrthoSlice( const int axis, const Types::Coordinate location ) const;

  /** Get date gradient vector at pixel using central differences.
   * This function cannot be called for pixels on the volume boundaries, i.e.,
   * we require that 0 < i,j,k < [Dims[0]-1,Dims[1]-1,Dims[2]-1].
   *
   *\warning No parameter range checking is currently performed!
   */
  virtual const Self::CoordinateVectorType GetGradientAt( const int i, const int j, const int k );

  /** Get data Hessian matrix at pixel using central differences.
   * This function cannot be called for pixels within a two-pixel distance from the
   * volume boundary, i.e., we require that 1 < i,j,k < [Dims[0]-2,Dims[1]-2,Dims[2]-2].
   *
   *\warning No parameter range checking is currently performed!
   */
  virtual void GetHessianAt( Matrix3x3<Types::DataItem>& H, const int i, const int j, const int k );

  /// Get data value at specified coordinate.
  template<class TData> inline bool ProbeData( TData& result, const TData* dataPtr, const Self::CoordinateVectorType& location ) const;

  /// Get data values from multiple arrays at specified coordinate.
  template<class TData,class TOutputIterator> inline bool ProbeData
  ( const TOutputIterator& result, const std::vector<TData*>& dataPtr, const Self::CoordinateVectorType& location ) const;

  /// Return linearly interpolated voxel without applying a transformation.
  inline bool ProbeNoXform ( ProbeInfo&, const Self::CoordinateVectorType& ) const;

  /** Find a voxel in the volume.
   *@param location Real-world coordinates of the location that is to be 
   * found.
   *@param idx This array is used to store the index components of the model
   * voxel containing the given location. Values range from 0 to 
   * ModelDims[dim-1].
   *@param from Real-world coordinate of the voxel's lower left corner.
   *@param to Real-world coordinate of the voxel's upper right corner.
   *@return A non-zero value is returned if and only if the given location
   * lies within the model volume. If zero is returned, the location is outside
   * and all other output values (see parameters) are invalid.
   */
  inline bool FindVoxel( const Self::CoordinateVectorType& location, int *const idx, Types::Coordinate *const from, Types::Coordinate *const to ) const;
  
  /** Find a voxel in the volume.
   *@param location Real-world coordinates of the location that is to be 
   * found.
   *@param idx This array is used to store the index components of the model
   * voxel containing the given location. Values range from 0 to 
   * ModelDims[dim-1].
   *@return A non-zero value is returned if and only if the given location
   * lies within the model volume. If zero is returned, the location is outside
   * and all other output values (see parameters) are invalid.
   */
  inline bool FindVoxel( const Self::CoordinateVectorType& location, int *const idx ) const;

  /** Find a grid index inside or outside the volume.
   *@param location Real-world coordinates of the location that is to be 
   * found.
   *@param idx This array is used to store the index components of the model
   * voxel containing the given location. Values range from 0 to 
   * ModelDims[dim-1].
   */
  inline void GetVoxelIndexNoBounds( const Self::CoordinateVectorType& location, int *const idx ) const;

  /** Find a voxel in the volume by fractional index.
   *@param fracIndex Fractional 3D voxel index.
   *@param idx Output: integer 3D voxel index.
   *@param frac Output: fractional within-voxel location.
   *@return True value is returned if and only if the given location lies
   * within the model volume. If false is returned, the location is outside
   * and all other output values (see parameters) are invalid.
   */
  inline bool FindVoxelByIndex( const Self::CoordinateVectorType& fracIndex, int *const idx, Types::Coordinate *const frac ) const;

  /** Find a voxel in the volume by fractional index without range checking.
   *@param fracIndex Fractional 3D voxel index.
   *@param idx Output: integer 3D voxel index.
   *@param frac Output: fractional within-voxel location.
   */
  inline void FindVoxelByIndexUnsafe( const Self::CoordinateVectorType& fracIndex, int *const idx, Types::Coordinate *const frac ) const;

  /** Find voxel without range checking.
   * This function is identical to FindVoxel except that no range checking is
   * performed. Therefore, it provides additional performance if for every call
   * to this function it is known that the respective voxel actually exists.
   */
  inline void FindVoxelUnsafe( const Self::CoordinateVectorType& location, int *const idx, Types::Coordinate *const from, Types::Coordinate *const to ) const;

  /// Get 3D grid region from continuous lower and upper corner.
  const UniformVolume::RegionType GetGridRange( const Self::CoordinateVectorType& fromVOI, const Self::CoordinateVectorType& toVOI ) const;
 
  /// Get plane coordinate.
  virtual Types::Coordinate GetPlaneCoord( const int axis, const int plane ) const 
  {
    return this->m_Offset[axis] + plane * this->m_Delta[axis];
  }

  /** Get grid index of slice with highest coordinate smaller than given.
   * @param dim specifies the dimension the location parameters refers to.
   * @param location specifies the location in the range from 0 to Size[plane].
   */
  virtual int GetCoordIndex( const int axis, const Types::Coordinate location ) const 
  {
    return std::max<int>( 0, std::min<int>( (int) ((location-this->m_Offset[axis]) / this->m_Delta[axis]), this->m_Dims[axis]-1 ) );
  }
  
  /** Get grid index corresponding (as close as possible) to coordinate.
   * @param dim specifies the dimension the location parameters refers to.
   * @param location specifies the location in the range from 0 to Size[plane].
   */
  virtual int GetClosestCoordIndex( const int axis, const Types::Coordinate location ) const 
  {
    const int idx = (int)MathUtil::Round((location-this->m_Offset[axis]) / this->m_Delta[axis]);
    return std::max<int>( 0, std::min<int>( idx, this->m_Dims[axis]-1 ) );
  }

  /** Get grid index corresponding (as close as possible) to coordinate.
   *\return True if given point is inside image, false if outside.
   */
  virtual bool GetClosestGridPointIndex( const Self::CoordinateVectorType v, int *const xyz ) const 
  {
    for ( int dim = 0; dim < 3; ++dim )
      {
      xyz[dim] = (int) MathUtil::Round((v[dim]-this->m_Offset[dim]) / this->m_Delta[dim]);
      if ( (xyz[dim] < 0) || ( xyz[dim] > this->m_Dims[dim]-1) )
	return false;
      }
    return true;
  }

  /** Get grid index corresponding to coordinate by truncation, not rounding..
   * @param dim specifies the dimension the location parameters refers to.
   * @param location specifies the location in the range from 0 to Size[plane].
   */
  virtual int GetTruncCoordIndex( const int axis, const Types::Coordinate location ) const 
  {
    const int idx = static_cast<int>((location-this->m_Offset[axis]) / this->m_Delta[axis]);
    return std::max<int>( 0, std::min<int>( idx, this->m_Dims[axis]-1 ) );
  }

  /** Get grid index corresponding to coordinate by truncation, not rounding.
   *\return True if given point is inside image, false if outside.
   */
  virtual bool GetTruncGridPointIndex( const Self::CoordinateVectorType v, int *const xyz ) const 
  {
    for ( int dim = 0; dim < 3; ++dim )
      {
      xyz[dim] = static_cast<int>((v[dim]-this->m_Offset[dim]) / this->m_Delta[dim]);
      if ( (xyz[dim] < 0) || ( xyz[dim] > this->m_Dims[dim]-1) )
	return false;
      }
    return true;
  }

  /** Get a grid coordinate.
   * This function directly calculates the grid location from the volume's
   * grid deltas.
   *@param x,y,z The indices of the intended grid element with respect to the
   * three coordinate axes. Valid range is from 0 to Dims[...]-1.
   *@return The location of the given grid element as a Self::CoordinateVectorType.
   */
  virtual const Self::CoordinateVectorType GetGridLocation( const int x, const int y, const int z ) const 
  {
    const Types::Coordinate loc[3] = { this->m_Offset[0] + x * this->m_Delta[0], this->m_Offset[1] + y * this->m_Delta[1], this->m_Offset[2] + z * this->m_Delta[2] };
    return Self::CoordinateVectorType( loc );
  }
  
  /** Get a grid coordinate by continuous pixel index.
   * This function directly calculates the grid location from the volume's
   * grid deltas.
   *@param idx The index of the intended grid element. 
   * Valid range is from 0 to (Dims[0]*Dims[1]*Dims[2])-1.
   *@return The location of the given grid element as a Self::CoordinateVectorType.
   */
  virtual const Self::CoordinateVectorType GetGridLocation( const size_t idx ) const 
  {
    const Types::Coordinate loc[3] = { this->m_Offset[0] +  (idx % this->nextJ) * this->m_Delta[0], 
				       this->m_Offset[1] +  (idx % this->nextK) / this->nextJ * this->m_Delta[1], 
				       this->m_Offset[2] +  (idx / this->nextK) * this->m_Delta[2] };
    return Self::CoordinateVectorType( loc );
  }

  //@}

  /** Set cropping region in real-world coordinates.
   */
  void SetHighResCropRegion( const Self::CoordinateRegionType& crop );

  /** Get cropped volume in real-world coordinates.
   */
  const Self::CoordinateRegionType GetHighResCropRegion() const;

  /// Catch calls to inherited SetCropRegion() and reset high-res crop region.
  virtual void SetCropRegion( const Self::RegionType& region )
  {
    this->m_HighResCropRegion = Self::CoordinateRegionType::SmartPtr( NULL );
    Superclass::SetCropRegion( region );
  }

  /** Calculate volume center.
   *@return Returned is the center of the bounding box.
   */
  Self::CoordinateVectorType GetCenterCropRegion() const 
  {
    const Self::CoordinateRegionType region = this->GetHighResCropRegion();
    return 0.5 * ( region.From() + region.To() );
  }
  
  /** Return cropped uniform volume.
   */
  Self::SmartPtr GetCroppedVolume() const;

  /** Return projection (e.g., MIP, sum) along one axis.
   * This function calls its equivalent in DataGrid and adds calibration
   * info (i.e., correct pixel sizes) to the resulting image.
   *\NOTE The pixel size if taken from the size of the first grid element along
   * each axis -- non-uniform spacings will lead to incorrect results.
   */
  template<class TAccumulator>
  ScalarImage* ComputeProjection( const int axis ) const;
  
  /// Get center of mass of pixel data.
  virtual Self::CoordinateVectorType GetCenterOfMass() const
  {
    Self::CoordinateVectorType com = this->Superclass::GetCenterOfMassGrid();
    for ( int dim = 0; dim < 3; ++dim )
      (com[dim] *= this->m_Delta[dim]) += this->m_Offset[dim];
    return com;
  }
  
  /// Get center of mass of pixel data.
  virtual Self::CoordinateVectorType GetCenterOfMass( Self::CoordinateVectorType& firstOrderMoment ) const
  {
    Self::CoordinateVectorType com = this->Superclass::GetCenterOfMassGrid( firstOrderMoment );
    for ( int dim = 0; dim < 3; ++dim )
      {
      (com[dim] *= this->m_Delta[dim]) += this->m_Offset[dim];
      firstOrderMoment[dim] *= this->m_Delta[dim];
      }
    return com;
  }

  /** Get principal axes.
   *\return This function returns two types of information: the principal directions (axes) of the
   * image data, and the relative scales. The principal axes are returned as a 3x3 array, such that
   * directions[0] is the 3D array of the major principal direction, i.e. the one with the largest
   * norm. Then, directions[0][1] is the y-component of that vector, etc. Accordingly, directions[1]
   * and directions[2] are the two remaining principal directions, with the norm of directions[1]
   * larger than, or equal to, the norm of directions[2].
   */
  void GetPrincipalAxes( Matrix3x3<Types::Coordinate>& directions, Self::CoordinateVectorType& centerOfMass ) const;

protected:
  /** Create a physical copy of this object.
   *@param copyData If true, the associated data array is also copied.
   */
  virtual Self* CloneVirtual( const bool copyData );
  virtual Self* CloneVirtual() const;

  /// Virtual grid cloning constructor.
  virtual Self* CloneGridVirtual() const;

private:
  /** Optional high-resolution crop region.
   * If this is unset (i.e., a NULL pointer), then calls to GetHighResCropRegion() will simply convert the grid-based crop region.
   *
   *\note The crop region includes the volume offset, m_Offset. Its upper limit can, therefore, be larger
   *  than m_Size, which does NOT include m_Offset.
   */
  Self::CoordinateRegionType::SmartPtr m_HighResCropRegion;

  /** Friend declaration of WarpXform class.
   * This allows direct access to Dims and Delta fields in RegisterVolumePoints
   * member function of WarpXform.
   */
  friend class WarpXform;
  
  /// Make axes hash a friend.
  friend class VolumeAxesHash;

  /// Grid grid-to-grid lookup friend.
  friend class VolumeGridToGridLookup;

  class ResampleThreadInfo;

  /// Give thread parameter class access to local types.
  friend class ResampleThreadInfo;

  /** Thread parameter block for volume resampling.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  class ResampleTaskInfo :
    public ThreadParameters<const Self> 
  {
  public:
    /// Array that takes up resamples data.
    Types::DataItem *ResampledData;
    /// Lookup for grid-to-grid cell translation.
    const VolumeGridToGridLookup *GridLookup;
    /// The other volume that has the original data.
    const Self* OtherVolume;
    /// The other volume that has the original data.
    const TypedArray* FromData;
  };

  /// Multi-threaded resampling for grey data.
  static void ResampleThreadPoolExecuteGrey( void *const arg, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCnt );

  /// Multi-threaded resampling for label data (using partial volume averaging).
  static void ResampleThreadPoolExecuteLabels( void *const arg, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCnt );

  /** Extract orthogonal slice as a new volume.
   */
  virtual Self* ExtractSliceVirtual( const int axis, const unsigned int plane ) const;
};

//@}

} // namespace cmtk

#include "cmtkUniformVolume.txx"

#endif // #ifndef __cmtkUniformVolume_h_included_
