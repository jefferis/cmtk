/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkDataGrid_h_included_
#define __cmtkDataGrid_h_included_

#include <cmtkconfig.h>

#include <cmtkTypes.h>
#include <cmtkTypedArray.h>
#include <cmtkVector3D.h>
#include <cmtkScalarImage.h>
#include <cmtkRectangle.h>
#include <cmtkInformationObject.h>
#include <cmtkAnatomicalOrientation.h>

#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Grid topology of data arranged in a 3D lattice.
 * This class extends the plain data handling functions of TypedArray
 * with a 3D topology. Real world coordinates, however, are not considered and
 * need to be handled by derived classes. Thus, this class is used by both
 * Volume and igsVirtualUniformGrid to provide them with coordinate
 * independent services such as median filtering and, to a certain extent,
 * interpolation.
 */
class DataGrid :
  /// Inherit class that handles meta information.
  public InformationObject
{
  /// Number of grid samples in the three spatial dimensions
  igsGetMacro3Array(int,Dims);

  /// Data array (element type is variable)
  igsGetSetMacro(TypedArray::SmartPtr,Data);

public:
  /// This class.
  typedef DataGrid Self;

  /// Smart pointer to DataGrid.
  typedef SmartPointer<Self> SmartPtr;

  /// Default constructor.
  DataGrid() : Data( NULL )
  {
    memset( Dims, 0, sizeof(Dims) );
  }
  
  /// Constructor.
  DataGrid( const int* dims ) : Data( NULL )
  {
    memcpy( Dims, dims, sizeof(Dims) );
  }
  
  /// Virtual destructor.
  virtual ~DataGrid() {}

  /// Downsampling constructor function.
  virtual DataGrid* GetDownsampled( const int (&downsample)[3] ) const;

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
  virtual DataGrid* GetReoriented( const char* newOrientation = AnatomicalOrientation::ORIENTATION_STANDARD ) const;
  
  /** Set dimensions array.
   * This function updates the internal offsets for fast access to adjacent
   * voxel rows, columns, planes etc.
   */
  void SetDims( const int dims[3] );

  /** Create data array.
   *\param dataType ID of scalar data type for the array. This is the image pixel type.
   *\param setToZero If this flag is set, all values in the newly created array will be initialized to zero.
   */
  virtual TypedArray::SmartPtr CreateDataArray( const ScalarDataType dataType, const bool setToZero = false );

  /// Get number of data items in the volume.
  size_t GetNumberOfPixels () const { return Dims[0]*Dims[1]*Dims[2]; }

  /// Check whether given pixel index is inside grid.
  bool IndexIsInRange( const int x, const int y, const int z ) const
  {
    return (x>=0) && (x<Dims[0]) && (y>=0) && (y<Dims[1]) && (z>=0) && (z<Dims[2]);
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
    return Data->Get( data, offset );
  }

  /// Set data at specified offset
  void SetDataAt ( const Types::DataItem data, const size_t offset )
  {
    Data->Set( data, offset );
  }
  
  /// Return data at specified grid point.
  bool GetDataAt ( Types::DataItem& data, const int x, const int y, const int z ) const
  {
    return this->GetDataAt( data, x+Dims[0]*(y+Dims[1]*z) );
  }

  /// Set data at specified grid point.
  void SetDataAt ( const Types::DataItem data, const int x, const int y, const int z )
  {
    this->SetDataAt( data, x+Dims[0]*(y+Dims[1]*z) );
  }

  /// Return data at specified grid point, or a given default value if no data exists there.
  Types::DataItem GetDataAt ( const int x, const int y, const int z, const Types::DataItem defaultValue = 0.0 ) const
  {
    Types::DataItem value;
    if ( this->GetDataAt( value, x+Dims[0]*(y+Dims[1]*z) ) )
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
  TypedArray* GetDataMirrorPlane( const int axis = AXIS_X ) const;

  /// Replace data with mirrored version.
  void ApplyMirrorPlane( const int axis = AXIS_X );

private:
  /// Mirror about plane without allocating additional memory.
  static void MirrorPlaneInPlace( TypedArray *const data, const int dims[3], const int axis = AXIS_X );

public:
  /// Accessor functions for protected member variables
  int GetNextI() const { return nextI; }
  int GetNextJ() const { return nextJ; }
  int GetNextK() const { return nextK; }
  int GetNextIJ() const { return nextIJ; }
  int GetNextIK() const { return nextIK; }
  int GetNextJK() const { return nextJK; }
  int GetNextIJK() const { return nextIJK; }
  
  /// Return data after median-filtering.
  TypedArray* GetDataMedianFiltered( const int radius = 1 ) const;

  /// Replace data with median-filtered version.
  void ApplyMedianFilter( const int radius = 3 );

  /// Return data after median-filtering.
  TypedArray* GetDataSobelFiltered() const;

  /// Replace data with Sobel-filtered version.
  void ApplySobelFilter();

  /// Return data after erosion operator.
  TypedArray* GetDataErode( const int iterations = 1 ) const;

  /// Replace data with eroded data.
  void ApplyErode( const int iterations = 1 );

  /// Return data after dilation operator.
  TypedArray* GetDataDilate( const int iterations = 1 ) const;

  /// Replace data with dilated data.
  void ApplyDilate( const int iterations = 1 );

  /// Return data after eliminating padding data by neighborhood voting.
  TypedArray* GetDataEliminatePaddingVoting() const;

  /** Return data after eliminating padding data by neighborhood voting.
   * This function also returns a flag that is set if and only if
   * a change was made to the data.
   */
  TypedArray* GetDataEliminatePaddingVoting( bool& changed ) const;

  /// Eliminate padding data by neighborhood voting.
  void ApplyEliminatePaddingVoting( const int iterations = 1 );

  /** Return map of region boundaries.
   * This function returns a byte data array where each pixel is one if it is
   * a boundary pixel, i.e., if one of its neighbours in this object has a
   * different value than it has itself. All other pixels are set to zero.
   *\param multiValue If this is set (default: false), then the resulting
   *  boundary map is multi valued, i.e., instead of setting boundary pixels
   *  to "1", they are set to the value present in the image at that location.
   *\note The boundary contours are at least 2 pixels wide since "boundaryness"
   * is a symmetric relationship.
   */
  TypedArray* GetBoundaryMap( const bool multiValued = false ) const;

  /// Get center of mass of pixel data.
  virtual Vector3D GetCenterOfMass() const;
  
  /// Get center of mass and first-order moments of pixel data.
  virtual Vector3D GetCenterOfMass( Vector3D& firstOrderMoment ) const;
  
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
  TypedArray* GetDataMirrored( const int axis = AXIS_X ) const;

  /// Print object.
  void Print() const;

  /** Draw sphere.
    */
  virtual void DrawSphere( const Vector3D& center, const Types::Coordinate radius, const Types::DataItem value );

  /** Draw rectangular box.
   */
  void DrawBox( const IntROI3D& box, const Types::DataItem value );

  /// Return after filtering with a separable kernel
  TypedArray* GetFilteredData( const std::vector<Types::DataItem>& filterX, const std::vector<Types::DataItem>& filterY, const std::vector<Types::DataItem>& filterZ ) const;

private:
  /// Thread parameter for entropy evaluation.
  class FilterThreadParameters : 
    /// Inherit from generic thread parameter class.
    public ThreadParameters<Self>
  {
  public:
    /// Filter kernel.
    const std::vector<Types::DataItem>* m_Filter;

    /// Pointer to result pixel data
    TypedArray* m_Result;
  };
  
  /// Thread function for separable filtering in x-direction.
  static CMTK_THREAD_RETURN_TYPE GetFilteredDataThreadX( void *args );

  /// Thread function for separable filtering in y-direction.
  static CMTK_THREAD_RETURN_TYPE GetFilteredDataThreadY( void *args );

  /// Thread function for separable filtering in z-direction.
  static CMTK_THREAD_RETURN_TYPE GetFilteredDataThreadZ( void *args );

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
  bool TrilinearInterpolation( Types::DataItem&, const int, const int, const int, const Vector3D&, const Types::Coordinate*, const Types::Coordinate* ) const;

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
  ( const TData* dataPtr, const int x, const int y, const int z, const Vector3D& gridPosition, const Types::Coordinate* cellFrom, 
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
};

//@}

} // namespace cmtk

#include <cmtkDataGrid.txx>

#endif // #ifndef __cmtkDataGrid_h_included_
