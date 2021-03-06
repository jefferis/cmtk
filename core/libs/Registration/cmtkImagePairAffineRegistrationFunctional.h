/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkImagePairAffineRegistrationFunctional_h_included_
#define __cmtkImagePairAffineRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkImagePairRegistrationFunctional.h>

#include <Base/cmtkVector.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkVolume.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkMathUtil.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkVolumeClipping.h>

#include <System/cmtkException.h>
#include <System/cmtkThreadPool.h>

#include <cassert>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base-class for affine registration functionals.
 */
class ImagePairAffineRegistrationFunctional : 
  /// Inherit from voxel matching functional.
  public ImagePairRegistrationFunctional
{
public:
  /// This class type.
  typedef ImagePairAffineRegistrationFunctional Self;

  /// Superclass.
  typedef ImagePairRegistrationFunctional Superclass;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v )  
  {
    this->m_AffineXform->GetParamVector( v );
  }

  /// Return parameter stepping.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const;

  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const 
  {
    return this->m_AffineXform->ParamVectorDim();
  }

  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const 
  {
    return this->m_AffineXform->VariableParamVectorDim();
  }

  /// Set optional restriction to axis-orthogonal in-plane transformations.
  void SetRestrictToInPlane( const int axis )
  {
    this->m_RestrictToInPlane = axis;
  }

protected:
  /// Current coordinate transformation.
  AffineXform::SmartPtr m_AffineXform;

  /** Restrict to in-plance transformations.
   * If this is -1 (default), then the full 3D transformation is optimized. If set to 0, 1, or 2, the
   * value is the index of a coordinate axis (x, y, or z), and the transformation is restricted to
   * components in-plane perpendicular to that axis, e.g., rotation around the axis and translations
   * perpendicular to it.
   */
  int m_RestrictToInPlane;

  /// Utility object for volume clipping.
  VolumeClipping Clipper;

  /** Perform clipping/cropping in z-direction.
   * This function computes the intersection of reference and floating data in
   * z-direction. It determines the range of indices of those planes in the
   * reference that intersect the floating. This is the range over which to 
   * for-loop during metric computation.
   *\param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *\param origin Starting point of the reference volume.
   *\param start Upon return, this reference is set to the index of first plane
   * in the reference that intersects the floating.
   *\param end Upon return, this reference is set to one plus the index of the
   * last plane in the reference that intersects the floating.
   *\return 1 if there is an intersection of reference and floating, 0 if there
   * isn't. The range of indices returned in "start" and "end" is only
   * guaranteed to be valid if 1 is the return value.
   */
  Types::GridIndexType ClipZ ( const VolumeClipping& clipper, const Vector3D& origin, DataGrid::IndexType::ValueType& start, DataGrid::IndexType::ValueType &end ) const
  {
    // perform clipping
    Types::Coordinate fromFactor, toFactor;
    if (! clipper.ClipZ( fromFactor, toFactor, origin ) )
      return 0;

    // there is an intersection: Look up the corresponding grid indices
    start = static_cast<DataGrid::IndexType::ValueType>( (this->m_ReferenceDims[2]-1)*fromFactor );
    end = 1+std::min( (Types::GridIndexType)(this->m_ReferenceDims[2]-1), (Types::GridIndexType)(1 + ((this->m_ReferenceDims[2]-1)*toFactor) ) );
    
    // finally, apply cropping boundaries of the reference volume
    start = std::max<DataGrid::IndexType::ValueType>( start, this->m_ReferenceCropRegion.From()[2] );
    end = std::min<DataGrid::IndexType::ValueType>( end, this->m_ReferenceCropRegion.To()[2] );
    
    // return 1 iff index range is non-empty.
    return (start < end );
  }

  /** Perform clipping/cropping in x-direction.
   * This function computes the intersection of reference and floating data in
   * x-direction. It determines the range of indices of those voxels in the
   * current reference row that intersect the floating image. This is the range
   * over which to for-loop during metric computation.
   *
   * Compared to ClipZ and ClipY, this step has to operate very exact as there
   * is no further level that would reduce remaining invalid voxels. Therefore,
   * clipper.ClipX() is called with an extended initial range of indices and an
   * explicitly open upper bound.
   *
   * This is necessary to discriminate inside-boundary from on-boundary voxels.
   * For the right, upper and back boundary, on-boundary voxels are already
   * outside the allowed range as the upper boundaries of the volume are open
   * in terms of interpolation.
   *\param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *\param origin Starting point of the current row in the reference volume.
   *\param start Upon return, this reference is set to the index of first voxel
   * in the reference that intersects the floating image.
   *\param end Upon return, this reference is set to one plus the index of the
   * last voxel in the reference that intersects the floating image.
   *\return 1 if there is an intersection of the current reference row and
   * the floating, 0 if there isn't. The range of indices returned in "start"
   * and "end" is only guaranteed to be valid if 1 is the return value.
   */
  Types::GridIndexType ClipX ( const VolumeClipping& clipper, const Vector3D& origin, DataGrid::IndexType::ValueType& start, DataGrid::IndexType::ValueType &end ) const
  {
    // perform clipping
    Types::Coordinate fromFactor, toFactor;
    if ( ! clipper.ClipX( fromFactor, toFactor, origin, 0, 2, false, true ) )
      return 0;

    fromFactor = std::min<Types::Coordinate>( 1.0, fromFactor );
	      
    // there is an intersection: Look up the corresponding grid indices
    start = std::max<Types::GridIndexType>( 0, (Types::GridIndexType)((this->m_ReferenceDims[0]-1)*fromFactor)-1 );
    while ( ( start*this->m_ReferenceGrid->m_Delta[0] < fromFactor*this->m_ReferenceSize[0]) && ( start < this->m_ReferenceDims[0] ) ) 
      ++start;
    
    if ( (toFactor > 1.0) || (start == this->m_ReferenceDims[0]) ) 
      {
      end = this->m_ReferenceDims[0];
      } 
    else
      {
      end = std::min( this->m_ReferenceDims[0]-2, (Types::GridIndexType)(1 + (this->m_ReferenceDims[0]-1)*toFactor));
      while ( end*this->m_ReferenceGrid->m_Delta[0] > toFactor*this->m_ReferenceSize[0] ) // 'if' not sufficient!	
	--end;
      ++end; // otherwise end=1+min(...) and ...[0][end-1] above!!
      }
    
    // finally, apply cropping boundaries of the reference volume
    start = std::max<DataGrid::IndexType::ValueType>( start, this->m_ReferenceCropRegion.From()[0] );
    end = std::min<DataGrid::IndexType::ValueType>( end, this->m_ReferenceCropRegion.To()[0] );
    
    // return 1 iff index range is non-empty.
    return (start < end );
  }

  /** Perform clipping/cropping in y-direction.
   * This function computes the intersection of reference and floating data in
   * y-direction. It determines the range of indices of those rows in the
   * current reference plane that intersect the floating image. This is the
   * range over which to for-loop during metric computation.
   *\param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *\param origin Starting point of the current plane in the reference volume.
   *\param start Upon return, this reference is set to the index of first row
   * in the reference that intersects the floating image.
   *\param end Upon return, this reference is set to one plus the index of the
   * last row in the reference that intersects the floating image.
   *\return 1 if there is an intersection of the current reference plane and
   * the floating, 0 if there isn't. The range of indices returned in "start" 
   * and "end" is only guaranteed to be valid if 1 is the return value.
   */
  Types::GridIndexType ClipY ( const VolumeClipping& clipper, const Vector3D& origin, DataGrid::IndexType::ValueType& start, DataGrid::IndexType::ValueType &end ) const
  {
    // perform clipping
    Types::Coordinate fromFactor, toFactor;
    if ( !clipper.ClipY( fromFactor, toFactor, origin ) )
      return 0;

    // there is an intersection: Look up the corresponding grid indices
    start = static_cast<DataGrid::IndexType::ValueType>( (this->m_ReferenceDims[1]-1)*fromFactor );
    
    if ( toFactor > 1.0 ) 
      {
      end = this->m_ReferenceDims[1];
      } 
    else
      {
      end = 1+std::min( this->m_ReferenceDims[1]-1, (Types::GridIndexType)(1+(this->m_ReferenceDims[1]-1)*toFactor ) );
      }
    // finally, apply cropping boundaries of the reference volume
    start = std::max<DataGrid::IndexType::ValueType>( start, this->m_ReferenceCropRegion.From()[1] );
    end = std::min<DataGrid::IndexType::ValueType>( end, this->m_ReferenceCropRegion.To()[1] );
    
    // return 1 iff index range is non-empty.
    return (start < end );
  }

public:
  /// Constructor.
  ImagePairAffineRegistrationFunctional( UniformVolume::SmartPtr refVolume, UniformVolume::SmartPtr modVolume, AffineXform::SmartPtr& affineXform ) 
    : ImagePairRegistrationFunctional( refVolume, modVolume ),
      m_AffineXform( affineXform ),
      m_RestrictToInPlane( -1 )
  {}

  /// Destructor.
  virtual ~ImagePairAffineRegistrationFunctional() {}

  /** Constructor function for affine voxel registration functionals.
   * This function takes the index of a metric in the list of available voxel
   * similarity measures plus all required objects. It the creates an appropriate
   * instance of ImagePairAffineRegistrationFunctional with the correct metric class as template
   * parameter.
   */
  static ImagePairAffineRegistrationFunctional* Create( const int metric /*!< Index of image similarity measure.*/,
							UniformVolume::SmartPtr& refVolume /*!< Reference volume.*/,
							UniformVolume::SmartPtr& fltVolume /*!< Floating volume*/,
							const Interpolators::InterpolationEnum interpolation /*!< Floating volume intnerpolation.*/,
							AffineXform::SmartPtr& affineXform /*!< Use this affine transformation.*/ );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairAffineRegistrationFunctional_h_included_
