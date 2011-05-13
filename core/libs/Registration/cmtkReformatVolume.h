/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkReformatVolume_h_included_
#define __cmtkReformatVolume_h_included_

#include <cmtkconfig.h>

#ifdef HAVE_VALUES_H
#  include <values.h>
#endif

#include <Base/cmtkVolume.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkInterpolator.h>
#include <Base/cmtkUniformVolumeInterpolator.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkMathUtil.h>
#include <Base/cmtkMacros.h>
#include <Base/cmtkBitVector.h>
#include <Base/cmtkXformList.h>

#include <System/cmtkThreads.h>
#include <System/cmtkProgress.h>

#include <cstddef>
#include <math.h>
#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Utility class for reformatting volume data.
 * This class takes two volume data sets, and affine and an (optional) local
 * deformation. It provides member functions to reformat from one of the
 * images a slice plane corresponding exactly to one of the original planes
 * from the respective other images. The class provides optional checkerboard filling 
 * of areas with no valid original image data present.
 */
class ReformatVolume 
{
public:
  /// This class.
  typedef ReformatVolume Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Threshold for the reference image.
  cmtkGetSetMacroDefault(float,LowerThresholdReference,FLT_MIN);

  /// Threshold for the floating image.
  cmtkGetSetMacroDefault(float,LowerThresholdFloating,FLT_MIN);

  /// Threshold for the reference image.
  cmtkGetSetMacroDefault(float,UpperThresholdReference,FLT_MAX);

  /// Threshold for the floating image.
  cmtkGetSetMacroDefault(float,UpperThresholdFloating,FLT_MAX);

  /** Flag whether the PADDING value is defined by client.
   */
  cmtkGetSetMacroDefault(bool,UsePaddingValue,true);

  /** PADDING value as defined by client.
   */
  cmtkGetSetMacro(Types::DataItem,PaddingValue);

  /// Default constructor.
  ReformatVolume();

  /// Set interpolation mode.
  void SetInterpolation( const cmtk::Interpolators::InterpolationEnum interpolation ) 
  {
    Interpolation = interpolation;
  }

  /// Create interpolator object for given volume according to interpolation mode set in this object.
  UniformVolumeInterpolatorBase::SmartPtr CreateInterpolator( const UniformVolume::SmartConstPtr& volume );

  /// Create interpolator object for given volume according to interpolation mode set in this object.
  static UniformVolumeInterpolatorBase::SmartPtr CreateInterpolator( const cmtk::Interpolators::InterpolationEnum interpolation, const UniformVolume::SmartConstPtr& volume );

  /// Set user-defined data type.
  void SetUserDataType( const ScalarDataType dataType ) 
  {
    this->m_UserDataType = dataType;
  }

  /// Set the reference volume for reformatting.
  void SetReferenceVolume( const UniformVolume::SmartConstPtr& referenceVolume );

  /// Set the floating (transformed) volume for reformatting.
  void SetFloatingVolume( const UniformVolume::SmartConstPtr& floatingVolume );

  /// Set affine transformation to be applied to the floating volume.
  void SetAffineXform( const AffineXform::SmartPtr& affineXform );

  /// Set the local deformation to be applied to the reference grid.
  void SetWarpXform( const WarpXform::SmartPtr& warpXform );

  /** Plain reformatting.
   */
  const UniformVolume::SmartPtr PlainReformat();

  /** Plain reformatting of a single plane.
   * This function reformats the floating data to a plane spatially identical 
   * to the given plane in the reference image. This is useful for interactive
   * reformatting, where we want a single plane reformatted as fast as possible.
   */
  TypedArray::SmartPtr PlainReformat( const int plane, TypedArray::SmartPtr& target = TypedArray::SmartPtr::Null(), const size_t targetOffset = 0 );
  
  /// Constants for transformation field mode.
  typedef enum 
  {
    MODE_MEAN,
    MODE_MEDIAN,
    MODE_ROBUST90
  } AveragingMode;
  
  /// Apply forward warp transformation to reference volume.
  UniformVolume* GetTransformedReference
  ( const std::vector<SplineWarpXform::SmartPtr>* xformList, std::vector<UniformVolume::SmartPtr>* volumeList, 
    Types::Coordinate *const volumeOffset = NULL, const bool includeReferenceData = true );

  /// Average Jacobians into deformed reference coordinate system.
  UniformVolume* GetTransformedReferenceJacobianAvg
  ( const std::vector<SplineWarpXform::SmartPtr>* xformList, Types::Coordinate *const volumeOffset = NULL, const bool includeReferenceData = true );
  
  /// Complex reformat using target data as the mask.
  template<class TInterpolator, class Fct> static TypedArray::SmartPtr ReformatMasked
  ( const UniformVolume* target, const cmtk::XformList& targetToRef, const cmtk::XformList& refToFloat, Fct& fct, const UniformVolume* floating = NULL, TInterpolator& interpolator = TInterpolator::Null );
  
  /// Complex reformat without mask.
  template<class TInterpolator, class Fct> static TypedArray::SmartPtr ReformatUnmasked
  ( const UniformVolume* target, const cmtk::XformList& targetToRef, const cmtk::XformList& refToFloat, Fct& fct, const UniformVolume* floating = NULL, TInterpolator& interpolator = TInterpolator::Null );
  
  /// Constants for extended reformatting mode.
  typedef enum 
  {
    REFORMAT_PLAIN,
    REFORMAT_JACOBIAN
  } Mode;
  
  /// Function class for plain reformating.
  class Plain
  {
  public:
    /** Constructor. */
    Plain( const ScalarDataType dataType = TYPE_NONE ) 
      : DataType( dataType ),
	PaddingValue( 0 ),
	UsePaddingValue( false )
      {};
    
    /// Set output padding value.
    void SetPaddingValue( const Types::DataItem paddingValue )
    {
      this->PaddingValue = paddingValue;
      this->UsePaddingValue = true;
    }

    /** Query operator. */
    template <class TInterpolatorInstantiationPtr>
    bool operator()( Types::DataItem& value, const Vector3D& inRef, const cmtk::XformList& refToFloat, TInterpolatorInstantiationPtr& interpolator );
    
    /// Return preferred data type for reformatted data.
    ScalarDataType GetDataType( const UniformVolume& fltImage ) const
    { 
      if ( this->DataType != TYPE_NONE )
	return this->DataType;
      else
	return fltImage.GetData()->GetType(); 
    }

  protected:
    /** Data type for reformated image.
     */
    ScalarDataType DataType;

  public:
    /** Padding value for output. */
    Types::DataItem PaddingValue;

    /** Use padding value. */
    bool UsePaddingValue;
  };
  
  /// Function class for reformating of a Jacobian map.
  class Jacobian :
    /// Inherit from plain reformatter for interpolation etc.
    public ReformatVolume::Plain
  {
  public:
    /** Constructor. */
    Jacobian( const ScalarDataType dataType = TYPE_FLOAT, const bool correctGlobalScale = true ) 
      : Plain( dataType ), CorrectGlobalScale( correctGlobalScale ) {};
    
    /** Query operator. */
    template <class TInterpolatorInstantiationPtr>
    bool operator()( Types::DataItem& value, const Vector3D& inRef, const cmtk::XformList& refToFloat, TInterpolatorInstantiationPtr& interpolator );
    
    /** Return preferred data type for reformatted data. */
    ScalarDataType GetDataType( const UniformVolume& ) const
    { 
      if ( this->DataType != TYPE_NONE )
	return this->DataType;
      else
	return TYPE_FLOAT;
    }
    
  private:
    /** Flag for correction of global scale.
      * If this is true, all Jacobian determinants are divided (forward 
      * transformation) or multiplied (inverse transformation) by the
      * global affine scale factor of the transformation.
      */
    bool CorrectGlobalScale;
  };

private:
  /// Interpolation mode.
  cmtk::Interpolators::InterpolationEnum Interpolation;

  /// User-selected data type for reformatted data.
  ScalarDataType m_UserDataType;

  /// Pointer to the reference volume.
  UniformVolume::SmartConstPtr ReferenceVolume;

  /// Pointer to the floating volume.
  UniformVolume::SmartConstPtr FloatingVolume;

  /// Make target volume matching (reoriented) reference volume.
  const UniformVolume::SmartPtr MakeTargetVolume() const;

  /// Pointer to the affine transformation of the floating volume.
  AffineXform::SmartConstPtr m_AffineXform;
  
  /// Pointer to the local deformation of the reference grid.
  WarpXform::SmartConstPtr m_WarpXform;

  class GetTransformedReferenceTP : 
    public ThreadParameters<ReformatVolume>
  {
  public:
    GetTransformedReferenceTP() : m_Offset( 0 ), m_Stride( 1 ) {};
    TypedArray::SmartPtr dataArray;
    const SplineWarpXform* splineXform;
    DataGrid::IndexType dims;
    /// Offset for blockwise distributed computation.
    size_t m_Offset;
    size_t m_Stride;
    const Types::Coordinate* delta;
    const Types::Coordinate* bbFrom;
    unsigned int numberOfImages;
    const std::vector<SplineWarpXform::SmartPtr>* xformList;
    const std::vector<UniformVolume::SmartPtr>* volumeList;
    const std::vector<UniformVolumeInterpolatorBase::SmartConstPtr>* interpolatorList;
    const UniformVolumeInterpolatorBase* referenceInterpolator;
    int maxLabel;
    AveragingMode avgMode;
    bool IncludeReferenceData;
  };

  /// Apply forward warp transformation to grey-level reference volume.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceGrey( void *const arg );

  /// Apply forward warp transformation to average grey-level reference volume.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceGreyAvg( void *const arg );

  /// Apply forward warp transformation to label reference volume.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceLabel( void *const arg );

  /// Thread function: average Jacobians into deformed reference system.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceJacobianAvgThread( void *const arg );

  /// Create uniform volume with correct dimensions for reformatted reference.
  UniformVolume* CreateTransformedReference( Types::Coordinate *const bbFrom, Types::Coordinate *const delta, Types::Coordinate *const volumeOffset = NULL );

};

//@}

} // namespace cmtk

#include "cmtkReformatVolumeReformat.txx"
#include "cmtkReformatVolumePlain.txx"
#include "cmtkReformatVolumeJacobian.txx"

#endif // #ifndef __cmtkReformatVolume_h_included_
