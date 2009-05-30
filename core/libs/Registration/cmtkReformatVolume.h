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

#include <math.h>

#ifdef HAVE_VALUES_H
#  include <values.h>
#endif

#ifndef HAVE_STDDEF_H
#  include <stddef.h>
#endif

#include <vector>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>
#include <cmtkInterpolator.h>
#include <cmtkUniformVolumeInterpolator.h>
#include <cmtkTypedArray.h>
#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>
#include <cmtkSplineWarpXform.h>
#include <cmtkMathUtil.h>
#include <cmtkProgress.h>

#include <cmtkMacros.h>
#include <cmtkThreads.h>
#include <cmtkBitVector.h>

#include <cmtkXformList.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Utility class for reformatting volume data.
 * This class takes two volume data sets, and affine and an (optional) local
 * deformation. It provides member functions to reformat from one of the
 * images a slice plane corresponding exactly to one of the original planes
 * from the respective other images. The class provides several different
 * modes, such as standard, overlay, and subtraction. It can also handle 
 * rescaling of data and checkerboard filling of areas with no valid original
 * image data present.
 */
class ReformatVolume 
{
public:
  /// This class.
  typedef ReformatVolume Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Threshold for the reference image.
  igsGetSetMacroDefault(float,LowerThresholdReference,FLT_MIN);

  /// Threshold for the floating image.
  igsGetSetMacroDefault(float,LowerThresholdFloating,FLT_MIN);

  /// Threshold for the reference image.
  igsGetSetMacroDefault(float,UpperThresholdReference,FLT_MAX);

  /// Threshold for the floating image.
  igsGetSetMacroDefault(float,UpperThresholdFloating,FLT_MAX);

  /// Flag for rescaling.
  igsGetSetMacroDefault(bool,Rescale,false);

  /** Flag whether the PADDING value is defined by client.
   */
  igsGetSetMacroDefault(bool,UsePaddingValue,true);

  /** PADDING value as defined by client.
   */
  igsGetSetMacro(Types::DataItem,PaddingValue);

  /// Default constructor.
  ReformatVolume();

  /// Set interpolation mode.
  void SetInterpolation( const cmtk::Interpolators::InterpolationEnum interpolation ) 
  {
    Interpolation = interpolation;
  }

  /// Create interpolator object for given volume according to interpolation mode set in this object.
  UniformVolumeInterpolatorBase* CreateInterpolator( UniformVolume::SmartPtr& volume );

  /// Create interpolator object for given volume according to interpolation mode set in this object.
  static UniformVolumeInterpolatorBase* CreateInterpolator( const cmtk::Interpolators::InterpolationEnum interpolation, UniformVolume::SmartPtr& volume );

  /// Set user-defined data type.
  void SetUserDataType( const ScalarDataType dataType ) 
  {
    this->m_UserDataType = dataType;
  }

  /// Set the reference volume for reformatting.
  void SetReferenceVolume( UniformVolume::SmartPtr& referenceVolume );

  /// Set the floating (transformed) volume for reformatting.
  void SetFloatingVolume( UniformVolume::SmartPtr& floatingVolume );

  /// Set affine transformation to be applied to the floating volume.
  void SetAffineXform( AffineXform::SmartPtr& affineXform );

  /// Set the local deformation to be applied to the reference grid.
  void SetWarpXform( WarpXform::SmartPtr& warpXform );

  /** Set parameters for image data rescaling.
   * As a side effect to setting the scaling parameters, this object is also
   * set into "rescaling mode". No rescaling takes place otherwise.
   * Rescaling is performed AFTER thresholding, if any is done.
   *@param rescaleReference If this parameter is non-zero, rescaling is applied
   * to the reference image. This has no effect for the "Standard" mode of
   * reformatting, as only the floating data is copied to the output in the 
   * mode.
   *@param rescaleOffset The offset value of rescaling. This value is added to
   * the transformed data values after multiplying them with "rescaleSlope".
   *@param rescaleSlope Scaling factor for image data rescaling. The original
   * image data is multiplied by this value before adding "rescaleOffset".
   *@see #UnsetRescale
   */
  void SetRescale( const int rescaleReference, const Types::DataItem rescaleOffset, const Types::DataItem rescaleSlope );

  /// Set flag for checkerboard mode.
  void SetCheckerboardMode( const bool checkerboardMode = true ) 
  {
    CheckerboardMode = checkerboardMode;
  }

  /** Plain reformatting.
   */
  UniformVolume* PlainReformat();

  /** Plain reformatting of a single plane.
   * This function reformats the floating data to a plane spatially identical 
   * to the given plane in the reference image.
   */
  TypedArray* PlainReformat( const int plane, TypedArray *const target = NULL, const size_t targetOffset = 0 );

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

  /// Apply forward warp transformation to reference volume.
  UniformVolume* GetTransformedReference( Types::Coordinate *const volumeOffset = NULL );

  /// Average Jacobians into deformed reference coordinate system.
  UniformVolume* GetTransformedReferenceJacobianAvg
  ( const std::vector<SplineWarpXform::SmartPtr>* xformList, Types::Coordinate *const volumeOffset = NULL, const bool includeReferenceData = true );
  
  /// Push-forward reformat.
  static TypedArray* ReformatPushForward( const UniformVolume* floating, cmtk::XformList& targetToRef, const UniformVolume* reference );
  
  /// Push-forward reformat with value accumulation.
  static TypedArray* ReformatPushForwardAccumulate( const UniformVolume* floating, cmtk::XformList& targetToRef, const UniformVolume* reference );
  
  /// Complex reformat.
  template<class TInterpolator, class Fct> static TypedArray* Reformat
  ( const UniformVolume* target, cmtk::XformList& targetToRef, const UniformVolume* reference, cmtk::XformList& refToFloat, TInterpolator& interpolator, Fct& fct );
  
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
    bool operator()( const Vector3D& inRef, cmtk::XformList& refToFloat, TInterpolatorInstantiationPtr& interpolator, Types::DataItem& value );
    
    /// Return preferred data type for reformatted data.
    ScalarDataType GetDataType( const UniformVolume*, const UniformVolume* flt ) const
    { 
      if ( this->DataType != TYPE_NONE )
	return this->DataType;
      else
	return flt->GetData()->GetType(); 
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
    bool operator()( const Vector3D& inRef, cmtk::XformList& refToFloat, TInterpolatorInstantiationPtr& interpolator, Types::DataItem& value );
    
    /** Return preferred data type for reformatted data. */
    ScalarDataType GetDataType( const UniformVolume*, const UniformVolume* ) const
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
  UniformVolume::SmartPtr ReferenceVolume;

  /// Pointer to the floating volume.
  UniformVolume::SmartPtr FloatingVolume;

  /// Make target volume matching (reoriented) reference volume.
  UniformVolume* MakeTargetVolume() const;

  /// Pointer to the affine transformation of the floating volume.
  AffineXform::SmartPtr AffineXform;
  
  /// Pointer to the local deformation of the reference grid.
  WarpXform::SmartPtr WarpXform;

  /** Maximum value in reference and floating image.
   * This field is updated by both SetXXXVolume methods. Its value is used for
   * the optional generation of a checkerboard pattern.
   */
  Types::DataItem MaximumValue;

  /** Rescale reference flag.
   * If this flag is non-zero, rescaling is applied to the data from the 
   * reference image. Otherwise, the data from the floating image is rescaled.
   */
  int RescaleReference;

  /** Offset for image data rescaling.
   * This value is added to the original data values after multiplying them 
   * with "RescaleSlope".
   */
  Types::DataItem RescaleOffset;

  /** Factor for image data rescaling.
   * The original data values are multiplied with this value first during
   * rescaling.
   */
  Types::DataItem RescaleSlope;

  /// Flag for checkerboard mode.
  bool CheckerboardMode;

  class GetTransformedReferenceTP : 
    public ThreadParameters<ReformatVolume>
  {
  public:
    GetTransformedReferenceTP() : m_Offset( 0 ), m_Stride( 1 ) {};
    TypedArray::SmartPtr dataArray;
    const SplineWarpXform* splineXform;
    const int* dims;
    /// Offset for blockwise distributed computation.
    size_t m_Offset;
    size_t m_Stride;
    const Types::Coordinate* delta;
    const Types::Coordinate* bbFrom;
    unsigned int numberOfImages;
    const std::vector<SplineWarpXform::SmartPtr>* xformList;
    const std::vector<const UniformVolumeInterpolatorBase*>* interpolatorList;
    const UniformVolumeInterpolatorBase* referenceInterpolator;
    int maxLabel;
    AveragingMode avgMode;
    bool IncludeReferenceData;
  };

  /// Apply forward warp transformation to grey-level reference volume.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceGrey( void *const arg );

  /// Apply forward warp transformation to average grey-level reference volume.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceGreyAvg( void *const arg );

  class GetTransformedReferenceLabelTP : 
    public ThreadParameters<ReformatVolume>
  {
  public:
    GetTransformedReferenceLabelTP() : m_Offset( 0 ), m_Stride( 1 ) {};
    TypedArray::SmartPtr dataArray;
    const SplineWarpXform* splineXform;
    const int* dims;
    /// Offset for blockwise distributed computation.
    size_t m_Offset;
    size_t m_Stride;
    const Types::Coordinate* delta;
    const Types::Coordinate* bbFrom;
    unsigned int numberOfImages;
    const std::vector<SplineWarpXform::SmartPtr>* xformList;
    const std::vector<UniformVolume::SmartPtr>* volumeList;
    int maxLabel;
    AveragingMode avgMode;
    bool IncludeReferenceData;
  };

  /// Apply forward warp transformation to label reference volume.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceLabel( void *const arg );

  /// Thread function: average Jacobians into deformed reference system.
  static CMTK_THREAD_RETURN_TYPE GetTransformedReferenceJacobianAvgThread( void *const arg );

  /// Create uniform volume with correct dimensions for reformatted reference.
  UniformVolume* CreateTransformedReference( Types::Coordinate *const bbFrom, Types::Coordinate *const delta, Types::Coordinate *const volumeOffset = NULL );

};

//@}

} // namespace cmtk

#include <cmtkReformatVolumeReformat.txx>

#endif // #ifndef __cmtkReformatVolume_h_included_
