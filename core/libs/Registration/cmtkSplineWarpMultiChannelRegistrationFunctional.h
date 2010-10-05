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

#ifndef __cmtkSplineWarpMultiChannelRegistrationFunctional_h_included_
#define __cmtkSplineWarpMultiChannelRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkMultiChannelRMIRegistrationFunctional.h>
#include <Registration/cmtkTemplateMultiChannelRegistrationFunctional.h>
#include <Registration/cmtkAffineMultiChannelRegistrationFunctional.h>

#include <Base/cmtkSplineWarpXform.h>

#include <System/cmtkThreads.h>
#include <System/cmtkMutexLock.h>

#include <list>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Forward declaration of class template. */
template<class MetricFunctionalClass> 
class SplineWarpMultiChannelIntensityCorrectionRegistrationFunctional;

/** Class for spline warp multi-channel registration functional. */
template<class TMetricFunctional = MultiChannelRMIRegistrationFunctional<> >
class SplineWarpMultiChannelRegistrationFunctional :
  /** Inherit from multi-channel registration functional base class. */
  public TemplateMultiChannelRegistrationFunctional<SplineWarpXform, TMetricFunctional >
{
public:
  /** This class. */
  typedef SplineWarpMultiChannelRegistrationFunctional Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** This class. */
  typedef TemplateMultiChannelRegistrationFunctional<SplineWarpXform, TMetricFunctional > Superclass;

  /** Metric data subclass */
  typedef typename Superclass::MetricData MetricData;

  /** Default constructor. */
  SplineWarpMultiChannelRegistrationFunctional();

  /** Constructor from affine multi-channel functional. */
  template<class TAffineMetricFunctional>
  SplineWarpMultiChannelRegistrationFunctional
  ( AffineMultiChannelRegistrationFunctional<TAffineMetricFunctional>& affineFunctional );

  /** Set initial affine transformation. */
  virtual void SetInitialAffineTransformation( const AffineXform& initialAffine )
  {
    this->m_InitialAffineTransformation = initialAffine;
  }

  /** Set flag for entropy vs. intensity thresholding. */
  void SetAdaptiveFixEntropyThreshold( const bool flag )
  {
    this->m_AdaptiveFixEntropyThreshold = flag;
  }

  /** Set adaptive fixing entropy threshold factor. */
  void SetAdaptiveFixThreshFactor( const float factor )
  {
    this->m_AdaptiveFixThreshFactor = factor;
  }

  /** Clear list of fixed coordinate dimensions. */
  void ClearFixedCoordinateDimensions()
  {
    this->m_FixedCoordinateDimensions.clear();
  }

  /** Add a coordinate list of fixed coordinate dimensions. */
  void AddFixedCoordinateDimension( const int dim )
  {
    this->m_FixedCoordinateDimensions.push_back( dim );
  }

  /** Set Jacobian volume preservation constraint weight. */
  void SetJacobianConstraintWeight( const float weight )
  {
    this->m_JacobianConstraintWeight = weight;
  }

  /** Initialize transformation. */
  virtual void InitTransformation( const Vector3D& domain, const Types::Coordinate gridSpacing, const bool exact );

  /** Refine transformation control point grid by factor 2. */
  virtual void RefineTransformation();

  /** Reset channels, clear all images. */
  virtual void ClearAllChannels()
  {
    this->ClearReformattedFloatingChannels();
    Superclass::ClearAllChannels();
  }

  /** Compute functional value. */
  virtual typename Self::ReturnType Evaluate();
  
  /** Compute functional value and gradient. */
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 );

protected:
  /** Update all transformation-related data after init, refine, or image change. */
  virtual void NewReferenceChannelGeometry() 
  {
    this->UpdateTransformationData();
  }

private:
  /** Initial affine transformation. */
  AffineXform m_InitialAffineTransformation;

  /** Update all transformation-related data after init, refine, or image change. */
  virtual void UpdateTransformationData();

  /** Floating channels reformatted under current baseline transformation. */
  std::vector< std::vector<float> > m_ReformattedFloatingChannels;
  
  /** Allocate reformatted floating channel memory. */
  virtual void AllocateReformattedFloatingChannels();

  /** Clear and free reformatted floating channel memory. */
  virtual void ClearReformattedFloatingChannels();

  /** Evaluate metric after a local transformation change. */
  typename Self::ReturnType EvaluateIncremental( const SplineWarpXform* warp, MetricData& metricData, const DataGrid::RegionType& region );

  /** Continue metric computation and store reformatted floating channels for local recomputation. */
  virtual void ContinueMetricStoreReformatted( MetricData& metricData, const size_t rindex, const Vector3D& fvector );

  /** Locally undo metric computation. */
  virtual void BacktraceMetric( MetricData& metricData, const DataGrid::RegionType& voi );

  /** Parameter step scale vector. */
  std::vector<Types::Coordinate> m_StepScaleVector;

  /** Table of volumes of influence per warp control point. */
  std::vector<DataGrid::RegionType> m_VolumeOfInfluenceVector;

  /** Flag for entropy vs. intensity thresholding */
  bool m_AdaptiveFixEntropyThreshold;

  /** Threshold for adaptive parameter fixing. */
  float m_AdaptiveFixThreshFactor;

  /** List of coordinate dimensions (x=0, y=1, z=2) that are fixed during optimization. */
  std::list<int> m_FixedCoordinateDimensions;

  /** Weight for Jacobian volume preservation constraint. */
  float m_JacobianConstraintWeight;

  /** Update fixed and active control points in transformation. */
  virtual void UpdateTransformationFixedControlPoints();

  /** Flag whether fixed control points need to be updated. */
  bool m_UpdateTransformationFixedControlPointsRequired;

  /** Number of parallel threads. */
  const size_t m_NumberOfThreads;

  /** Thread function for gradient computation. */
  static void EvaluateThreadFunction( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /** Mutex lock for shared metric data object. */
  MutexLock m_MetricDataMutex;

  /** Separate transformations for threaded gradient computation. */
  std::vector<SplineWarpXform::SmartPtr> m_ThreadTransformations;

  /** Parameters for threaded gradient computation. */
  class EvaluateGradientThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParameters<Self>
  {
  public:
    /// Global parameter step scale factor.
    Types::Coordinate m_Step;

    /// Pointer to array with computed local gradient components.
    Types::Coordinate* m_Gradient;

    /** Current baseline parameter vector. */
    const CoordinateVector* m_ParameterVector;

    /// Current metric value.
    typename Self::ReturnType m_MetricBaseValue;
  };

  /** Thread function for gradient computation. */
  static void EvaluateWithGradientThreadFunction( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /** Make functional class with Jacobian intensity correction a friend. */
  template<class MetricFunctionalClass> friend class SplineWarpMultiChannelIntensityCorrectionRegistrationFunctional;
};

//@}

} // namespace cmtk

#include "cmtkSplineWarpMultiChannelRegistrationFunctional.txx"
#include "cmtkSplineWarpMultiChannelRegistrationFunctionalThreadFunctions.txx"

#endif // #ifndef __cmtkSplineWarpMultiChannelRegistrationFunctional_h_included_
