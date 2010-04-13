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

#ifndef __cmtkMultiChannelRMIRegistrationFunctional_h_included_
#define __cmtkMultiChannelRMIRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkMultiChannelRegistrationFunctional.h>

#include <cmtkUniformVolume.h>
#include <cmtkSmartPtr.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/** Base class for multi-channel registration functionals using the RMI metric. */
template<class TRealType = double,
	 class TDataType = float,
	 class TInterpolator = UniformVolumeInterpolator<Interpolators::Linear> >
class MultiChannelRMIRegistrationFunctional :
  /** Inherit functional interface. */
  public MultiChannelRegistrationFunctional<TInterpolator>
{
public:
  /** This class. */
  typedef MultiChannelRMIRegistrationFunctional<TRealType,TDataType> Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** Superclass. */
  typedef MultiChannelRegistrationFunctional<TInterpolator> Superclass;

  /** Real value type for internal computations. */
  typedef TRealType RealType;

  /** Real value type for data representation. */
  typedef TDataType DataType;

protected:
  /** Local class for data needed to compute similarity metric. */
  class MetricData
  {
  private:
    /** Typedef of parent class. */
    typedef MultiChannelRMIRegistrationFunctional<TRealType,TDataType> Parent;
    
  public:
    /** This class type. */
    typedef MetricData Self;

    /** Initialize metric object and local storage. */
    void Init( Parent *const parent );

    /** Vector of pixel value sums. */
    std::vector<RealType> m_Sums;
    
    /** Vector (actually matrix) of pairwise pixel value products. */
    std::vector<RealType> m_Products;
    
    /** Covariance matrix for joint entropy computation. */
    Matrix2D<RealType> m_CovarianceMatrix;
    
    /** Covariance matrix for reference channels entropy computation. */
    Matrix2D<RealType> m_CovarianceMatrixRef;
    
    /** Covariance matrix for floating channels entropy computation. */
    Matrix2D<RealType> m_CovarianceMatrixFlt;
    
    /** Total number of samples (pixels) under current transformation. */
    size_t m_TotalNumberOfSamples;

    /** Assignment operator. */
    Self& operator=( const Self& source );

    /** In-place addition operator. */
    Self& operator+=( const Self& other );

    /** In-place subtraction operator. */
    Self& operator-=( const Self& other );

    /** In-place single sample addition operator. */
    void operator+=( const Types::DataItem* values );

    /** In-place single sample subtraction operator. */
    void operator-=( const Types::DataItem* values );
  };

  /// Global data structure for metric computation.
  MetricData m_MetricData;

  /** Continue metric computation. */
  virtual void ContinueMetric( MetricData& metricData, const size_t rindex, const Vector3D& fvector );

  /** Get metric value. */
  virtual RealType GetMetric( const MetricData& metricData ) const;
};

//@}

} // namespace cmtk

#include <cmtkMultiChannelRMIRegistrationFunctional.txx>
#include <cmtkMultiChannelRMIRegistrationFunctionalMetricData.txx>

#endif // #ifndef __cmtkMultiChannelRMIRegistrationFunctional_h_included_
