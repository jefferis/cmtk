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

#ifndef __cmtkMultiChannelHistogramRegistrationFunctional_h_included_
#define __cmtkMultiChannelHistogramRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkMultiChannelRegistrationFunctional.h>

#include <cmtkHashMapSTL.h>
#include <cmtkUniformVolume.h>
#include <cmtkSmartPtr.h>
#include <cmtkUniformVolumeInterpolator.h>
#include <cmtkLinearInterpolator.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/** Base class for multi-channel registration functionals using the Histogram metric. */
template<class TDataType = float, 
	 class TInterpolator = cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear>,
	 class THashKeyType = unsigned int, 
	 char NBitsPerChannel=6>
class MultiChannelHistogramRegistrationFunctional :
  /** Inherit functional interface. */
  public MultiChannelRegistrationFunctional<TInterpolator>
{
public:
  /** Number of bits per channel in the histogram bin index type. */
  static const char m_HistogramBitsPerChannel = NBitsPerChannel;

  /** This class. */
  typedef MultiChannelHistogramRegistrationFunctional<TDataType,THashKeyType> Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** Superclass. */
  typedef MultiChannelRegistrationFunctional<TInterpolator> Superclass;

  /** Real value type for data representation. */
  typedef TDataType DataType;

  /** Hash key type. */
  typedef THashKeyType HashKeyType;

  /** Add reference channel. */
  virtual void AddReferenceChannel( UniformVolume::SmartPtr& channel );

  /** Add floating channel. */
  virtual void AddFloatingChannel( UniformVolume::SmartPtr& channel );

  /** Reset channels, clear all images. */
  virtual void ClearAllChannels();

protected:
  /** Local class for data needed to compute similarity metric. */
  class MetricData
  {
  private:
    /** Typedef of parent class. */
    typedef MultiChannelHistogramRegistrationFunctional<TDataType,TInterpolator,THashKeyType,NBitsPerChannel> Parent;

    /** Parent object. */
    Parent* m_Parent;

  public:
    /** This class type. */
    typedef MetricData Self;

    /** Initialize metric object and local storage. */
    void Init( Parent *const parent );

    /** Hash table type. */
    typedef HashMapSTL<HashKeyType, int> HashTableType;

    /** Joint haash table for reference and floating channels. */
    HashTableType m_JointHash;

    /** Hash table for reference channels. */
    HashTableType m_ReferenceHash;

    /** Hash table for floating channels. */
    HashTableType m_FloatingHash;

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
  virtual Functional::ReturnType GetMetric( const MetricData& metricData ) const;

private:
  /** Scale factors for reference values to histogram "bins". */
  std::vector<TDataType> m_HashKeyScaleRef;

  /** Scale factors for floating values to histogram "bins". */
  std::vector<TDataType> m_HashKeyScaleFlt;

  /** Offsets for reference values to histogram "bins". */
  std::vector<TDataType> m_HashKeyOffsRef;

  /** Offsets for reference values to histogram "bins". */
  std::vector<TDataType> m_HashKeyOffsFlt;

  /** Bit shift count to combine floating and reference hash keys into joint keys. */
  size_t m_HashKeyShiftRef;
};

//@}

} // namespace cmtk

#include <cmtkMultiChannelHistogramRegistrationFunctional.txx>
#include <cmtkMultiChannelHistogramRegistrationFunctionalMetricData.txx>

#endif // #ifndef __cmtkMultiChannelHistogramRegistrationFunctional_h_included_
