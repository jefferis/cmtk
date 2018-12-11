/*
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

#ifndef __cmtkUniformDistanceMap_h_included_
#define __cmtkUniformDistanceMap_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkDistanceMap.h>
#include <Base/cmtkUniformVolume.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkThreads.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/** Distance map on a uniform grid.
 *\author This class is based on code originally written by Calvin R. Maurer, Jr.
 */
template<class TDistanceDataType>
class UniformDistanceMap : public DistanceMap
{
public:
  /** This class. */
  typedef UniformDistanceMap<TDistanceDataType> Self;

  /** Superclass. */
  typedef DistanceMap Superclass;

  /// Smart pointer to distance map.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to distance map.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /** Distance data type. */
  typedef TDistanceDataType DistanceDataType;

  /** Construct map from given volume.
   *\param volume 3-D feature image.
   *\param flags Computation flags.
   *\param value Feature value
   *\param window Window radius around feature value.
   */
  UniformDistanceMap( const UniformVolume& volume, const byte flags = Self::DEFAULT, const Types::DataItem value = 0, const Types::DataItem window = 0 );

  // Get the computed distance map.
  UniformVolume::SmartPtr Get()
  {
    return this->m_DistanceMap;
  }

private:
  /// Compute distance map.
  void BuildDistanceMap( const UniformVolume& volume, const byte flags, const Types::DataItem value=0, const Types::DataItem window = 0 );
  
  /// Compute 3-D Euclidean Distance Transformation.
  void ComputeEDT( DistanceDataType *const distance );

  /// Compute 2-D Euclidean Distance Transformation for one image plane.
  void ComputeEDT2D( DistanceDataType *const plane, std::vector<DistanceDataType>& gTemp, std::vector<DistanceDataType>& hTemp );

  /// Compute 1-D Voronoi Euclidean Distance Transform.
  bool VoronoiEDT( DistanceDataType *const lpY, const int nSize, const DistanceDataType delta, std::vector<DistanceDataType>& gTemp, std::vector<DistanceDataType>& hTemp );

  /// Internal: pointer to row storage.
  std::vector< std::vector<DistanceDataType> > m_G;

  /// Internal: pointer to row storage.
  std::vector< std::vector<DistanceDataType> > m_H;

  /** Thread parameters. */
  class ThreadParametersEDT :
    /** Inherit from standard parameters. */
    public ThreadParameters<Self>
  {
  public:
    /** Distance map pointer. */
    DistanceDataType* m_Distance;
  };
  
  /** Thread function for first phase (xy) of EDT computation. */
  static void ComputeEDTThreadPhase1( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /** Thread function for second phase (z) of EDT computation. */
  static void ComputeEDTThreadPhase2( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// The computed distance map.
  UniformVolume::SmartPtr m_DistanceMap;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDistanceMap_h_included_
