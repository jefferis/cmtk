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

#ifndef __cmtkUniformDistanceMap_h_included_
#define __cmtkUniformDistanceMap_h_included_

#include <cmtkconfig.h>

#include <cmtkUniformVolume.h>
#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

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
class UniformDistanceMap : 
  /// We consider the distance map a volume image as well.
  public UniformVolume
{
public:
  static const long int EDT_MAX_DISTANCE_SQUARED = 2147329548;

  /** This class. */
  typedef UniformDistanceMap<TDistanceDataType> Self;

  /** Superclass. */
  typedef UniformVolume Superclass;

  /// Smart pointer to distance map.
  typedef SmartPointer<Self> SmartPtr;

  /** Distance data type. */
  typedef TDistanceDataType DistanceDataType;

  /** Enumeration with binary flags that control distance map computation.
   * The defined values can be combined by arithmetic "or".
   */
  typedef enum 
  {
    /** No special functions. 
     * This flag will create a distance-from-feature map where any non-zero 
     * voxel in the feature image will be considered.
     */
    DEFAULT = 0,
    /** Compute distance-from-background map.
     * The resulting distance values will be non-zero at all feature voxels,
     * specifying the distance to the nearest non-feature voxel.
     */
    INSIDE = 1,
    /** Use specific feature value.
     * If this flag is set, only voxels in the feature image the value of which
     * matches a given constant will be considered feature voxels. All voxels 
     * with different values, zero or non-zero, will be considered background 
     * voxels.
     */
    VALUE_EXACT = 2,
    /** Use specific feature value.
     * If this flag is set, only voxels in the feature image with values equal
     * to or larger than a given constant will be considered feature voxels.
     * All voxels with lower values will be considered background voxels.
     */
    VALUE_WINDOW = 4,
    /** Use window around specific feature value.
     * If this flag is set, only voxels in the feature image with values that are
     * within a range from a given constant will be considered feature voxels.
     * All voxels with lower values will be considered background voxels.
     */
    VALUE_THRESHOLD = 8
  } Flags;

  /** Construct map from given volume.
   *@param volume 3-D feature image.
   *@param flags Computation flags.
   *@param value Feature value
   *@param window Window radius around feature value.
   */
  UniformDistanceMap( const UniformVolume* volume, const byte flags = DEFAULT, const Types::DataItem value = 0, const Types::DataItem window = 0 );

private:
  /// Compute distance map.
  void BuildDistanceMap( const UniformVolume *volume, const byte flags, const Types::DataItem value=0, const Types::DataItem window = 0 );
  
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
  static CMTK_THREAD_RETURN_TYPE ComputeEDTThreadPhase1( void* args );

  /** Thread function for second phase (z) of EDT computation. */
  static CMTK_THREAD_RETURN_TYPE ComputeEDTThreadPhase2( void* args );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDistanceMap_h_included_
