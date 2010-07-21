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

#ifndef __cmtkVolume_h_included_
#define __cmtkVolume_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkDataGrid.h"

#include "Base/cmtkMacros.h"
#include "Base/cmtkMathUtil.h"
#include "Base/cmtkFixedVector.h"
#include "Base/cmtkVector3D.h"
#include "Base/cmtkAffineXform.h"
#include "Base/cmtkProbeInfo.h"
#include "Base/cmtkLandmarkList.h"
#include "Base/cmtkAnatomicalOrientation.h"

#include "System/cmtkSmartPtr.h"
#include "System/cmtkThreads.h"
#include "System/cmtkMemory.h"

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** General 3D volume.
 * This class handles three-dimensional volume data with a coordinate
 * transformation and associated distance measure. Methods to retrieve data and
 * general structural information are provided.
 *@author Torsten Rohlfing
 */
class Volume : 
  /// Inherit from 3-D data grid.
  public DataGrid 
{
public:
  /// This class type.
  typedef Volume Self;

  /// Superclass.
  typedef DataGrid Superclass;

  /// Smart pointer to Volume.
  typedef SmartPointer<Self> SmartPtr;
  
  /// Region type.
  typedef Region<3,Types::Coordinate> CoordinateRegionType;

  /// Index type.
  typedef CoordinateRegionType::IndexType CoordinateVectorType;

  /** Volume offset (coordinate of first voxel in RAS standard space).
   *\note This offset is NOT included in the volume size, m_Size.
   */
  CoordinateVectorType m_Offset;

  /// Set volume offset.
  void SetOffset( const Vector3D& o )
  {
    this->m_Offset = o;
  }

  /// List of landmarks defined in this volume.
  LandmarkList::SmartPtr m_LandmarkList;

  /** Spatial extent of the volume in world coordinates
   *\note This is the actual size of the volume between first and last pixel. Therefore,
   * if a non-zero volume coordinate offset is set in m_Offset, this does not affect 
   * this field, because the volume size as the product of pixel size times number of pixels
   * per dimension minus one remains unaffected.
   */
  FixedVector<3,Types::Coordinate> Size;

  /// Default constructor.
  Volume() : m_Offset( CoordinateVectorType::Init( 0.0 ) ) {}

  /** Destructor.
   * Do nothing really; just be present and virtual.
   */
  virtual ~Volume () {};

  /** Get minumum extent.
   *@return Minimum volume extent among the three spatial dimensions.
   */
  virtual Types::Coordinate MinSize () const 
  {
    return std::min<Types::Coordinate>( Size[0], std::min<Types::Coordinate>( Size[1], Size[2] ) );
  }
  
  /** Get maximum extent.
   *@return Maximum volume extent among the three spatial dimensions.
   */
  virtual Types::Coordinate MaxSize () const 
  {
    return std::max<Types::Coordinate>( Size[0], std::max<Types::Coordinate>( Size[1], Size[2] ) );
  }
  
  /** Get total volume.
   *@return Product of the spatial extents in all three coordinate directions.
   */
  virtual Types::Coordinate TotalVolume () const 
  {
    return Size[0] * Size[1] * Size[2];
  }
  
  /// Return average volume of all voxels.
  virtual Types::Coordinate AverageVoxelVolume () const 
  {
    return ( (Size[0]*Size[1]*Size[2]) / ((this->m_Dims[0]-1)*(this->m_Dims[1]-1)*(this->m_Dims[2]-1)) );
  }
  
  /** Calculate volume center.
   *@return Returned is the center of the bounding box.
   */
  Vector3D GetCenter () const;

protected:
  /** Get information needed for trilinear interpolation.
   *@return 1 if operation was successful, 0 if no valid data could be found
   * at the given location.
   */
  bool GetTrilinear ( ProbeInfo&, const int, const int, const int, const Vector3D&, const Types::Coordinate*, const Types::Coordinate* ) const;
};

//@}

} // namespace cmtk

#endif
