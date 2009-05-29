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

#ifndef __cmtkVolume_h_included_
#define __cmtkVolume_h_included_

#include <cmtkconfig.h>

#include <cmtkDataGrid.h>

#include <cmtkMacros.h>
#include <cmtkMathUtil.h>
#include <cmtkVector3D.h>
#include <cmtkAffineXform.h>
#include <cmtkInfinitePlane.h>
#include <cmtkProbeInfo.h>
#include <cmtkLandmarkList.h>
#include <cmtkSmartPtr.h>
#include <cmtkRect3D.h>
#include <cmtkAnatomicalOrientation.h>

#include <cmtkThreads.h>
#include <cmtkMemory.h>

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
  
  /// Volume origin (coordinate of first voxel).
  Vector3D m_Origin;

  /// Set volume origin.
  void SetOrigin( const Vector3D& o )
  {
    this->m_Origin = o;
  }

  /// List of landmarks defined in this volume.
  LandmarkList::SmartPtr m_LandmarkList;

  /// Spatial extent of the volume in world coordinates
  Types::Coordinate Size[3];

  /// Constructor.
  Volume () 
  { 
    Memory::Set<Types::Coordinate>( Size, 0, 3 );
    m_Origin.Set( 0, 0, 0 );
    Memory::Set<int>( CropFrom, 0, 3 );
    Memory::Set<int>( CropTo, 0, 3 );
    Memory::Set<Types::Coordinate>( CropFromReal, 0, 3 );
    Memory::Set<Types::Coordinate>( CropToReal, 0, 3 );
  };
  
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
    return ( (Size[0]*Size[1]*Size[2]) / ((Dims[0]-1)*(Dims[1]-1)*(Dims[2]-1)) );
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

protected:
  /** Start of cropped volume in grid elements.
   */
  int CropFrom[3];

  /** End of cropped volume in grid elements.
   */
  int CropTo[3];

  /** Start of cropped volume in world coordinates.
   */
  Types::Coordinate CropFromReal[3];

  /** End of cropped volume in world coordinates.
   */
  Types::Coordinate CropToReal[3];
};

//@}

} // namespace cmtk

#endif
