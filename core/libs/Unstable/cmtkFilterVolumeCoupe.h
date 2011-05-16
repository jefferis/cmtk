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

#ifndef __cmtkFilterVolume_h_included_
#define __cmtkFilterVolume_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkUnits.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class for filtering volume images.
class FilterVolumeCoupe
{
public:
  /** Apply Coupe denoising filter.
   *\param volume Input 3D image.
   *\param beta Smoothing adjustment parameter
   *\param windowRadius Distance from center voxel to outer edge of window 
   *\return A newly allocated TypedArray object that can be used, for 
   * example, to replace the one held by the input image. The data type of the
   * array is identical to the input array.
   */
  static TypedArray::SmartPtr CoupeFilter
  ( const UniformVolume* volume, 
    const int windowRadius,
    const float beta = 0.5 );

private:
  /** Return the mean value of a vector of Types::DataItems.
   *\param items Block of Types::DataItems. 
   */
  static Types::DataItem Mean( TypedArray::SmartPtr items, const int numItems );

  /** Return the variance of a vector of Types::DataItems.
   *\param items Block of Types::DataItems.
   *\param mean The (pre-computed) mean value of the items vector.
   */
  static Types::DataItem Variance( TypedArray::SmartPtr items, const int numItems, const Types::DataItem mean );

  /** Fill an array with the neighborhood centered around a given voxel.
   * Note that this routine will only include existing
   * voxels.  So asking for a block near the edge
   * of a volume will produce a truncated block.
   *\param block To store the newly-fetched block
   *\param data Data array from input 3D image.
   *\param dims Dimensions of the input 3D image.
   *\param x X-coordinate of block center.
   *\param y Y-coordinate of block center.
   *\param z Z-coordinate of block center.
   */
  static void GetNeighborhood( TypedArray::SmartPtr neighborhood, const int radius, const TypedArray* data, 
			       const int* dims, const int x, const int y, const int z );
  
  /** Fill a CoupeBlock with the block centered around a given voxel.
   * Note that this routine will only include existing
   * voxels.  So asking for a block near the edge
   * of a volume will produce a truncated block.
   *\param block To store the newly-fetched block
   *\param data Data array from input 3D image.
   *\param dims Dimensions of the input 3D image.
   *\param x X-coordinate of block center.
   *\param y Y-coordinate of block center.
   *\param z Z-coordinate of block center.
   */
  static void GetCoupeBlock( TypedArray::SmartPtr block, const TypedArray* data, const int* dims, const int x, const int y, const int z );

  /** Compute the NL-weighting factor between two blocks.
   *\param smoothingParam Smoothing parameter for the weighting function
   *\param centerBlock The block to be restored using NL-means
   *\param outerBlock A block to contribute to the restoration of centerBlock
   */
  static double ComputeCoupeWeight
  ( const Types::DataItem smoothingParam,
    TypedArray::SmartPtr centerBlock, 
    TypedArray::SmartPtr outerBlock );

  /** Add two blocks, in-place on the first block.
   *\param v1 First vector addend, becomes the sum
   *\param v2 Second vector addend
   */
  static void BlockAddInPlace
  ( TypedArray::SmartPtr v1, TypedArray::SmartPtr v2, const int blockSize );

  /** Subtract two blocks.
   *\param diff Block to contain the difference between v1 and v2
   *\param v1 First block subtrahend
   *\param v2 Second block subtrahend
   */
  static void BlockSubtract( TypedArray::SmartPtr diff, TypedArray::SmartPtr v1, TypedArray::SmartPtr v2, const int blockSize );

  /** Return the product of the input vector multiplied by a constant float.
   *\param prod An block to store the product
   *\param items An block containing the items to be multiplied
   *\param mult The constant to multiply by
   */
  static void BlockConstMult( TypedArray::SmartPtr prod, TypedArray::SmartPtr items, const Types::DataItem mult, const int blockSize );
  
  /** Return the squared Euclidean distance between two blocks.
   */
  static double BlockSquaredDistance( TypedArray::SmartPtr centerBlock, TypedArray::SmartPtr outerBlock, const int blockSize );
 
  /** Process a window for NL-means accumulation.
   *\param NL Block in which to put the NL computation
   *\param blockLocations UNDOCUMENTED
   *\param data Data array from input 3D image.
   *\param dims Dimensions of the input 3D image.
   *\param smoothingParam Smoothing parameter for the weighting function
   *\param x X-coordinate of volume center.
   *\param y Y-coordinate of volume center.
   *\param z Z-coordinate of volume center.
   *\param windowRadius Distance from center voxel to outer edge of window 
   *\param beta Smoothing adjustment parameter
   *\param localMeansMap Mean intensity values at each block
   *\param localVariancesMap Variance of the intensities of the center block
   *\param centerBlock The block to be restored using NL-means
   */
   static void ComputeNLWithinWindow
  ( TypedArray::SmartPtr NL,
    const TypedArray* blockLocations,
    const TypedArray* data, const int* dims, const Types::DataItem smoothingParam,
    const int x, const int y, const int z, 
    const int windowRadius, 
    const float beta, 
    const TypedArray* localMeansMap, 
    const TypedArray* localVariancesMap, 
    TypedArray::SmartPtr centerBlock );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFilterVolume_h_included_
