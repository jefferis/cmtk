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

#ifndef __cmtkAnatomicalOrientation_h_included_
#define __cmtkAnatomicalOrientation_h_included_

#include <cmtkconfig.h>

#include <cmtkSmartPtr.h>
#include <cmtkTypes.h>
#include <cmtkInformationObject.h>
#include <cmtkAffineXform.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class with helper functions for handling anatomical image orientation.
class AnatomicalOrientation
{
public:
  /// This class.
  typedef AnatomicalOrientation Self;

  /// Orientation of STANDARD space (LR/AP/IS).
  static const char *const ORIENTATION_STANDARD;

  /** Get closest anatomical orientation based on space dimensions and image direction vectors.
   *\param orientation The resulting orientation string (three characters plus terminating zero) will
   * be put into this string.
   *\param directions The image direction vectors. "directions[0]" is the 3-D vector of the in-plane x image
   * direction, "directions[1]" is in-plane y, and "directions[2]" is the vector from one plane origin to the
   * origin of the next plane.
   *\param spacing Three-dimensional vector of pixel spacings; directions[i] needs to be divided by spacing[i] to obtain normalized vector (i=0,1,2).
   *\param spaceAxes Six-character string that defines the orientation of the underlying space. The first two characters
   * relate to the x coordinate, the second two to the y coordinate, and the final two to the z coordinate.
   */
  static void GetOrientationFromDirections( char* orientation, const AffineXform::MatrixType& directions, const char* spaceAxes );

  /// Get permutation table of coordinate axes from space axes and image orientation.
  static void GetImageToSpaceAxesPermutation( int (&imageToSpaceAxesPermutation)[3][3], const char* orientation, const char* spaceAxes );

  /** Get closest orientation from a list.
   * This function is used to determine which orientation to bring an image into so it can be written to a file
   * format with limited orientation support (e.g., Analyze).
   */
  static const char* GetClosestOrientation( const char* desiredOrientation, const char *const availableOrientations[] );

  /** Return true if the direction corresponding to the 
   * character 'from' is on the same axis as that corresponding
   * to 'to'.
   *\param from Either L, R, A, P, I, or S
   *\param to Either L, R, A, P, I, or S 
   */
  static bool OnSameAxis( const char from, const char to );

  /// Class for permutation matrix that can be applied to pixel indexes, as well as dimension, and pixel size arrays.
  class PermutationMatrix
  {
  public:
    /// This class.
    typedef PermutationMatrix Self;

    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Constructor: determine axes permutation, flipping, and store local copy of reoriented dimensions array.
    PermutationMatrix( const int* sourceDims, const Types::Coordinate* sourceSize, const std::string curOrientation, const char newOrientation[3] );
   
    /** Takes the dimensions of a volume (e.g. a grid, a voxel)
     * and returns the dimensions of that volume after the
     * reorientation described by this permutation matrix 
     *\param source Original array 
     *\param target A three-element array to contain the
     * permuted results.
     */ 
    template<class T>
    void GetPermutedArray( const T* source, T* target ) const
    {
      for ( int i = 0; i < 3; i++ )
        {
        target[i] = source[this->m_Axes[i]];
        }
    }

    /** Permute index-to-physical matrix
     */
    AffineXform::MatrixType GetPermutedMatrix( const AffineXform::MatrixType& inMatrix ) const;
    
    /** Get new point index from old point index.
     *\param origPoint The input pixel index.
     *\return The offset of the corresponding pixel in the reoriented volume.
     */
    void GetReorientedIndex( const int* origPoint, int* newPoint ) const
    {
      for ( int i = 0; i < 3; ++i )
	newPoint[i] = this->m_Multipliers[i] * origPoint[this->m_Axes[i]] + this->m_Offsets[i];
    }

    /** Applies the permutation described by this matrix to a pixel index.
     *\param origPoint The input pixel index.
     *\return The offset of the corresponding pixel in the reoriented volume.
     */
    size_t NewOffsetFromPoint( const int* origPoint ) const
    {
      return ( this->m_Multipliers[0] * origPoint[this->m_Axes[0]] + this->m_Offsets[0] ) + 
              this->m_NewDims[0] * 
              ( ( this->m_Multipliers[1] * origPoint[this->m_Axes[1]] + this->m_Offsets[1] ) + 
                  this->m_NewDims[1] * 
                  ( this->m_Multipliers[2] * origPoint[this->m_Axes[2]] + this->m_Offsets[2] ) );
    }

  private:
    /** Input-to-output axis index assignments.
     * m_Axes[i] contains the index (0, 1, or 2) of the 
     * axis of the input orientation which gets moved to
     * as the i'th axis of the output orientation.
     * (where i is 0, 1, or 2, corresponding to X, Y, 
     *  and Z, respectively)
     */
    int m_Axes[3];

    /** Multiplies (flip direction) table.
     * m_Multipliers[i] contains -1 if the axis at 
     * m_Axes[i] is to be reversed from its direction
     * prior to re-orientation, and 1 otherwise
     */
    int m_Multipliers[3];
    
    /** Dimension of the reoriented 
     * image in the direction of m_Axes[i]
     */
    int m_NewDims[3];
    
    /** Size of the reoriented 
     * image in the direction of m_Axes[i]
     */
    Types::Coordinate m_NewSize[3];
    
    /** m_Offsets[i] contains 0 if m_Multipliers[i] == 1,
     * and (m_NewDims[i] - 1) if m_Multipliers[i] == -1.
     * (Since in the latter case, we will be writing from
     * the highest-valued pixel in direction i and iterating
     * backwards.)
     */
    int m_Offsets[3];
  };

private:
  /// Get inverse of axis orientation.
  static char OppositeDirection( const char direction )
  {
    const char table[27] = "PbcdefghSjkRmnoAqLItuvwxyz";
    return table[direction-'A'];    
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkAnatomicalOrientation_h_included_
