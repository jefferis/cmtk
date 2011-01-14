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

#ifndef __cmtkAnatomicalOrientation_h_included_
#define __cmtkAnatomicalOrientation_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkAnatomicalOrientationBase.h>

#include <System/cmtkSmartPtr.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkFixedVector.h>
#include <Base/cmtkMetaInformationObject.h>
#include <Base/cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class with helper functions for handling anatomical image orientation.
class AnatomicalOrientation
  : public AnatomicalOrientationBase
{
public:
  /// This class.
  typedef AnatomicalOrientation Self;

  /// Parent class.
  typedef AnatomicalOrientationBase Superclass;

  /** Get closest anatomical orientation based on space dimensions and image direction vectors.
   *\param orientation The resulting orientation string (three characters plus terminating zero) will
   * be put into this string.
   *\param directions The image direction vectors. "directions[0]" is the 3-D vector of the in-plane x image
   * direction, "directions[1]" is in-plane y, and "directions[2]" is the vector from one plane origin to the
   * origin of the next plane.
   *\param spaceAxes Six-character string that defines the orientation of the underlying space. The first two characters
   * relate to the x coordinate, the second two to the y coordinate, and the final two to the z coordinate.
   */
  static void GetOrientationFromDirections( char* orientation, const AffineXform::MatrixType& directions, const char* spaceAxes );

  /// Get permutation table of coordinate axes from space axes and image orientation.
  static void GetImageToSpaceAxesPermutation( int (&imageToSpaceAxesPermutation)[3][3], const char* orientation, const char* spaceAxes );

  /// Class for permutation matrix that can be applied to pixel indexes, as well as dimension, and pixel size arrays.
  class PermutationMatrix
  {
  public:
    /// This class.
    typedef PermutationMatrix Self;

    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Constructor: determine axes permutation, flipping, and store local copy of reoriented dimensions array.
    PermutationMatrix( const FixedVector<3,int>& sourceDims, const std::string& curOrientation, const char newOrientation[3] );
   
    /** Takes the dimensions of a volume (e.g. a grid, a voxel)
     * and returns the dimensions of that volume after the
     * reorientation described by this permutation matrix 
     *\param source Original array 
     *\param target A three-element array to contain the
     * permuted results.
     */ 
    template<class T>
    const FixedVector<3,T> GetPermutedArray( const FixedVector<3,T>& source ) const
    {
      FixedVector<3,T> target;
      for ( int i = 0; i < 3; i++ )
        {
        target[i] = source[this->m_Axes[i]];
        }
      return target;
    }

    /** Permute index-to-physical matrix
     */
    AffineXform::MatrixType GetPermutedMatrix( const AffineXform::MatrixType& inMatrix ) const;
    
    /** Get new point index from old point index.
     *\param origPoint The input pixel index in the original image grid.
     *\param newPoint The output pixel index in the reoriented image grid.
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
    FixedVector<3,int> m_Axes;

    /** Multiplies (flip direction) table.
     * m_Multipliers[i] contains -1 if the axis at 
     * m_Axes[i] is to be reversed from its direction
     * prior to re-orientation, and 1 otherwise
     */
    FixedVector<3,int> m_Multipliers;

    /** Dimension of the reoriented 
     * image in the direction of m_Axes[i]
     */
    FixedVector<3,int> m_NewDims;
    
    /** m_Offsets[i] contains 0 if m_Multipliers[i] == 1,
     * and (m_NewDims[i] - 1) if m_Multipliers[i] == -1.
     * (Since in the latter case, we will be writing from
     * the highest-valued pixel in direction i and iterating
     * backwards.)
     */
    FixedVector<3,int> m_Offsets;
  };
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkAnatomicalOrientation_h_included_
