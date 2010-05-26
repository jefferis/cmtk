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

#ifndef __cmtkAffineGroupwiseRegistrationRMIFunctional_h_included_
#define __cmtkAffineGroupwiseRegistrationRMIFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkGroupwiseRegistrationRMIFunctional.h>

#include <cmtkSmartPtr.h>

#include <cmtkUniformVolume.h>
#include <cmtkAffineXform.h>
#include <cmtkHistogram.h>

#include <vector>

#include <cmtkClassStream.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for affine groupwise registration using "RMI" metric..
 */
class AffineGroupwiseRegistrationRMIFunctional : 
  /// Inherit general RMI-based groupwise registration.
  public GroupwiseRegistrationRMIFunctional<AffineXform>
{
public:
  /// Type of parent class.
  typedef GroupwiseRegistrationRMIFunctional<AffineXform> Superclass;

  /// Type of this class.
  typedef AffineGroupwiseRegistrationRMIFunctional Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  AffineGroupwiseRegistrationRMIFunctional();

  /// Destructor.
  virtual ~AffineGroupwiseRegistrationRMIFunctional();

  /// Set number of degrees of freedom per transformation.
  void SetXformNumberDOFs( const int numberDOFs );

  /** Set affine transformations.
   */
  void SetXforms( const std::vector<AffineXform::SmartPtr>& xformVector );

protected:
  /// Number of DOFs per transformation.
  int m_XformNumberDOFs;

  /** Interpolate given moving image to template.
   * This function overrides the interpolation function provided by the base
   * class. It makes use of the fact that affine transformations preserve
   * parallel lines for more efficient computation.
   *\param idx Index of of to reformat to template. This also determines which
   *  transformation is used.
   *\param destination The reformatted pixel data is stored in this array.
   *  Sufficient memory (for as many pixels as there are in the template grid)
   *  must be allocated there.
   */
  virtual void InterpolateImage( const size_t idx, byte* const destination );

private:
  /// Thread function parameters for image interpolation.
  class InterpolateImageThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParameters<Self>
  {
  public:
    /// Index of the image to be interpolated.
    size_t m_Idx;

    /// Pointer to storage that will hold the reformatted pixel data.
    byte* m_Destination;

    const Vector3D* m_HashX;
    const Vector3D* m_HashY;
    const Vector3D* m_HashZ;
  };

  /// Image interpolation thread function.
  static void InterpolateImageThread( void *const threadParameters, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t );
  
  /// Image interpolation with probabilistic sampling thread function.
  static void InterpolateImageProbabilisticThread( void *const threadParameters, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t );

  friend ClassStream& operator<<( ClassStream& stream, const AffineGroupwiseRegistrationRMIFunctional& func );
  friend ClassStream& operator>>( ClassStream& stream, AffineGroupwiseRegistrationRMIFunctional& func );
};

/// Class stream write function.
ClassStream& operator<<( ClassStream& stream, const AffineGroupwiseRegistrationRMIFunctional& func );

/// Class stream read function.
ClassStream& operator>>( ClassStream& stream, AffineGroupwiseRegistrationRMIFunctional& func );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkAffineGroupwiseRegistrationRMIFunctional_h_included_
