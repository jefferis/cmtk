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

#ifndef __cmtkImagePairAffineRegistration_h_included_
#define __cmtkImagePairAffineRegistration_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkInterpolator.h>

#include <Registration/cmtkImagePairRegistration.h>
#include <Registration/cmtkMakeInitialAffineTransformation.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for affine multi-resolution voxel registration.
 */
class ImagePairAffineRegistration : 
  /// Inherit general voxel registration interface and functionality.
  public ImagePairRegistration 
{
protected:
  /** Flag for initial alignment of volume centers.
   */
  cmtkGetSetMacro(MakeInitialAffineTransformation::Mode,Initializer);
  
  /// Flag: use symmetric registration functional for simultaneous forward/inverse transformation computation.
  bool m_SymmetricFwdBwd;

  /// Flag whether to adjust floating image histogram to match reference image.
  cmtkGetSetMacro(bool,MatchFltToRefHistogram);

  /** Numbers of degrees of freedom.
   * This list contains the numbers of degrees of freedom for every resolution
   * level. Registration is repeated with the same data as many times as there
   * are entries in this list. If the derived classes do not set any entries,
   * InitRegistration() will push a "6"  into the list, resulting in an affine
   * registration.
   */
  std::vector<short> NumberDOFs;

  /** Numbers of degrees of freedom for the final resolution level.
   * Just as "NumberDOFs", this list defines the sequence of numbers of degrees
   * of freedom for the finest resolution level.
   */
  std::vector<short> NumberDOFsFinal;

  /** Restrict to in-plane transformations.
   * If this is set to a number between 0 and 2, it indexes a coordinate axis,
   * and computed transformations are restricted to be in-plane with respect to
   * planes perpendicular to that axis.
   */
  int m_RestrictToInPlane;

  /** Initialize registration.
   * This function is called by Register before any other operations. It can
   * be overloaded to open status dialog windows, etc. Derived implementations
   * should call their base class' InitRegistration first.
   *\return Overriding functions should return a value other than 
   * CALLBACK_OK if the registration is to be interrupted.
   */
  virtual CallbackResult InitRegistration ();

  /** Enter a resolution level.
   * This function mainly determines the next effective number of degrees of
   * freedom of the optimization. It sets the transformation object accordingly
   * and writes a comment to the "Callback" object. Afterwards, the inherited 
   * EnterResolution() function is called.
   */
  virtual void EnterResolution( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, const int level, const int total );

  /** Finish resolution level.
   * This function determines whether there are any remaining numbers in the
   * effective NumberDOFs list. If so, the list iterator is advanced and 0 is
   * returned, indicating that the current resolution level is to be repeated.
   * With no number left, 1 is returned indicating that the current level has
   * been completed.
   */
  virtual int DoneResolution( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, const int level, const int total );

public:
  /// This class.
  typedef ImagePairAffineRegistration Self;

  /// Parent class.
  typedef ImagePairRegistration Superclass;

  /** Default constructor. 
   */
  ImagePairAffineRegistration ();

  /** Destructor.
   */
  virtual ~ImagePairAffineRegistration ();

  /// Return final transformation.
  AffineXform::SmartPtr GetTransformation() const;
  
  /// Get reformatted floating image.
  const UniformVolume::SmartPtr GetReformattedFloatingImage( Interpolators::InterpolationEnum interpolator = Interpolators::LINEAR ) const;

  /// Add a number to the general list of numbers of DOFs.
  void AddNumberDOFs( const int numberDOFs ) 
  {
    NumberDOFs.push_back( numberDOFs );
  }
  
  /// Add a number to the list of numbers of DOFs for the last level.
  void AddNumberDOFsFinal( const int numberDOFs ) 
  {
    NumberDOFsFinal.push_back( numberDOFs );
  }

private:
  /// Iterator for NumberDOFs and NumberDOFsFinal
  std::vector<short>::iterator NumberDOFsIterator;

  /// Base class for registration level parameters.
  class LevelParameters
    /// Inherit from superclass parameters.
    : public Superclass::LevelParameters
  {
  public:
    /// Constructor: take image resolution.
    LevelParameters( const Types::Coordinate resolution ) : m_Resolution( resolution ) {}

    /// Image resolution for this level.
    Types::Coordinate m_Resolution;
  };

  /** Create functional with all settings for next level.
   */
  virtual Functional* MakeFunctional( const int level, const Superclass::LevelParameters* );
};

//@}

} // namespace cmtk

#endif // __cmtkImagePairAffineRegistration_h_included_
