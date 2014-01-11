/*
//
//  Copyright 2012, 2014 SRI International
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

#ifndef __cmtkFitPolynomialToLandmarks_h_included_
#define __cmtkFitPolynomialToLandmarks_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkPolynomialXform.h>
#include <Base/cmtkLandmarkPairList.h>
#include <Base/cmtkMatrix3x3.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Fit polynomial transformation to a set of landmark pairs.
 * The fitting algorithm uses SVD and works incrementally. That is, lower polynomial degrees are fitted first, followed by higher ones.
 */
class FitPolynomialToLandmarks
{
public:
  /// This class.
  typedef FitPolynomialToLandmarks Self;

  /// Constructor.
  FitPolynomialToLandmarks( const LandmarkPairList& landmarkPairs /*!< Landmark pairs to fit transformation to */, const byte degree /*!< Degree of the fitted polynomial */ );

  /// Return the affine transformation.
  PolynomialXform::SmartPtr GetPolynomialXform()
  {
    return this->m_PolynomialXform;
  }
  
  /// Return the constant affine transformation.
  PolynomialXform::SmartConstPtr GetPolynomialXform() const
  {
    return this->m_PolynomialXform;
  }
  
private:
  /// The fitted transformation.
  PolynomialXform::SmartPtr m_PolynomialXform;
};

} // namespace

#endif // #ifndef __cmtkFitPolynomialToLandmarks_h_included_
