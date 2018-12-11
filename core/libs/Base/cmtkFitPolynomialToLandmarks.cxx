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

#include "cmtkFitPolynomialToLandmarks.h"

cmtk::FitPolynomialToLandmarks::FitPolynomialToLandmarks( const LandmarkPairList& landmarkPairs, const byte degree )
{
  // first, get the centroids in "from" and "to" space
  PolynomialXform::SpaceVectorType cFrom( 0.0 );
  PolynomialXform::SpaceVectorType cTo( 0.0 );
  
  size_t nLandmarks = 0;
  for ( LandmarkPairList::const_iterator it = landmarkPairs.begin(); it != landmarkPairs.end(); ++it )
    {
    cFrom += it->m_Location;
    cTo += it->m_TargetLocation;
    ++nLandmarks;
    }
  
  cFrom /= nLandmarks;
  cTo /= nLandmarks;
  
  // put everything together
  this->m_PolynomialXform = PolynomialXform::SmartPtr( new PolynomialXform( degree ) );
  this->m_PolynomialXform->SetCenter( cFrom );

  // Fit incrementally - start with lower degrees.
  // For degree = 0 we actually get relative translation of the centers of mass of the landmarks in the two spaces
  for ( byte fitDegree = 0; fitDegree <= degree; ++fitDegree )
    {
    // what is the first parameter at the current degree?
    const size_t firstParameter = PolynomialHelper::GetNumberOfMonomials( fitDegree-1 );

    // how many parameters in addition to prior degree?
    const size_t nParameters = PolynomialHelper::GetNumberOfMonomials( fitDegree ) - firstParameter;
    
    // set up matrix for SVD-based linear regression
    Matrix2D<double> U( nLandmarks, nParameters ); 
    
    std::vector<PolynomialXform::SpaceVectorType> residuals( nLandmarks );
    size_t lmIdx = 0;
    for ( LandmarkPairList::const_iterator it = landmarkPairs.begin(); it != landmarkPairs.end(); ++it, ++lmIdx )
      {
      // current residual for this landmark
      residuals[lmIdx] = it->m_TargetLocation - this->m_PolynomialXform->Apply( it->m_Location );
      
      // monomials for this landmark and each of the currently fitted parameters
      for ( size_t paramIdx = 0; paramIdx < nParameters; ++paramIdx )
	{
	U[lmIdx][paramIdx] = this->m_PolynomialXform->GetMonomialAt( firstParameter+paramIdx, it->m_Location );
	}
      }

    // do SVD of the monomial matrix
    Matrix2D<double> V( nParameters, nParameters );
    std::vector<double> W( nParameters );
    MathUtil::SVD( U, W, V );
    
    // for each dimension separately, use the SVD to get parameters that minimize the residuals (keeping in mind that poly xform is additive, i.e., adding monomials makes an additive contribution to existing transformation)
    std::vector<double> params;
    for ( int dim = 0; dim < 3; ++dim )
      {
      // first, extract residuals for current spatial dimension, dim
      std::vector<double> dimResidual( nLandmarks );
      for ( size_t lmIdx = 0; lmIdx < nLandmarks; ++lmIdx )
	dimResidual[lmIdx] = residuals[lmIdx][dim];

      // now apply the previously computed SVD to get coefficients for the additional monomials at this level
      MathUtil::SVDLinearRegression( U, W, V, dimResidual, params );
      for ( size_t paramIdx = 0; paramIdx < nParameters; ++paramIdx )
	{
	this->m_PolynomialXform->m_Parameters[3*(firstParameter+paramIdx) + dim] = params[paramIdx];
	}
      }
    }
}
