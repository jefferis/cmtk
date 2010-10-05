/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkHypothesisTests_h_included_
#define __cmtkHypothesisTests_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMacros.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkVolume.h>
#include <Base/cmtkTypedArray.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/// Statistical hypothesis testing between groups of images.
class HypothesisTests
{
public:
  /// Test Jacobian maps of two populations for statistical independence.
  static TypedArray::SmartPtr GetUnpairedTwoTailedTTest
  ( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, TypedArray::SmartPtr* tstatData, TypedArray::SmartPtr* avgXData, TypedArray::SmartPtr* avgYData,
    const TypedArray* mask = NULL );

  /// Test parameter maps of two populations for statistical independence.
  static TypedArray::SmartPtr GetPairedTwoTailedTTest
  ( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, TypedArray::SmartPtr* tstatData, TypedArray::SmartPtr* avgXData, TypedArray::SmartPtr* avgYData, 
    const TypedArray* mask = NULL );

  /// Get pixel-wise correlation between two sets of input images.
  static TypedArray::SmartPtr 
  GetPairedCorrelation( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, TypedArray::SmartPtr* pData = NULL, const TypedArray* mask = NULL );

  /// Test mean of Jacobian map of a single population for difference from zero.
  static TypedArray::SmartPtr 
  GetOneSampleTTest( std::vector<TypedArray::SmartPtr>& dataX, TypedArray::SmartPtr* tstatData, TypedArray::SmartPtr* avgXData, const TypedArray* mask = NULL );

  /// Get pixelwise heritability of two populations.
  static TypedArray::SmartPtr 
  GetHeritability( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, const TypedArray* mask = NULL );
  
  /** Get pixelwise z-scores.
    * The X distribution is taken as the "true" or "reference" distribution.
    * The Y distribution is taken as the "test" or "sample" distribution.
    */
  static TypedArray::SmartPtr GetZScores( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, const TypedArray* mask = NULL );
  
  /// Get pixelwise genetic covariance from MZ and DZ twin data.
  static TypedArray::SmartPtr GetGeneticCovariance( std::vector<TypedArray::SmartPtr>& dataMZ, std::vector<TypedArray::SmartPtr>& dataDZ, const TypedArray* mask = NULL );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkHypothesisTests_h_included_
