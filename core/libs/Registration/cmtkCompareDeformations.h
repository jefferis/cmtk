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

#ifndef __cmtkCompareDeformations_h_included_
#define __cmtkCompareDeformations_h_included_

#include <cmtkconfig.h>

#include <cmtkMacros.h>
#include <cmtkWarpXform.h>
#include <cmtkVolume.h>
#include <cmtkTypedArray.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

class CompareDeformations
{
public:
  void AddDeformation0( WarpXform::SmartPtr& warpXform, WarpXform::SmartPtr& warpToRef )
  {
    Deformations0.push_back( warpXform );
    WarpToRef0.push_back( warpToRef );
  }
  
  void AddDeformation1( WarpXform::SmartPtr& warpXform, WarpXform::SmartPtr& warpToRef  ) 
  {
    Deformations1.push_back( warpXform );
    WarpToRef1.push_back( warpToRef );
  }
  
  void SetDeformations0( std::vector<WarpXform::SmartPtr>& deformations, std::vector<WarpXform::SmartPtr>& warpToRef ) 
  {
    Deformations0 = deformations;
    WarpToRef0 = warpToRef;
  }
  
  void SetDeformations1( std::vector<WarpXform::SmartPtr>& deformations, std::vector<WarpXform::SmartPtr>& warpToRef ) 
  {
    Deformations1 = deformations;
    WarpToRef1 = warpToRef;
  }
  
  /// Test Jacobian maps of two populations for statistical independence.
  static TypedArray* GetJacobianTTest( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, TypedArray** tstatData, TypedArray** avgXData, TypedArray** avgYData,
				       const TypedArray* mask = NULL );

  /// Test parameter maps of two populations for statistical independence.
  static TypedArray* GetPairedTTest( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, TypedArray** tstatData, TypedArray** avgXData, TypedArray** avgYData, 
				     const TypedArray* mask = NULL );

  /// Get pixel-wise correlation between two sets of input images.
  static TypedArray* GetPairedCorrelation( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, TypedArray** pData = NULL, const TypedArray* mask = NULL );

  /// Test mean of Jacobian map of a single population for difference from zero.
  static TypedArray* GetOneSampleTTest( std::vector<TypedArray::SmartPtr>& dataX, TypedArray** tstatData, TypedArray** avgXData, const TypedArray* mask = NULL );

  /// Get pixelwise heritability of two populations.
  static TypedArray* GetHeritability( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, const TypedArray* mask = NULL );

  /** Get pixelwise z-scores.
    * The X distribution is taken as the "true" or "reference" distribution.
    * The Y distribution is taken as the "test" or "sample" distribution.
    */
  static TypedArray* GetZScores( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, const TypedArray* mask = NULL );

  /// Get pixelwise genetic covariance from MZ and DZ twin data.
  static TypedArray* GetGeneticCovariance( std::vector<TypedArray::SmartPtr>& dataMZ, std::vector<TypedArray::SmartPtr>& dataDZ, const TypedArray* mask = NULL );

private:
  std::vector<WarpXform::SmartPtr> Deformations0;
  std::vector<WarpXform::SmartPtr> WarpToRef0;

  std::vector<WarpXform::SmartPtr> Deformations1;
  std::vector<WarpXform::SmartPtr> WarpToRef1;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkCompareDeformations_h_included_
