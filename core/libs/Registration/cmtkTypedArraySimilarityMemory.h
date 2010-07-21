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

#ifndef __cmtkTypedArraySimilarityMemory_h_included_
#define __cmtkTypedArraySimilarityMemory_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkJointHistogram.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Memory for typed array similarity computation.
 * This class provides for optional persistent memory between similarity 
 * computations. This ensures, for example, identical value ranges in 
 * histogram-based similarity measures.
 */
class TypedArraySimilarityMemory
{
public:
  /** Initialize a similarity instance with memory.
   * By instantiating a class object, the otherwise static member functions
   * can be given a memory that coordinates their behaviour between calls. For
   * example, we can make sure that all evaluations of Mutual Information use
   * the same histogram resolution.
   *\param repeatCheck If this flag is set, the object will repeat the range
   * check for every call to GetRangeX or GetRangeY. If the current data
   * range exceeds the one stored in this object, the latter will be adapted
   * accordingly.
   *\param safetyMargin The percentage of the actual value range that will be 
   * stored in this object.
   */
  TypedArraySimilarityMemory( const bool repeatCheck = true )
    : ValidX( false ), ValidY( false ), MinNumBins( 8 ), MaxNumBins( 128 )
  {  
    RepeatCheck = repeatCheck;
  }

  /** Get range of X distribution.
   * If this object is not yet initialized, the given array is queried for
   * its value range, and this object is initialized accordingly.
   */
  const Types::DataItemRange GetRangeX( const TypedArray* array, const size_t numBins );

  /** Get range of Y distribution.
   * If this object is not yet initialized, the given array is queried for
   * its value range, and this object is initialized accordingly.
   */
  const Types::DataItemRange GetRangeY( const TypedArray* array, const size_t numBins );
  
  /// Set minimum number of histogram bins.
  void SetMinNumBins( const size_t minNumBins ) { MinNumBins = minNumBins; }

  /// Set maximum number of histogram bins.
  void SetMaxNumBins( const size_t maxNumBins ) { MaxNumBins = maxNumBins; }

  /// Create histogram based on memorized settings.
  JointHistogram<unsigned int>::SmartPtr CreateHistogram( const TypedArray* array0, const TypedArray* array1 );

private:
  /// Repeat range check with each call to GetRangeX and GetRangeY.
  bool RepeatCheck;

  /// Flag whether memory for X distribution is already initialized.
  bool ValidX;

  /// Remembered range of X values.
  Types::DataItemRange RangeX;

  /// Remembered number of bins for the X distribution.
  size_t NumberBinsX;

  /// Flag whether memory for X distribution is already initialized.
  bool ValidY;

  /// Remembered range of Y values.
  Types::DataItemRange RangeY;

  /// Remembered number of bins for the Y distribution.
  size_t NumberBinsY;

  /// Minimum number of histogram bins.
  size_t MinNumBins;

  /// Maximum number of histogram bins.
  size_t MaxNumBins;

  /// Allow similarity computation class access.
  friend class TypedArraySimilarity;
};
  
//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypedArraySimilarityMemory_h_included_
