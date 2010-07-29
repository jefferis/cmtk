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

#ifndef __cmtkFilterMask_h_included_
#define __cmtkFilterMask_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkTypes.h"
#include "Base/cmtkUnits.h"
#include "Base/cmtkFixedVector.h"
#include "Base/cmtkMathUtil.h"

#include "System/cmtkSmartPtr.h"

#include <cmath>
#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Filter mask pixel entry.
 * This class handles a single entry in a pre-computed filter mask for
 * multidimensional images with relative coordinates and filter coeffiecient 
 * for a single pixel.
 */
template<int DIM>
class FilterMaskPixel
{
public:
  /// This class.
  typedef FilterMaskPixel Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Default constructor.
  FilterMaskPixel() {}

  /// Explicit constructor.
  FilterMaskPixel
  ( const FixedVector<DIM,int>& location, const int relativeIndex, const Types::DataItem coefficient ) 
    : Location( location ),
      RelativeIndex ( relativeIndex ),
      Coefficient( coefficient )
  {}

  /// Relative location of this pixel.
  FixedVector<DIM,int> Location;

  /// Relative index of this pixel in source image from center of kernel.
  int RelativeIndex;

  /// Filter coefficient.
  Types::DataItem Coefficient;

  /// Cached source index of image pixel for this element.
  int PixelIndex;

  /** Active flag.
   * This flag can be used in client code to flag pixels in the filter mask as
   * either valid (i.e., inside domain) or invalid (outsside domain). This
   * eliminates repeated range checking in cases where there are multiple
   * iterations over the filter mask at one location.
   */
  bool Valid;
};

/** Filter mask.
 * This class handles pre-computed filter masks for multidimensional images
 * with relative pixel coordinates and filter coeffiecients.
 */
template<int DIM>
class FilterMask : 
  /// Inherit from STL container.
  public std::vector< FilterMaskPixel<DIM> >
{
public:
  /// This class type.
  typedef FilterMask<DIM> Self;

  /// Direct parent class.
  typedef std::vector< FilterMaskPixel<DIM> > Superclass;
  
  /// Smert pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  template<class F>
  FilterMask( const FixedVector<DIM,int>& dims, const FixedVector<DIM,Types::Coordinate>& deltas, const Types::Coordinate radius, F filter ) 
  {
    FixedVector<DIM,int> pixel;
    FixedVector<DIM,int> width;
    FixedVector<DIM,Types::Coordinate> position;

    for ( int dim = 0; dim < DIM; ++dim ) 
      {
      pixel[dim] = - (width[dim] = 1+static_cast<int>( radius / deltas[dim] ));
      position[dim] = pixel[dim] * deltas[dim];
      }

    bool done = false;
    while ( ! done ) 
      {
      // increment the DIM-digit pixel index counter including overflow.
      for ( int dim = 0; dim < DIM; ++dim ) 
	{
	++pixel[dim];
	if ( pixel[dim] <= width[dim] ) 
	  {
	  // no overflow, leave for loop since we're done
	  dim = DIM;
	  } 
	else
	  {
	  if ( dim+1 == DIM ) 
	    // was this the last dimension? if so, leave while() loop
	    done = true;
	  else 
	    { 
	    // no, then reset this dimension and repeat loop to increment next
	    pixel[dim] = -width[dim];
	    }
	  }
	}
      // are we done with the kernel?
      if ( ! done ) 
	{
	// no, then compute Euclidean distance from center
	Types::Coordinate distance = 0.0;
	for ( int dim = 0; dim < DIM; ++dim ) 
	  {
	  position[dim] = pixel[dim] * deltas[dim];
	  distance += position[dim] * position[dim];
	  }
	distance = sqrt( distance );
	// if distance is within radius then add a pixel to the filter mask
	if ( distance < radius ) 
	  {
	  const int index = pixel[0] + dims[0] * (pixel[1] + dims[1] * pixel[2] );
	  this->push_back( FilterMaskPixel<DIM>( pixel, index, filter( position ) ) );
	  }
	}
      }
  }
  
  /// Gaussian filter as an example of a concrete filter implementation.
  class Gaussian 
  {
  public:
    /// Constructor.
    Gaussian( const Units::GaussianSigma& standardDeviation ) 
    {
      InvStandardDeviation = 1.0 / standardDeviation.Value();
      NormFactor = 1.0 / (sqrt(2.0 * M_PI) * standardDeviation.Value());
    }
    
    /// Get filter coefficient at relative location from filter center.
    Types::DataItem operator() ( const FixedVector<DIM,Types::Coordinate>& relativePosition ) 
    {
      Types::Coordinate distance = 0;
      for ( int i = 0; i < DIM; ++i ) 
	distance += relativePosition[i] * relativePosition[i];
      return static_cast<Types::DataItem>( NormFactor * exp( -distance * MathUtil::Square( InvStandardDeviation ) / 2 ) );
    }
    
  private:
    /// Standard deviation.
    Types::Coordinate InvStandardDeviation;

    /// Gaussian normalization factor.
    Types::Coordinate NormFactor;
  };
};

} // namespace  cmtk

#endif // #ifndef __cmtkFilterMask_h_included_
