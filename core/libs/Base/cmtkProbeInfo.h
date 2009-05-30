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

#ifndef __cmtkProbeInfo_h_included_
#define __cmtkProbeInfo_h_included_

#include <cmtkconfig.h>

#include <cmtkTypes.h>
#include <cmtkVector3D.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Information on volume node data.
 * This class represents all information needed for trilinear interpolation and
 * partial derivative computation. It holds a grid location and the node data
 * on the corners of the grid cube containing that location.
 */
class ProbeInfo 
{
public:
  /** Node data on cube corners around given location.
   * The order of the data elements is as follows: (x0,y0,z0), (x1,y0,z0),
   * (x0,y1,z0), (x1,y1,z0), (x0,y0,z1), (x1,y0,z1), (x0,y1,z1), (x1,y1,z1).
   */
  Types::DataItem Values[8];

  /** Dimensions of the grid cube.
   */
  Types::Coordinate Deltas[3];

  /** Relative location of the probed coordinate.
   * Represented are the relative distances of the probed coordinate from the
   * surrounding cube's faces. Values range from 0 to 1. The order of the
   * coefficients is dX0, dY0, dZ0, dX1, dY1, dZ1. Lower x-coordinates for
   * example are weighted with dX0, others respectively.
   */
  Types::Coordinate Offsets[6];

  /** The real location that was probed.
   */
  Vector3D Location;

  /** Return data value at desired location.
   * Trlinear interpolation is performed using the pre-calculated coefficients.
   */
  Types::DataItem GetValueTrilinear () const 
  {
    return static_cast<Types::DataItem>( Offsets[2]*(Offsets[1]*(Offsets[0]*Values[0]+Offsets[3]*Values[1])+ 
					     Offsets[4]*(Offsets[0]*Values[2]+Offsets[3]*Values[3]))+
				 Offsets[5]*(Offsets[1]*(Offsets[0]*Values[4]+Offsets[3]*Values[5])+ 
					     Offsets[4]*(Offsets[0]*Values[6]+Offsets[3]*Values[7]) ) );
  }

  /** Return data value at desired location using only values above threshold.
   */
  Types::DataItem GetValueTrilinearAbove ( const Types::DataItem threshold ) const 
  {
    Types::DataItem result = 0;
    for ( int idx=0; idx<8; ++idx ) 
      {
      if ( Values[idx] >= threshold )
	result += static_cast<Types::DataItem>( GetWeight( idx ) * Values[idx] );
      }
    return result;
  }

  /** Return data value at desired location using only values below threshold.
   */
  Types::DataItem GetValueTrilinearBelow ( const Types::DataItem threshold ) const 
  {
    Types::DataItem result = 0;
    for ( int idx=0; idx<8; ++idx ) {
      if ( Values[idx] < threshold )
	result += static_cast<Types::DataItem>( GetWeight( idx ) * Values[idx] );
    }
    return result;
  }

  /** Return number of values above (or equal to) threshold.
   */
  byte GetNumberOfValuesAbove ( const Types::DataItem threshold ) const 
  {
    byte result = 0;
    for ( int idx=0; idx<8; ++idx )
      if ( Values[idx] >= threshold )
	++result;
    return result;
  }

  /** Return number of values below threshold.
   */
  byte GetNumberOfValuesBelow ( const Types::DataItem threshold ) const 
  {
    byte result = 0;
    for ( int idx=0; idx<8; ++idx )
      if ( Values[idx] < threshold )
	++result;
    return result;
  }

  /// Return relative weight of given index.
  Types::Coordinate GetWeight( const int index ) const 
  {
    switch ( index ) 
      {
      case 0: return Offsets[2] * Offsets[1] * Offsets[0];
      case 1: return Offsets[2] * Offsets[1] * Offsets[3];
      case 2: return Offsets[2] * Offsets[4] * Offsets[0];
      case 3: return Offsets[2] * Offsets[4] * Offsets[3];
      case 4: return Offsets[5] * Offsets[1] * Offsets[0];
      case 5: return Offsets[5] * Offsets[1] * Offsets[3];
      case 6: return Offsets[5] * Offsets[4] * Offsets[0];
      case 7: return Offsets[5] * Offsets[4] * Offsets[3];
      }
    return 0;
  }

  /** Return partial derivatives of node data w.r.t. grid dimensions.
   * The data on the cube faces is used for a finite-difference approximation
   * of the first-order derivatives w.r.t. the grid dimensions x, y, and z.
   *@param d The object of type Vector3D in which the result is to be stored.
   *@return A reference to the destination parameter.
   */
  Vector3D& GetPartialDerivatives ( Vector3D& d ) const 
  {
    d.XYZ[0] = static_cast<Types::Coordinate>( Deltas[0] * 
					   ( Offsets[2]*(Offsets[1]*Values[1]+Offsets[4]*Values[3])+
					     Offsets[5]*(Offsets[1]*Values[5]+Offsets[4]*Values[7])-
					     Offsets[2]*(Offsets[1]*Values[0]+Offsets[4]*Values[2])-
					     Offsets[5]*(Offsets[1]*Values[4]+Offsets[4]*Values[6])) );
    
    d.XYZ[1] = static_cast<Types::Coordinate>( Deltas[1] * 
					   ( Offsets[2]*(Offsets[0]*Values[2]+Offsets[3]*Values[3])+
					     Offsets[5]*(Offsets[0]*Values[6]+Offsets[3]*Values[7])-
					     Offsets[2]*(Offsets[0]*Values[0]+Offsets[3]*Values[1])-
					     Offsets[5]*(Offsets[0]*Values[4]+Offsets[3]*Values[5])) );
    
    d.XYZ[2] = static_cast<Types::Coordinate>( Deltas[2] *
					   ( Offsets[1]*(Offsets[0]*Values[4]+Offsets[3]*Values[5])+
					     Offsets[4]*(Offsets[0]*Values[6]+Offsets[3]*Values[7])-
					     Offsets[1]*(Offsets[0]*Values[0]+Offsets[3]*Values[1])-
					     Offsets[4]*(Offsets[0]*Values[2]+Offsets[3]*Values[3])) );
    return d;
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkProbeInfo_h_included_
