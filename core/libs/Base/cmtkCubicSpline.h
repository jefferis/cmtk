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

#ifndef __cmtkCubicSpline_h_included_
#define __cmtkCubicSpline_h_included_

#include <cmtkconfig.h>

#include <cmtkTypes.h>
#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class computing cubic splines.
 * This class is used for example by tri-cubic intensity interpolation and
 * B-spline deformations. It supports both approximating splines (for
 * deformations) and interpolating splines (for intensity interpolation).
 */
class CubicSpline 
{
public:
  /// Compute a value of the 0th approximating spline function.
  static Types::Coordinate ApproxSpline0 ( const Types::Coordinate t ) 
  {
    return ( MathUtil::Square(1-t) * (1-t) ) / 6;
  }

  /// Compute a value of the 1st approximating spline function.
  static Types::Coordinate ApproxSpline1 ( const Types::Coordinate t ) 
  {
    return ( 4 + MathUtil::Square(t) * ( 3 * t - 6 ) ) / 6;
  }
  
  /// Compute a value of the 2nd approximating spline function.
  static Types::Coordinate ApproxSpline2 ( const Types::Coordinate t ) 
  {
    return ( 1 + t * (3 + t * (3 - 3*t))) / 6;
  }

  /// Compute a value of the 3rd approximating spline function.
  static Types::Coordinate ApproxSpline3 ( const Types::Coordinate t ) 
  {
    return t*t*t/6;
  }  

  /// Compute the value of a given approximating spline function.
  static Types::Coordinate ApproxSpline ( const int k, const Types::Coordinate t ) 
  {
    switch (k) 
      {
      case 0: return ApproxSpline0( t );
      case 1: return ApproxSpline1( t );
      case 2: return ApproxSpline2( t );
      case 3: return ApproxSpline3( t );
      default: return 0;
      }
#ifdef MSDOS
    return 0;
#endif
  }

  /// Compute the derivative of the 0th approximating spline function.
  static Types::Coordinate DerivApproxSpline0 ( const Types::Coordinate t ) 
  {
    return  -MathUtil::Square(1-t) / 2;
  }
  
  /// Compute derivative of the 1st approximating spline function.
  static Types::Coordinate DerivApproxSpline1 ( const Types::Coordinate t ) 
  {
    return 3*t*t/2-2*t;
  }
  
  /// Compute derivative of the 2nd approximating spline function.
  static Types::Coordinate DerivApproxSpline2 ( const Types::Coordinate t ) 
  {
    return ( 1 + 2*t - 3*t*t ) / 2;
  }

  /// Compute derivative of the 3rd approximating spline function.
  static Types::Coordinate DerivApproxSpline3 ( const Types::Coordinate t ) 
  {
    return t*t/2;
  }  
  
  /// Compute the derivative of a given approximating spline function.
  static Types::Coordinate DerivApproxSpline ( const int k, const Types::Coordinate t ) 
  {
    switch (k) 
      {
      case 0: return DerivApproxSpline0( t );
      case 1: return DerivApproxSpline1( t );
      case 2: return DerivApproxSpline2( t );
      case 3: return DerivApproxSpline3( t );
      default: return 0;
      }
#ifdef MSDOS
    return 0;
#endif
  }
  
  /// Compute the second derivative of the 0th approximating spline function.
  static Types::Coordinate SecondDerivApproxSpline0 ( const Types::Coordinate t ) 
  {
    return  1 - t;
  }
  
  /// Compute second derivative of the 1st approximating spline function.
  static Types::Coordinate SecondDerivApproxSpline1 ( const Types::Coordinate t ) 
  {
    return 3 * t - 2;
  }

  /// Compute second derivative of the 2nd approximating spline function.
  static Types::Coordinate SecondDerivApproxSpline2 ( const Types::Coordinate t ) 
  {
    return 1 - 3 * t;
  }

  /// Compute second derivative of the 3rd approximating spline function.
  static Types::Coordinate SecondDerivApproxSpline3 ( const Types::Coordinate t ) 
  {
    return t;
  }  

  /// Compute the second derivative of a given approximating spline function.
  static Types::Coordinate SecondDerivApproxSpline ( const int k, const Types::Coordinate t ) 
  {
    switch (k) 
      {
      case 0: return SecondDerivApproxSpline0( t );
      case 1: return SecondDerivApproxSpline1( t );
      case 2: return SecondDerivApproxSpline2( t );
      case 3: return SecondDerivApproxSpline3( t );
      default: return 0;
      }
#ifdef MSDOS
    return 0;
#endif
  }
  
  /// Compute a value of the 0th interpolating spline function.
  static Types::Coordinate InterpSpline0 ( const Types::Coordinate t ) 
  {
    return (Types::Coordinate)( t * ( 0.5 * ( -1 + t * ( 2 - t ) ) ) );
  }

  /// Compute a value of the 1st interpolating spline function.
  static Types::Coordinate InterpSpline1 ( const Types::Coordinate t ) 
  {
    return (Types::Coordinate)( 0.5 * ( 3 * t - 5 ) * t*t + 1 );
  }
  
  /// Compute a value of the 2nd interpolating spline function.
  static Types::Coordinate InterpSpline2 ( const Types::Coordinate t ) 
  {
    return (Types::Coordinate)( 0.5 * t * ( 1 + t * ( 4  - 3 * t ) ) );
  }

  /// Compute a value of the 3rd interpolating spline function.
  static Types::Coordinate InterpSpline3 ( const Types::Coordinate t ) 
  {
    return (Types::Coordinate)( 0.5 * ( t*t * ( t-1 ) ) );
  }  

  /// Compute the value of a given interpolating spline function.
  static Types::Coordinate InterpSpline ( const int k, const Types::Coordinate t ) 
  {
    switch (k) 
      {
      case 0: return InterpSpline0( t );
      case 1: return InterpSpline1( t );
      case 2: return InterpSpline2( t );
      case 3: return InterpSpline3( t );
      default: return 0;
      }
#ifdef MSDOS
    return 0;
#endif
  }
  
  /// Compute derivative of the 0th interpolating spline function.
  static Types::Coordinate DerivInterpSpline0 ( const Types::Coordinate t ) 
  {
    return ((-3 * t + 4) * t - 1) / 2;
  }
  
  /// Compute derivative of the 1st interpolating spline function.
  static Types::Coordinate DerivInterpSpline1 ( const Types::Coordinate t ) 
  {
    return ( ( 9 * t - 10 ) * t ) / 2;
  }
  
  /// Compute derivative of the 2nd interpolating spline function.
  static Types::Coordinate DerivInterpSpline2 ( const Types::Coordinate t ) 
  {
    return ( ( -9 * t + 8 ) * t + 1 ) / 2;
  }

  /// Compute derivative of the 3rd interpolating spline function.
  static Types::Coordinate DerivInterpSpline3 ( const Types::Coordinate t ) 
  {
    return (( 3 * t - 2 ) * t) / 2;
  }  
  
  /// Compute derivative of a given interpolating spline function.
  static Types::Coordinate DerivInterpSpline ( const int k, const Types::Coordinate t ) 
  {
    switch (k) 
      {
      case 0: return DerivInterpSpline0( t );
      case 1: return DerivInterpSpline1( t );
      case 2: return DerivInterpSpline2( t );
      case 3: return DerivInterpSpline3( t );
      default: return 0;
      }
#ifdef MSDOS
    return 0;
#endif
  }
  
  /// Compute second derivative of the 0th interpolating spline function.
  static Types::Coordinate SecondDerivInterpSpline0 ( const Types::Coordinate t ) 
  {
    return -3 * t + 2;
  }
  
  /// Compute second derivative of the 1st interpolating spline function.
  static Types::Coordinate SecondDerivInterpSpline1 ( const Types::Coordinate t ) 
  {
    return 9 * t - 5;
  }

  /// Compute second derivative of the 2nd interpolating spline function.
  static Types::Coordinate SecondDerivInterpSpline2 ( const Types::Coordinate t )
  {
    return -9 * t + 4;
  }

  /// Compute second derivative of the 3rd interpolating spline function.
  static Types::Coordinate SecondDerivInterpSpline3 ( const Types::Coordinate t ) 
  {
    return 3 * t - 1;
  }  
  
  /// Compute second derivative of a given interpolating spline function.
  static Types::Coordinate SecondDerivInterpSpline ( const int k, const Types::Coordinate t ) 
  {
    switch (k) 
      {
      case 0: return SecondDerivInterpSpline0( t );
      case 1: return SecondDerivInterpSpline1( t );
      case 2: return SecondDerivInterpSpline2( t );
      case 3: return SecondDerivInterpSpline3( t );
      default: return 0;
      }
#ifdef MSDOS
    return 0;
#endif
  }
};

} // namespace

#endif // #ifndef __cmtkCubicSpline_h_included_
