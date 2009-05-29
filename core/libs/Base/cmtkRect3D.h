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

#ifndef __cmtkRect3D_h_included_
#define __cmtkRect3D_h_included_

#include <cmtkconfig.h>
#include <cmtkTypes.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Index type for grid elements.
typedef unsigned short GridIndexType;

/// 3-D rectangular box range.
template<class T,class T2>
class TemplateRect3D_Base
{
public:
  /** x coordinate range */
  T startX, endX;
  
  /** y coordinate range */
  T startY, endY;
  
  /** z coordinate range */
  T startZ, endZ;
  
  /** Default constructor.
   * Do not do anything.
   */
  TemplateRect3D_Base() {};
  
  /** Initializing constructor.
   * Set all fields to zero.
   */
  TemplateRect3D_Base( int )
  { startX = endX = startY = endY = startZ = endZ = 0; }
  
  /** Compare operator.
   *@param other The Rect3D object to compare this to.
   *@return Returns 1 if and only if other object defines the same box.
   */
  int operator== ( const TemplateRect3D_Base<T,T2>& other ) const 
  {
    return ( (startX==other.startX) && (endX==other.endX) &&
	     (startY==other.startY) && (endY==other.endY) &&
	     (startZ==other.startZ) && (endZ==other.endZ) );
  }
  
  /** Assignment operator.
   * Copy all fields to target object.
   */
  TemplateRect3D_Base<T,T2>& operator= ( const TemplateRect3D_Base<T,T2>& other ) 
  {
    startX = other.startX; endX = other.endX;
    startY = other.startY; endY = other.endY;
    startZ = other.startZ; endZ = other.endZ;
    return *this;
  }
  
  /** Size operator.
   *@return The size of this box. Edges DO count, ie. if all boundaries
   * are identical the size is 1.
   */
  T2 Size () const;
};

template<class T, class T2>
class TemplateRect3D :
  public TemplateRect3D_Base<T,T2>
{
public:
  T2 Size () const { return 0; }
};

template<>
class TemplateRect3D<GridIndexType,int> :
  public TemplateRect3D_Base<GridIndexType,int>
{
public:
  int Size () const
  {
    return (endX-startX+1) * (endY-startY+1) * (endZ-startZ+1);
  }
};

template<>
class TemplateRect3D<float,float> :
  public TemplateRect3D_Base<float,float>
{
public:
  float Size () const
  {
    return (endX-startX) * (endY-startY) * (endZ-startZ);
  }
};

template<>
class TemplateRect3D<double,double> :
  public TemplateRect3D_Base<double,double>
{
public:
  double Size () const
  {
    return (endX-startX) * (endY-startY) * (endZ-startZ);
  }
};

/// Convenience typedef for grid boxes.
typedef TemplateRect3D<GridIndexType,int> Rect3D;

/// Convenience typedef for coordinate boxes.
typedef TemplateRect3D<Types::Coordinate,Types::Coordinate> CoordinateRect3D;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRect3D_h_included_

