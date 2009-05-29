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

#ifndef __cmtkImageEdgeOperator_h_included_
#define __cmtkImageEdgeOperator_h_included_

#include <cmtkconfig.h>

#include <cmtkFilter.h>
#include <cmtkImage.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Filter class to apply edge operators to Image objects.
 * As a result, a new Image object of identical geometrical structure is 
 * created. Its primitive data type is "float" as some edge operators may
 * produce real-valued data even on when operating on integer values.
 */
class ImageEdgeOperator : 
  /// This is a filter from a scalar image to another.
  public Filter<Image,Image> 
{
public:
  /// Create new object.
  static ImageEdgeOperator* New() { return new ImageEdgeOperator; }

  /// Return virtual class name.
  virtual const char* GetClassName() const { return "ImageEdgeOperator"; }

  /// Execute edge operator.
  virtual void Execute();

  /** Flag selecting which edge operator to use.
   * Supported values are: 0 - Laplacian, Sum of absolute orthogonal Sobel 
   * operators.
   */
  igsClassParameter(int,Operator);

  /** Smoothing flag.
   * If this flag is set, a Gaussian smoothing operator is applied to the image
   * before edge computation.
   */
  igsClassParameter(bool,SmoothBeforeEdge);

  /** Gaussian smoothing kernel width (standard deviation in pixel units).
   */
  igsClassParameter(float,GaussianWidth);

protected:
  /// Default constructor.
  ImageEdgeOperator() 
  { Operator = 0; SmoothBeforeEdge = 0; GaussianWidth = 1; }

  /// Destructor.
  virtual ~ImageEdgeOperator() {}
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageEdgeOperator_h_included_
