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

#ifndef __cmtkSubtractFilter_h_included_
#define __cmtkSubtractFilter_h_included_

#include <cmtkArrayFilter.h>

#include <cmtkImage.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Filter class for computation of subtraction images.
 */
class SubtractFilter : 
  public ArrayFilter<Image,Image,2> 
{
public:
  /// Create new object.
  static SubtractFilter* New() { return new SubtractFilter; }

  /// Return virtual class name.
  virtual const char* GetClassName() const { return "SubtractFilter"; }

  /// Execute subtraction operator.
  virtual void Execute();

  /// Index of image to subtract from.
  igsClassParameter(int,SubtractImageIndex);

  /// Flag for absolute value subtraction.
  igsClassParameter(bool,Absolute);
  
  /// Flag for threshold application.
  igsClassParameter(bool,UseThresholds);

  /// Thresholding range.
  igsClassParameter2Array(float,Thresholds);
  
protected:
  /// Default constructor.
  SubtractFilter() {}

  /// Destructor.
  virtual ~SubtractFilter() {}
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSubtractFilter_h_included_
