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

#ifndef __cmtkDirectionSetOptimizer_h_included_
#define __cmtkDirectionSetOptimizer_h_included_

#include <cmtkconfig.h>

#include <cmtkOptimizer.h>
#include <cmtkDirectionSet.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Optimizer that works along a pre-defined set of search dimensions.
 * This class is particularly useful and intended for optimization during
 * a norigid registration according to the principal modes of variation in
 * an active deformation model.
 */
class DirectionSetOptimizer :
  /// Inherit interface from generic optimizer.
  public Optimizer
{
public:
  /// This class.
  typedef DirectionSetOptimizer Self;
  
  /// Superclass.
  typedef Optimizer Superclass;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;
  
  /// Direction set.
  cmtkGetSetMacro(DirectionSet::SmartPtr,DirectionSet);

  /// Optimize target functional.
  virtual CallbackResult Optimize( CoordinateVector& v, const Self::ParameterType exploration, const Self::ParameterType accuracy );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDirectionSetOptimizer_h_included_
