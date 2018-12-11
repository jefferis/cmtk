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

#ifndef __cmtkImageSymmetryPlaneCommandLine_h_included_
#define __cmtkImageSymmetryPlaneCommandLine_h_included_

#include <cmtkconfig.h>

#include "cmtkImageSymmetryPlaneCommandLineBase.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class template for symmetry plane computation command line tools.
 * This is templated over the symmetry plane functional, which can be either
 * CPU-based or GPU-based.
 */
template<class TFunctional>
class ImageSymmetryPlaneCommandLine :
  /// Inherit from non-template base class.
  public ImageSymmetryPlaneCommandLineBase
{
public:
  /// The functional type.
  typedef TFunctional FunctionalType;

  /// This class.
  typedef ImageSymmetryPlaneCommandLine<TFunctional> Self;

  /// Parent class.
  typedef ImageSymmetryPlaneCommandLineBase Superclass;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const.
  typedef SmartConstPointer<Self> SmartConstPtr;

protected:
  /// Make a functional of the given template type using an image volume.
  virtual ImageSymmetryPlaneFunctionalBase::SmartPtr CreateFunctional( UniformVolume::SmartPtr& volume )
  {
    return ImageSymmetryPlaneFunctionalBase::SmartPtr( new FunctionalType( volume ) );
  }

  /// Make a functional of the given template type using an image volume and a value range.
  virtual ImageSymmetryPlaneFunctionalBase::SmartPtr CreateFunctional( UniformVolume::SmartPtr& volume, const cmtk::Types::DataItemRange& range )
  {
    return ImageSymmetryPlaneFunctionalBase::SmartPtr( new FunctionalType( volume, range ) );
  }
};

} // namespace cmtk

#endif // #ifndef  __cmtkImageSymmetryPlaneCommandLine_h_included_
