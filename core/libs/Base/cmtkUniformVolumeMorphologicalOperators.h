/*
//
//  Copyright 2004-2013 SRI International
//
//  Copyright 1997-2010 Torsten Rohlfing
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
//  $Revision: 2731 $
//
//  $LastChangedDate: 2011-01-13 16:22:47 -0800 (Thu, 13 Jan 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkUniformVolumeMorphologicalOperators_h_included_
#define __cmtkUniformVolumeMorphologicalOperators_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>
#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Morphological operators acting on a 3D volume with known grid spacings.
 * This class provides implementations of morphological operators that take into account
 * the spacing of grid points (i.e., pixel size). These operators may be preferable to the
 * ones in cmtk::DataGridMorphologicalOperators in cases where the data grid is anisotropic.
 */
class UniformVolumeMorphologicalOperators :
  /// Prevent copying by inheritance.
  private CannotBeCopied 
{
public:
  /// This class.
  typedef UniformVolumeMorphologicalOperators Self;

  /// Constructor: link to UniformVolume object.
  UniformVolumeMorphologicalOperators( const UniformVolume::SmartConstPtr& uniformVolume );

  /// Get data after erosion operator using the Euclidean distance transform.
  TypedArray::SmartPtr GetErodedByDistance( const Types::Coordinate erodeBy /*!< Erosion distance. */ ) const;
  
  /// Get data after dilation operator using the Euclidean distance transform.
  TypedArray::SmartPtr GetDilatedByDistance( const Types::Coordinate dilateBy /*!< Dilation distance. */ ) const;

private:
  /// The UniformVolume object we're working on.
  UniformVolume::SmartConstPtr m_UniformVolume;
};

} // namespace cmtk

#endif // #ifndef __cmtkUniformVolumeMorphologicalOperators_h_included_
