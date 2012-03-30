/*
//
//  Copyright 2004-2012 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkUniformVolumeLaplaceFilter_h_included_
#define __cmtkUniformVolumeLaplaceFilter_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkDataGridFilter.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUnits.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Laplacian edge enhancement filter for uniform 3D volume images.
 */
class UniformVolumeLaplaceFilter :
    /// Prevent copying by inheritance.
  private DataGridFilter
{
public:
  /// This class.
  typedef UniformVolumeLaplaceFilter Self;

  /// Constructor: link to UniformVolume object.
  explicit UniformVolumeLaplaceFilter( UniformVolume::SmartConstPtr& volume ) : DataGridFilter( volume ), m_UniformVolume( volume ) {}

  /// Laplacian (edge enhancing) filter.
  TypedArray::SmartPtr Get() const;

private:
  /// The UniformVolume object we're working on.
  UniformVolume::SmartConstPtr m_UniformVolume;
};

} // namespace cmtk

#endif // #ifndef __cmtkUniformVolumeLaplaceFilter_h_included_
