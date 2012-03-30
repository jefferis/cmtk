/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkUniformVolumeLaplaceFilter.h"

#include <Base/cmtkDataGridFilter.h>
#include <Base/cmtkGaussianKernel.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

TypedArray::SmartPtr
UniformVolumeLaplaceFilter::Get() const
{
  // Apply Laplacian symmetric kernel (-1, 2, -1) to each dimension separately.
  std::vector<Types::DataItem> kernel( 2 );
  kernel[0] = 2;
  kernel[1] = -1;

  return DataGridFilter::GetDataKernelFiltered( kernel, kernel, kernel, false /*do not normalize - sum of kernel elements is zero*/ );
}

} // namespace cmtk
