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

#ifndef __cmtkWarpVtkPolyData_h_included_
#define __cmtkWarpVtkPolyData_h_included_

#include <cmtkconfig.h>

#include <vtkPolyDataToPolyDataFilter.h>

#include <cmtkWarpXform.h>
#include <cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

/** Apply deformation to Vtk PolyData object.
 * This class acts as a filter that transforms one vtkPolyData object into
 * another. While the mesh connectivity, ie. the cell information, is
 * not modified, a local coordinate transformation is applied to the node
 * coordinates. This effectively generates a surface that is deformed according
 * to the local transformation in effect.
 */
class WarpVtkPolyData : 
  public vtkPolyDataToPolyDataFilter 
{
public:
  vtkTypeMacro(WarpVtkPolyData,vtkPolyDataToPolyDataFilter);

  /// Create new instance.
  static WarpVtkPolyData *New() 
  { return new WarpVtkPolyData; }

  /// Set deformation to apply.
  void SetWarpXform( WarpXform::SmartPtr& warpXform );

  /// Set affine transformation to apply.
  void SetAffineXform( AffineXform::SmartPtr& affineXform );

  /// Perform actual deformation.
  void Execute();
  
protected:
  /// Default constructor.
  WarpVtkPolyData() : m_WarpXform( NULL ), m_AffineXform( NULL ) {}

  /// Virtual destructor.
  virtual ~WarpVtkPolyData() {}

private:
  /// The deformation to apply.
  WarpXform::SmartPtr m_WarpXform;

  /// We may also want to support affine transformations.
  AffineXform::SmartPtr m_AffineXform;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkWarpVtkPolyData_h_included_
