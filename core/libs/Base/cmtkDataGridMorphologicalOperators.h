/*
//
//  Copyright 2004-2010 SRI International
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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkDataGridMorphologicalOperators_h_included_
#define __cmtkDataGridMorphologicalOperators_h_included_

#include <cmtkconfig.h>

#include "System/cmtkCannotBeCopied.h"
#include "Base/cmtkDataGrid.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Grid topology of data arranged in a 3D lattice.
 * This class extends the plain data handling functions of TypedArray
 * with a 3D topology. Real world coordinates, however, are not considered and
 * need to be handled by derived classes. Thus, this class provides the coordinate
 * independent services such as median filtering and, to a certain extent,
 * interpolation.
 */
class DataGridMorphologicalOperators :
  /// Prevent copying by inheritance.
  private CannotBeCopied 
{
public:
  /// This class.
  typedef DataGridMorphologicalOperators Self;

  /// Constructor: link to DataGrid object.
  DataGridMorphologicalOperators( const DataGrid::SmartConstPtr& dataGrid );

  /** Eliminating padding data by neighborhood voting.
   *\return Returns "true" if data was actually changed, "false" if no change
   *  was necessary or possible.
   */
  bool EliminatePaddingVoting( const int iterations = 1 /*!< Number of elimination iterations.*/ );

  /** Return map of region boundaries.
   * This function returns a byte data array where each pixel is one if it is
   * a boundary pixel, i.e., if one of its neighbours in this object has a
   * different value than it has itself. All other pixels are set to zero.
   *\param multiValue If this is set (default: false), then the resulting
   *  boundary map is multi valued, i.e., instead of setting boundary pixels
   *  to "1", they are set to the value present in the image at that location.
   *\note The boundary contours are at least 2 pixels wide since "boundaryness"
   * is a symmetric relationship.
   */
  TypedArray::SmartPtr GetBoundaryMap( const bool multiValued = false ) const;

  /// Get data after erosion operator.
  TypedArray::SmartPtr GetEroded( const int iterations = 1 /*!< Number of erosion iterations. */ ) const;
  
  /// Get data after dilation operator.
  TypedArray::SmartPtr GetDilated( const int iterations = 1 /*!< Number of dilation iterations. */ ) const;

  /** Get connected components of a binary image.
   * All pixels with non-zero values are considered "foreground," and the result
   * of this function is a partitioning of the foreground into connected components.
   * Connectivity is determined based on 8 neighbours in the 3D grid.
   */
  TypedArray::SmartPtr GetBinaryConnectedComponents() const;

  /** Get data with region labels renumbered by decreasing region size.
   * This is typically run after GetBinaryConnectedComponents() to order the
   * components and be able to easily select, say, the largest k components.
   */
  TypedArray::SmartPtr GetRegionsRenumberedBySize() const;

private:
  /// The DataGrid object we're working on.
  DataGrid::SmartConstPtr m_DataGrid;
};

} // namespace cmtk

#endif // #ifndef __cmtkDataGridMorphologicalOperators_h_included_
