/*
//
//  Copyright 2011 SRI International
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

#ifndef __cmtkDataGridLocalCorrelation_h_included_
#define __cmtkDataGridLocalCorrelation_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkDataGrid.h>
#include <Base/cmtkFixedVector.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Compute local correlation between two data grid objects.
 */
class DataGridLocalCorrelation
{
public:
  /// This class.
  typedef DataGridLocalCorrelation Self;

  /// Constructor.
  DataGridLocalCorrelation( const DataGrid& dg1 /*!< First input data grid. */, const DataGrid& dg2 /*!< Second input data grid */ ) : m_Grid1( dg1 ), m_Grid2( dg2 ), m_Result( NULL ) {}

  /// Get result.
  DataGrid::SmartPtr& GetResult()
  {
    if ( ! this->m_Result )
      this->ComputeResult();

    return this->m_Result;
  }

private:
  /// First input data grid.
  const DataGrid& m_Grid1;

  /// Second input data grid.
  const DataGrid& m_Grid2;

  /// Window radius.
  DataGrid::IndexType m_Radius;

  /// Result data grid.
  DataGrid::SmartPtr m_Result;

  /// Compute the result using the current inputs and parameters.
  void ComputeResult();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDataGridLocalCorrelation_h_included_
