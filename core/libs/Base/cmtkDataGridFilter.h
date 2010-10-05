/*
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkDataGridFilter_h_included_
#define __cmtkDataGridFilter_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>
#include <Base/cmtkDataGrid.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Filter operations for data on 3D grids.
 * The filters in this class operate independent of grid spacings.
 */
class DataGridFilter :
  /// Prevent copying by inheritance.
  private CannotBeCopied 
{
public:
  /// This class.
  typedef DataGridFilter Self;

  /// Constructor: link to DataGrid object.
  explicit DataGridFilter( DataGrid::SmartPtr dataGrid );

  /// Return data after median-filtering with global filter radius (convenience function).
  TypedArray::SmartPtr GetDataMedianFiltered( const int radius ) const
  {
    return this->GetDataMedianFiltered( radius, radius, radius );
  }

  /// Return data after median-filtering with per-dimension filter radius
  TypedArray::SmartPtr GetDataMedianFiltered( const int radiusX, const int radiusY, const int radiusZ ) const;
  
  /// Return data after median-filtering.
  TypedArray::SmartPtr GetDataSobelFiltered() const;

  /// Return after filtering with a separable kernel
  TypedArray::SmartPtr GetDataKernelFiltered( const std::vector<Types::DataItem>& filterX, const std::vector<Types::DataItem>& filterY, const std::vector<Types::DataItem>& filterZ ) const;

private:
  /// The DataGrid object we're working on.
  DataGrid::SmartPtr m_DataGrid;

  /// Thread parameter for entropy evaluation.
  class FilterThreadParameters : 
    /// Inherit from generic thread parameter class.
    public ThreadParameters<const Self>
  {
  public:
    /// Filter kernel.
    const std::vector<Types::DataItem>* m_Filter;

    /// Pointer to result pixel data
    TypedArray::SmartPtr m_Result;
  };
  
  /// Thread function for separable filtering in x-direction.
  static void GetFilteredDataThreadX( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Thread function for separable filtering in y-direction.
  static void GetFilteredDataThreadY( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Thread function for separable filtering in z-direction.
  static void GetFilteredDataThreadZ( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

};

} // namespace cmtk

#endif // #ifndef __cmtkDataGridFilter_h_included_
