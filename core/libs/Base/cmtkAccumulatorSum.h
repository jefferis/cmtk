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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkAccumulatorSum_h_included_
#define __cmtkAccumulatorSum_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Helper classes that project a series of values into a single value.
namespace
Accumulators
{
/// Sum accumulator.
template<class TDataType>
class Sum
{
public:
  /// Default constructor.
  Sum()
  {
    this->m_Sum = 0;
  }

  /// Reset.
  void Reset()
  {
    this->m_Sum = 0;
  }

  /// Add a value.
  void AddValue( const TDataType& value )
  {
    this->m_Sum += value;
  }
  
  /// Retrive accumulated result.
  const TDataType& GetResult() const
  {
    return this->m_Sum;
  }

private:
  /// Current sum
  TDataType m_Sum;
};
} // namespace Accumulators
} // namespace cmtk

#endif // #ifdef __cmtkAccumulatorSum_h_included_
