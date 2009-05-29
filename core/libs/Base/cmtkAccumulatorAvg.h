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

#ifndef __cmtkAccumulatorAvg_h_included_
#define __cmtkAccumulatorAvg_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{
namespace
Accumulators
{
/// Average-value accumulator.
template<class TDataType>
class Avg
{
public:
  /// Default constructor.
  Avg()
  {
    this->m_Sum = 0;
    this->m_Count = 0;
  }

  /// Reset.
  void Reset()
  {
    this->m_Sum = 0;
    this->m_Count = 0;
  }

  /// Add a value.
  void AddValue( const TDataType& value )
  {
    this->m_Sum += value;
    ++this->m_Count;
  }
  
  /// Retrive accumulated result.
  const TDataType GetResult() const
  {
    return this->m_Sum / this->m_Count;
  }

private:
  /// Current sum of values.
  TDataType m_Sum;

  /// Counter.
  unsigned int m_Count;
};
} // namespace Accumulators
} // namespace cmtk

#endif // #ifdef __cmtkAccumulatorAvg_h_included_
