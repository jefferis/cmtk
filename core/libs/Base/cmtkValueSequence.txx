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

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T>
ValueSequence<T>& 
ValueSequence<T>::operator=( const ValueSequence<T>& other )
{
  this->NValues = other.NValues;
  this->Sum = other.Sum; 
  this->SumAbs = other.SumAbs; 
  this->SumOfSquares = other.SumOfSquares; 
  this->Minimum = other.Minimum; 
  this->Maximum = other.Maximum;
  this->MinimumAbs = other.MinimumAbs;
  this->MaximumAbs = other.MaximumAbs; 
  
  return *this;
}

template<class T>
ValueSequence<T> operator+( const ValueSequence<T>& a, const ValueSequence<T>& b )
{
  ValueSequence<T> result;

  result.NValues = a.NValues + b.NValues;
  result.Sum = a.Sum + b.Sum; 
  result.SumAbs = a.SumAbs + b.SumAbs; 
  result.SumOfSquares = a.SumOfSquares + b.SumOfSquares; 
  result.Minimum = std::min( a.Minimum, b.Minimum );
  result.Maximum = std::max( a.Maximum, b.Maximum );
  result.MinimumAbs = std::min( a.MinimumAbs, b.MinimumAbs );
  result.MaximumAbs = std::max( a.MaximumAbs, b.MaximumAbs ); 

  return result;
}

} // namespace cmtk
