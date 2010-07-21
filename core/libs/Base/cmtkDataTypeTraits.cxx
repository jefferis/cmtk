/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkDataTypeTraits.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

const Types::DataItem DataTypeTraits<byte>::Min = 0;
const Types::DataItem DataTypeTraits<byte>::Max = UCHAR_MAX;

const Types::DataItem DataTypeTraits<char>::Min = SCHAR_MIN;
const Types::DataItem DataTypeTraits<char>::Max = SCHAR_MAX;

const Types::DataItem DataTypeTraits<signed short>::Min = SHRT_MIN;
const Types::DataItem DataTypeTraits<signed short>::Max = SHRT_MAX;

const Types::DataItem DataTypeTraits<unsigned short>::Min = 0;
const Types::DataItem DataTypeTraits<unsigned short>::Max = USHRT_MAX;

const Types::DataItem DataTypeTraits<int>::Min = INT_MIN;
const Types::DataItem DataTypeTraits<int>::Max = INT_MAX;

} // namespace cmtk
