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

#include "cmtkXform.h"

#include <Base/cmtkVolume.h>
#include <Base/cmtkLandmarkPairList.h>

#include <math.h>
#include <stdio.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
Xform::AllocateParameterVector
( const size_t numberOfParameters )
{
  this->m_NumberOfParameters = numberOfParameters;
  if ( this->m_NumberOfParameters )
    {
    this->m_ParameterVector = CoordinateVector::SmartPtr( new CoordinateVector( this->m_NumberOfParameters ) );
    this->m_Parameters = this->m_ParameterVector->Elements;
    }
  else
    {
    this->m_ParameterVector = CoordinateVector::SmartPtr::Null();
    this->m_Parameters = NULL;
    }
}

void
Xform::SetParamVector ( CoordinateVector& v ) 
{
  if ( this->m_ParameterVector ) 
    {
    *this->m_ParameterVector = v;
    } 
  else
    {
    this->m_ParameterVector = CoordinateVector::SmartPtr( new CoordinateVector( v ) );
    }
  this->m_Parameters = this->m_ParameterVector->Elements;
}

void
Xform::SetParamVector ( const CoordinateVector& v ) 
{
  if ( this->m_ParameterVector ) 
    {
    *this->m_ParameterVector = v;
    } 
  else
    {
    this->m_ParameterVector = CoordinateVector::SmartPtr( new CoordinateVector( v ) );
    }
  this->m_Parameters = this->m_ParameterVector->Elements;
}

CoordinateVector& 
Xform::GetParamVector
( CoordinateVector& v, const size_t targetOffset ) const 
{
  v.AdjustDimension( std::max<int>( v.Dim, targetOffset + this->ParamVectorDim() ) );
  v.CopyToOffset( *this->m_ParameterVector, targetOffset, this->ParamVectorDim() );
  return v;
}

Types::Coordinate
Xform::GetLandmarksMSD( const LandmarkPairList& ll ) const
{
  Types::Coordinate msd = 0;

  for ( LandmarkPairList::const_iterator it = ll.begin(); it != ll.end(); ++it )
    {
    msd += ( this->Apply( it->m_Location ) - it->m_TargetLocation ).SumOfSquares();
    }
  
  return msd / ll.size();
}

} // namespace cmtk
