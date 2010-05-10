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

#include <cmtkXform.h>

#include <cmtkVolume.h>

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
    this->m_ParameterVector = CoordinateVector::SmartPtr::Null;
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

void
Xform::GetVolumeOfInfluence
( const size_t, const Self::SpaceVectorType& fromVol, const Self::SpaceVectorType& toVol, Self::SpaceVectorType& fromVOI, Self::SpaceVectorType& toVOI, const int ) const
{
  fromVOI = fromVol;
  toVOI = toVol;
}

Types::Coordinate
Xform::GetLandmarksMSD( const MatchedLandmarkList* ll ) const
{
  double MSD = 0;

  MatchedLandmarkList::const_iterator it = ll->begin();
  while ( it != ll->end() )
    {
    Self::SpaceVectorType source( (*it)->GetLocation() );
    Self::SpaceVectorType target( (*it)->GetTargetLocation() );
    this->ApplyInPlace( source );
    MSD += (source - target).SumOfSquares();
    ++it;
    }
  
  MSD /= ll->size();

  return MSD;
}

} // namespace cmtk
