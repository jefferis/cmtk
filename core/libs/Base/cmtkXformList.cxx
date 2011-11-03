/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include "cmtkXformList.h"

void
cmtk::XformList::Add
( const Xform::SmartConstPtr& xform, const bool inverse, const Types::Coordinate globalScale  )
{
  this->push_back( XformListEntry::SmartConstPtr( new XformListEntry( xform, inverse, globalScale ) ) );
}

void
cmtk::XformList::AddToFront
( const Xform::SmartConstPtr& xform, const bool inverse, const Types::Coordinate globalScale  )
{
  this->push_front( XformListEntry::SmartConstPtr( new XformListEntry( xform, inverse, globalScale ) ) );
}

bool
cmtk::XformList::ApplyInPlace( Xform::SpaceVectorType& v ) const
{
  for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
    {
    if ( (*it)->Inverse ) 
      {
      if ( (*it)->m_WarpXform ) 
	{
	// yes: use approximate inverse
	if ( ! (*it)->m_WarpXform->ApplyInverseInPlace( v, this->m_Epsilon ) ) 
	  {
	  // if that fails, return failure flag
	  return false;
	  }
	} 
      else
	{
	// is this an affine transformation that has an inverse?
	if ( (*it)->InverseAffineXform ) 
	  // apply inverse
	  (*it)->InverseAffineXform->ApplyInPlace( v );
	else
	  // nothing else we can do: exit with failure flag
	  return false;
	}
      } 
    else
      {
      // are we outside xform domain? then return failure.
      if ( !(*it)->m_Xform->InDomain( v ) ) return false;
      (*it)->m_Xform->ApplyInPlace( v );
      }
    }
  return true;
}

bool
cmtk::XformList::GetJacobian
( const Xform::SpaceVectorType& v, Types::DataItem& jacobian, const bool correctGlobalScale ) const
{
  Xform::SpaceVectorType vv( v );

  jacobian = static_cast<Types::DataItem>( 1.0 );
  for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
    {
    if ( (*it)->Inverse ) 
      {
      if ( correctGlobalScale )
	jacobian *= static_cast<Types::DataItem>( (*it)->GlobalScale );

      // is this a spline transformation?
      if ( (*it)->m_WarpXform ) 
	{
	// yes: use approximate inverse
	if ( (*it)->m_WarpXform->ApplyInverseInPlace( vv, this->m_Epsilon ) ) 
	  // compute Jacobian at destination and invert
	  jacobian /= static_cast<Types::DataItem>( (*it)->m_Xform->GetJacobianDeterminant( vv ) );	
	else
	  // if that fails, return failure flag
	  return false;
	} 
      else
	{
	// is this an affine transformation that has an inverse?
	if ( (*it)->InverseAffineXform ) 
	  // apply inverse
	  (*it)->InverseAffineXform->ApplyInPlace( vv );
	else
	  // nothing else we can do: exit with failure flag
	  return false;
	}
      } 
    else 
      {
      // are we outside xform domain? then return failure.
      if ( !(*it)->m_Xform->InDomain( v ) ) return false;

      jacobian *= static_cast<Types::DataItem>( (*it)->m_Xform->GetJacobianDeterminant( vv ) );
      if ( correctGlobalScale )
	jacobian /= static_cast<Types::DataItem>( (*it)->GlobalScale );
      (*it)->m_Xform->ApplyInPlace( vv );
      }
    }
  return true;
}

bool
cmtk::XformList::AllAffine() const
{
  for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
    {
    if ( !(*it)->IsAffine() )
      return false;
    }
  return true;
}

cmtk::XformList
cmtk::XformList::MakeAllAffine() const
{
  cmtk::XformList allAffine;

  for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
    {
    allAffine.push_back( (*it)->CopyAsAffine() );
    }

  return allAffine;
}

std::string
cmtk::XformList::GetFixedImagePath() const
{
  const XformListEntry& first = **(this->begin());

  // if transformation is inverse, get original "moving" path instead.
  if ( first.Inverse )
    return first.m_Xform->GetMetaInfo( META_XFORM_MOVING_IMAGE_PATH, "" );
  else
    return first.m_Xform->GetMetaInfo( META_XFORM_FIXED_IMAGE_PATH, "" );
}

std::string
cmtk::XformList::GetMovingImagePath() const
{
  const XformListEntry& last = **(this->rbegin());
  
  // if transformation is inverse, get original "fixed" path instead.
  if ( last.Inverse )
    return last.m_Xform->GetMetaInfo( META_XFORM_FIXED_IMAGE_PATH, "" );
  else
    return last.m_Xform->GetMetaInfo( META_XFORM_MOVING_IMAGE_PATH, "" );
}
