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

#include "cmtkBitVector.h"

#include <System/cmtkMemory.h>

#include <memory.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

BitVector::BitVector( const size_t size, const bool initial )
{
  this->m_Size = (size+7) / 8; // +7 to allocate an extra byte for 8n+1...8n+7 bits
  this->m_BitVector = Memory::AllocateArray<byte>( this->m_Size );

  if ( initial )
    this->Set();
  else
    this->Reset();
}

BitVector::BitVector( const size_t size, byte *const bitSet )
{
  this->m_Size = (size+7) / 8; // +7 to allocate an extra byte for 8n+1...8n+7 bits
  this->m_BitVector = bitSet;
}

BitVector::~BitVector()
{
  delete[] this->m_BitVector;
}

BitVector* 
BitVector::Clone() const
{
  byte *newBitVector = Memory::AllocateArray<byte>( this->m_Size );
  memcpy( newBitVector, this->m_BitVector, this->m_Size );
  return new BitVector( 8*this->m_Size, newBitVector );
}

void
BitVector::Set()
{
  memset( this->m_BitVector, 255, sizeof( *this->m_BitVector ) * this->m_Size );
}

void
BitVector::Set( const size_t pos, const bool val )
{
  if ( val )
    {
    this->m_BitVector[pos/8] |= (1<<(pos%8));
    } 
  else
    {
    this->m_BitVector[pos/8] &= ~(1<<(pos%8));
    }
}

void
BitVector::Reset( const bool value )
{
  if ( value )
    memset( this->m_BitVector, 255, sizeof( *this->m_BitVector ) * this->m_Size );
  else
    memset( this->m_BitVector, 0, sizeof( *this->m_BitVector ) * this->m_Size );
}

void
BitVector::Reset( const size_t pos )
{
  this->m_BitVector[pos/8] &= ~(1<<(pos%8));
}

void
BitVector::Flip()
{
  for ( size_t i=0; i<this->m_Size; ++i )
    this->m_BitVector[i] = ~this->m_BitVector[i];
}

void
BitVector::Flip( const size_t pos )
{
  this->m_BitVector[pos/8] ^= (1<<(pos%8));
}

bool
BitVector::operator[]( const size_t pos ) const
{
  return ( this->m_BitVector[pos/8] >> (pos%8) ) & 1;
}

} // namespace cmtk
