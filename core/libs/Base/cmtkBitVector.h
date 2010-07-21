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

#ifndef __cmtkBitVector_h_included_
#define __cmtkBitVector_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkTypes.h"
#include "System/cmtkSmartPtr.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Set of binary values.
 * This class provides functions similar to the STL's bitset class. However,
 * our class does not require the set size to be known at compile time. It
 * therefore allows creation of different size bitsets as they are needed by
 * the program.
 */
class BitVector 
{
public:
  /// Smart pointer to BitVector.
  typedef SmartPointer<BitVector> SmartPtr;

  /** Constructor.
   * The bitset is allocated, but NOT initialized.
   *@param size Number of bits handled by this object.
   */
  BitVector( const size_t size, const bool initial = false );
  
  /** Constructor.
   * The bitset is allocated, but NOT initialized.
   *@param size Number of bits handled by this object.
   */
  BitVector( const size_t size, byte *const bitset );
  
  /** Destructor.
   */
  ~BitVector();

  /** Create copy of this object.
   */
  BitVector* Clone() const;

  /// Set all bits to 1.
  void Set();

  /// Set one bit to a given value.
  void Set( const size_t pos, const bool val = true );

  /// Set all bits to given flag (default: clear all).
  void Reset( const bool value = false );

  /// Set one bit to 0.
  void Reset( const size_t pos );

  /// Flip (invert) the whole bitset.
  void Flip();

  /// Flip (invert) one bit.
  void Flip( const size_t pos );

  /// Return a given bit.
  bool operator[]( const size_t pos ) const;

  /// Get pointer to bitset data.
  const byte* GetBitVector() const
  { 
    return this->m_BitVector; 
  }

private:
  /// The bitset.
  byte *m_BitVector;

  /// The size of the allocated bitset in BYTES (!!).
  size_t m_Size;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkBitVector_h_included_
