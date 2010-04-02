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

#ifndef __cmtkCannotBeCopied_h_included_
#define __cmtkCannotBeCopied_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Base class to prevent copying of derived classes.
 * This class can be inherited by derived classes and will thus prevent
 * the compiler from generating default copy constructor and assignment 
 * operator (see Effective C++, 3rd ed.)
 */
class CannotBeCopied 
{
protected:
  /// Default constructor.
  CannotBeCopied() {};

  /// Default destructor.
  ~CannotBeCopied() {};

private:
  /// Undefined copy constructor.
  CannotBeCopied( const CannotBeCopied& );

  /// Undefined assignment operator.
  CannotBeCopied& operator=( const CannotBeCopied& );
};

}

#endif // #ifndef __cmtkCannotBeCopied_h_included_
