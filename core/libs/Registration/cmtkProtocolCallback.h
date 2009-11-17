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

#ifndef __cmtkProtocolCallback_h_included_
#define __cmtkProtocolCallback_h_included_

#include <cmtkconfig.h>

#include <cmtkRegistrationCallback.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Callback object with protocol functionality.
class ProtocolCallback : 
  public RegistrationCallback 
{
public:
  /// Open protocol stream if filename is given.
  ProtocolCallback ( const char *filename = NULL, const bool debug = false );

  /** Destructor. 
   * Closes protocol stream.
   */
  virtual ~ProtocolCallback ();

  /// Get virtual class name.
  virtual const char* GetClassName() const { return "ProtocolCallback"; }

  /// Execute callback action.
  virtual CallbackResult ExecuteWithData ( const CoordinateVector& v, const double metric );

  /// Write comment to protocol file.
  virtual void Comment ( const char* comment = NULL );

private:
  /// Protocol stream.
  FILE *fp;

  /// Debug flag.
  bool Debug;

  /// Convenience type definition.
  typedef RegistrationCallback Superclass;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkProtocolCallback_h_included_
