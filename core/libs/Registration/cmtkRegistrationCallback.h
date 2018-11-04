/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2009, 2013 SRI International
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

#ifndef __cmtkRegistrationCallback_h_included_
#define __cmtkRegistrationCallback_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkVector.h>
#include <Base/cmtkVector3D.h>

#include <System/cmtkSmartPtr.h>

namespace cmtk {

/** \addtogroup Registration */
//@{

/// Status code returned by Execute() methods.
typedef enum {
  /// Everything okay; continue as usual.
  CALLBACK_OK = 0,
  /// User requests interrupt of operation.
  CALLBACK_INTERRUPT = 1,
  /// Interrupt generated by timeout.
  CALLBACK_TIMEOUT = 2,
  /// Something went wrong.
  CALLBACK_FAILED = 3
} CallbackResult;

/** Generic callback class.
 * Callbacks define user specific actions during optimization, i.e. update of
 * progress indicators, protocoling, etc. This particular class is a "dummy"
 * callback, providing a common interface but taking no actions on any
 * member function calls.
 */
class RegistrationCallback {
 public:
  /// SMart pointer.
  typedef SmartPointer<RegistrationCallback> SmartPtr;

  /// Interface: Execute callback action
  virtual CallbackResult ExecuteWithData(const CoordinateVector &v,
                                         const double metric);

  /// Execute callback action without interim result.
  virtual CallbackResult Execute();

  /// Notify callback of an annotation.
  virtual void Comment(const char *comment = NULL);

  /// Default constructor.
  RegistrationCallback();

  /// Virtual destructor.
  virtual ~RegistrationCallback();
};

//@}

}  // namespace cmtk

/// Handler function for SIGINT interrupt signal.
extern "C" void cmtkRegistrationCallbackDispatchSIGINT(int sig);

#endif  // #ifndef __cmtkRegistrationCallback_h_included_
