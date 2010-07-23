/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkDeviceArrayCUDA_h_included_
#define __cmtkDeviceArrayCUDA_h_included_

#include "System/cmtkCannotBeCopied.h"
#include "System/cmtkSmartConstPtr.h"
#include "System/cmtkSmartPtr.h"

namespace
cmtk
{

/// Forward declaration.
struct cudaArray;

/** \addtogroup GPU */
//@{

/// Resource managing class for raw memory allocated on a GPU device through CUDA.
class DeviceArrayCUDA
    /// Make sure this is never copied.
  : private CannotBeCopied
{
public:
  /// This class.
  typedef DeviceArrayCUDA Self;

  /// Smart pointer-to-const.
  typedef SmartConstPointer<Self> SmartConstPtr;
  
  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;
  
  /// Exception for failed allocation.
  class bad_alloc : public std::bad_alloc {};
  
  /// Destructor: free memory through CUDA.
  virtual ~DeviceArrayCUDA();

private:
  /// Opaque pointer to array on device.
  struct cudaArray* m_DeviceArrayPtr;
};

} // namespace cmtk

#endif // #ifndef __cmtkDeviceArrayCUDA_h_included_
