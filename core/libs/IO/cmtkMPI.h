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

#ifndef __cmtkMPI_h_included_
#define __cmtkMPI_h_included_

#include <cmtkconfig.h>

#ifdef CMTK_BUILD_MPI

#include <mpi.h>

#include <cmtkSmartPtr.h>

namespace cmtk
{

/** \addtogroup IO */
//@{

namespace mpi
{

/// Broadcast an object from a root process to all others.
template<class TClass>
void
Broadcast( MPI::Intracomm& comm, SmartPointer<TClass>& inOutPtr, const int root );

} // namespace MPI

//@}

} // namespace cmtk

#endif // #ifdef CMTK_BUILD_MPI

#endif // #ifndef __cmtkMPI_h_included_
