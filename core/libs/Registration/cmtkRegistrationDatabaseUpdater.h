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

#ifndef __cmtkRegistrationDatabaseUpdater_h_included_
#define __cmtkRegistrationDatabaseUpdater_h_included_

#include <cmtkconfig.h>

#include "Registration/cmtkImageXformDB.h"

#include <string>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Class for updating the image/transformation database with a newly computed registration.
class RegistrationDatabaseUpdater
{
public:
  /// Constructor.
  RegistrationDatabaseUpdater() {};

  /// Set database path.
  void SetDatabasePath( const std::string& path )
  {
    this->m_DatabasePath = path;
  }

  /// Set input transformation.
  void SetInputXformPath( const std::string& path )
  {
    this->m_InputXformPath = path;
  }

  /// Set output transformation.
  void SetOutputXform( const std::string& path, /**!< Path to the output transformation.*/
		       const bool invertible /**!< Flag whether the output transformation is invertible (i.e., affine). */ )
  {
    this->m_OutputXformPath = path;
    this->m_OutputXformInvertible = invertible;
  }

  /// Update the database.
  void UpdateDB() const;

private:
  /// Path to the database.
  std::string m_DatabasePath;

  /// The input transformation.
  std::string m_InputXformPath;

  /// The output transformation.
  std::string m_OutputXformPath;

  /// Flag for invertible (i.e., affine) output transformation.
  bool m_OutputXformInvertible;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRegistrationDatabaseUpdater_h_included_
