/*
//
//  Copyright 2012 SRI International
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

#ifndef __cmtkRegressionTracker_h_included_
#define __cmtkRegressionTracker_h_included_

#include <cmtkconfig.h>

#include <System/cmtkConsole.h>

#include <stdio.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Class for generating and following regression tracks.
class RegressionTracker
{
public:
  /// This class.
  typedef RegressionTracker Self;

  /// Constructor.
  RegressionTracker();

  /// Virtual destructor.
  virtual ~RegressionTracker();

  /// Instantiate and return static instance.
  static Self& Static()
  {
    static Self tracker;
    return tracker;
  }

  /// Advance tracker: compute checksum and either write to baseline file or compare against existing file.
  template<class T> void Advance( const T* data, const size_t nBytes )
  {
    this->CompareChecksum( reinterpret_cast<const unsigned char*>( data ), nBytes );
  }

private:
  /// Regression trace file.
  FILE* m_File;

  /// Flag for writing.
  bool m_WriteFlag;

  /// Compare using checksum
  void CompareChecksum( const unsigned char *const data, size_t nBytes );

protected:
  /// Member fuction that is called when divergence from previous trace is detected.
  virtual void Trap()
  {
    StdErr << "Detected regression divergence.\n";
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRegressionTracker_h_included_
