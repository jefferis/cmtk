/*
//
//  Copyright 2009 SRI International
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
//  $Revision: 389 $
//
//  $LastChangedDate: 2009-08-05 11:33:03 -0700 (Wed, 05 Aug 2009) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkSplineWarpXformITKIO_h_included__
#define __cmtkSplineWarpXformITKIO_h_included__

#include <cmtkconfig.h>

#include <cmtkSplineWarpXform.h>

#include <string>

namespace
cmtk
{

/** Class for reading and writing affine transformations from and to ITK's file format.
 * This should also be understood by Slicer3 for transformation exchange with CMTK tools
 * run as plugins.
 */
class SplineWarpXformITKIO
{
public:
  /// Write transformation to ITK file.
  static void Write( const std::string& filename, const SplineWarpXform& xform );
};

} // namespace cmtk

#endif // #ifndef __cmtkSplineWarpXformITKIO_h_included__
