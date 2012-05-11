/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkClassStreamOutput_h_included_
#define __cmtkClassStreamOutput_h_included_

#include <cmtkconfig.h>

#include <IO/cmtkTypedStreamOutput.h>
#include <IO/cmtkStudy.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkParametricPlane.h>

#include <string>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Class for writing various library classes to and from disk.
 */
class ClassStreamOutput : 
  /// Inherit basic functionality from typed stream.
  public TypedStreamOutput
{
public:
  /// This class.
  typedef ClassStreamOutput Self;

  /// Parent class.
  typedef TypedStreamOutput Superclass;

  /// Default constructor.
  ClassStreamOutput() : TypedStreamOutput() {}

  /** Open constructor.
   *\param filename Name of the archive to open.
   */
  ClassStreamOutput( const std::string& filename, const Superclass::Mode mode ) : TypedStreamOutput( filename, mode ) {}

  /** Open constructor for separate path and archive names.
   *\param dir Directory to open archive in.
   *\param archive Name of the archive to open.
   */
  ClassStreamOutput( const std::string& dir, const std::string& archive, const Superclass::Mode mode ) : TypedStreamOutput( dir, archive, mode ) {}

  /** Write generic transformation object.
   * This function determines the virtual type of the transformation object
   * (spline or linear deformation) using a dynamic_cast. It then calls the
   * appropriate specialized output function.
   */
  ClassStreamOutput& operator << ( const WarpXform *warpXform );

  /** Write spline transformation object.
   * This function works on a reference rather than a pointer. It immediately
   * calls the pointer-based function defined above for the actual writing.
   */
  ClassStreamOutput& operator << ( const SplineWarpXform& splineWarpXform )
  { return (*this) << &splineWarpXform; }
  
  /** Write parametric plane object.
   */
  ClassStreamOutput& operator << ( const ParametricPlane *parametricPlane );

  /** Write parametric plane object.
   * This function works on a reference rather than a pointer. It immediately
   * calls the pointer-based function defined above for the actual writing.
   */
  ClassStreamOutput& operator << ( const ParametricPlane& parametricPlane )
  { return (*this) << &parametricPlane; }

private:
  /// Write actual warp transformation object.
  ClassStreamOutput& PutWarp( const WarpXform* warpXform  );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkClassStreamOutput_h_included_
