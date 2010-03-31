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

#ifndef __cmtkFusionHV_h_included_
#define __cmtkFusionHV_h_included_

#include "cmtkArrayFilter.h"

#include "cmtkImageRGB.h"
#include "cmtkMacros.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/// Class for color/brightness encoding image fusion.
class FusionHV : 
    public ArrayFilter<ImageRGB,ImageRGB,2> 
{
public:
  /// Create new class instance.
  static FusionHV* New() { return new FusionHV; }

  /// Parameter selecting the weighting of both modalities.
  igsClassParameter(double,Contrast);

  /// Parameter selecting the weighting of both modalities.
  igsClassParameter(int,ImageIndexV);

  /// Perform fusion.
  virtual void Execute();

protected:
  /// Constructor.
  FusionHV();

  /// Virtual destructor.
  virtual ~FusionHV() {};
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFusionHV_h_included_
