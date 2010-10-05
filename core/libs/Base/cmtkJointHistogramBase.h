/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkJointHistogramBase_h_included_
#define __cmtkJointHistogramBase_h_included_

#include <cmtkconfig.h>

#include <System/cmtkSmartPtr.h>
#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Basic (non-template) 2-D histogram functions.
class JointHistogramBase
{
public:
  /// This class.
  typedef JointHistogramBase Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Calculate optimum number of histogram bins.
  static size_t CalcNumBins( const size_t numberOfSamples /*!< Number of data values. */, const Types::DataItemRange& valueRange /*!< Range of values in the data.*/ );

  /// Calculate optimum number of histogram bins for given volume.
  static size_t CalcNumBins ( const UniformVolume* volume );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkJointHistogramBase_h_included_
