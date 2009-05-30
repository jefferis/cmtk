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

#ifndef __cmtkFusionSlicers_h_included_
#define __cmtkFusionSlicers_h_included_

#include <cmtkconfig.h>

#include <cmtkStudy.h>
#include <cmtkStudyList.h>

#include <cmtkSlicer.h>
#include <cmtkVolumeWrapper.h>
#include <cmtkInterpolator.h>

#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>

#include <map>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

class SlicerPipeline 
{
  igsGetSetMacroObject(Slicer,m_Slicer);

  igsGetSetMacroObject(VolumeWrapper,m_VolumeWrapper);

public:
  SlicerPipeline( Plane *const plane, Study::SmartPtr& study, AffineXform::SmartPtr& affineXform, WarpXform::SmartPtr& warpXform );

  ~SlicerPipeline();

  Image* GetOutput();
};

/// Class that maps a study name to the corresponding Slicer filter.
class FusionSlicers :
  /// This is essentially a map.
  public std::map<Study::SmartPtr,SlicerPipeline*>,
  /// We also want reference counting.
  public Object
{
public:
  /// Location of the reference slice.
  Types::Coordinate ReferenceSlicePosition;

  /// Set location of the reference slice.
  void SetReferenceSlicePosition( const Types::Coordinate referenceSlicePosition )
  {
    ReferenceSlicePosition = referenceSlicePosition;
    this->UpdateSlicePlane();
  }

  /// Slice axis.
  int SliceNormal;

  /// Set location of the reference slice.
  void SetSliceNormal( const int sliceNormal ) 
  {
    SliceNormal = sliceNormal;
    this->UpdateSlicePlane();
  }

  /// Interpolation mode.
  cmtk::Interpolators::InterpolationEnum InterpolationMode;

  /// Set interpolation mode.
  void SetInterpolationMode( const cmtk::Interpolators::InterpolationEnum interpolationMode );

  /// Warp mode.
  bool ApplyWarp;

  /// Set warp mode.
  void SetApplyWarp( const bool applyWarp );

  /// Constructor.
  FusionSlicers();

  /// Virtual destructor.
  virtual ~FusionSlicers();

  /// Select the reference study.
  void SetReferenceStudy( Study::SmartPtr& referenceStudy );

  /// Set the studylist.
  void SetStudyList( StudyList::SmartPtr& studyList );

  /// Get slice plane from given study.
  Image* GetOutput( Study::SmartPtr& study );

  /// Set reslice plane in all stored slicers.
  void SetReslicePlane( const Plane* plane );

private:
  /// The reference study.
  Study::SmartPtr ReferenceStudy;

  /// The complete studylist.
  StudyList::SmartPtr m_StudyList;

  /// The reference slice.
  Image* ReferenceSlice;

  /// Delete all objects referenced in the map and clear map itself.
  void ClearAndDelete();

  /// Update slice plane object after parameter change.
  void UpdateSlicePlane();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFusionSlicers_h_included_
