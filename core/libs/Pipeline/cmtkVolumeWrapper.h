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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkVolumeWrapper_h_included_
#define __cmtkVolumeWrapper_h_included_

#include <cmtkconfig.h>

#include <cmtkPipelineObject.h>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>

#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>
#include <cmtkStudyInfo.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class to encapsulate volume objects.
 * Depending on the particular geometry, volume data read from image files,
 * for example, can be either uniform or non-uniform. In order to provide
 * a single object to access data with varying structures, this class holds
 * a pointer to a volume object that can be replaced by another object of
 * potentially different structure. The instance of the wrapper class, 
 * however, remains the same so that client objects need not update their 
 * pointers.
 */
class VolumeWrapper : 
  /// Inherit from standard pipeline object.
  public PipelineObject 
{
public:
  /// Create new wrapper object.
  static VolumeWrapper* New();

  /// Return virtual class name.
  virtual const char *GetClassName() const { return "VolumeWrapper"; }

  /// Get pointer to the actual current volume.
  UniformVolume::SmartPtr& GetVolume() { return Volume; }
  const UniformVolume* GetVolume() const { return Volume; }
  
  /// Set new volume to encapsulate.
  void SetVolume( UniformVolume::SmartPtr& volume );

  /// Set new affine transformation.
  void SetAffineXform( AffineXform::SmartPtr& affineXform );

  /// Set new deformation.
  void SetWarpXform( WarpXform::SmartPtr& warpXform );

  /// Set new study info (DICOM tags).
  void SetStudyInfo( StudyInfo *const studyInfo );

  /// Get pointer to the affine transformation.
  AffineXform::SmartPtr& GetAffineXform() 
  { 
    return this->m_AffineXform; 
  }
  const AffineXform* GetAffineXform() const 
  { 
    return this->m_AffineXform; 
  }
  
  /// Get pointer to the deformation.
  WarpXform::SmartPtr& GetWarpXform() 
  { 
    return this->m_WarpXform; 
  }
  
  /// Get pointer to the study information (DICOM tags).
  StudyInfo* GetStudyInfo() 
  { 
    //    if ( m_StudyInfo ) m_StudyInfo->Reference(); 
    return m_StudyInfo; 
  }
  
protected:
  /// Default constructor.
  VolumeWrapper();

  /// Destructor.
  virtual ~VolumeWrapper();

private:
  /// Pointer to the encapsulated volume object.
  UniformVolume::SmartPtr Volume;

  /// The associated affine transformation.
  AffineXform::SmartPtr m_AffineXform;

  /// The associated deformation.
  WarpXform::SmartPtr m_WarpXform;

  /// DICOM tags from the original image data.
  StudyInfo *m_StudyInfo;
};

//@}

} // namespace cmtk

#endif
