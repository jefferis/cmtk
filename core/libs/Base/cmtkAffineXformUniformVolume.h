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

#ifndef __cmtkAffineXformUniformVolume_h_included_
#define __cmtkAffineXformUniformVolume_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkXformUniformVolume.h"
#include "Base/cmtkUniformVolume.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Pre-compute transformation for grid locations in a uniform volume.
 */
class AffineXformUniformVolume :
  /// Inherit from class to prevent copying.
  public XformUniformVolume
{
public:
  /// This class.
  typedef AffineXformUniformVolume Self;

  /// Parent class.
  typedef XformUniformVolume Superclass;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Constructor.
  AffineXformUniformVolume( const UniformVolume& volume, const AffineXform& xform );
  
  /// Virtual destructor.
  virtual ~AffineXformUniformVolume() {}
  
  /** Get transformed location of linked grid pixel.
   */
  virtual void GetTransformedGrid( Vector3D& v, const int idxX, const int idxY, const int idxZ ) const
  {
    ( (v = this->m_VolumeAxesX[idxX]) 
      += this->m_VolumeAxesY[idxY]) 
      += this->m_VolumeAxesZ[idxZ]; 
  }
  
  /** Get transformed locations of a series (scanline) of linked grid pixels.
   */
  virtual void GetTransformedGridSequence( Vector3D *const v, const size_t numPoints, const int idxX, const int idxY, const int idxZ ) const
  {
    Vector3D v0 = this->m_VolumeAxesY[idxY];
    v0 += this->m_VolumeAxesZ[idxZ];
    
    int idx = idxX;
    for ( size_t n=0; n < numPoints; ++n, ++idx ) 
      {
      (v[n] = v0) += this->m_VolumeAxesX[idx];
      }
  }
  
private:
  /// Axes hash for the points of a registered Volume.
  std::vector<Vector3D> m_VolumeAxesX;

  std::vector<Vector3D> m_VolumeAxesY;

  std::vector<Vector3D> m_VolumeAxesZ;
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkAffineXformUniformVolume_h_included_
