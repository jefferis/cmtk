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
//  $Revision: 4052 $
//
//  $LastChangedDate: 2012-03-21 10:13:42 -0700 (Wed, 21 Mar 2012) $
//
//  $LastChangedBy: torsten_at_home $
//
*/

#ifndef __cmtkMagphanEMR051_h_included_
#define __cmtkMagphanEMR051_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Description of the Magphan EMR051 structural imaging phantom (a.k.a. ADNI Phantom).
class MagphanEMR051
{
public:
  /// This class.
  typedef MagphanEMR051 Self;

  /// Number of phantom spheres.
  static const unsigned int NumberOfSpheres = 165;

  /// Enumeration type for sphere colors.
  typedef enum
  {
    /// No color.
    SPHERE_COLOR_NONE,
    /// Green sphere.
    SPHERE_COLOR_GREEN,
    /// Yellow sphere.
    SPHERE_COLOR_YELLOW,
    /// Red sphere.
    SPHERE_COLOR_RED,
    /// Orange sphere.
    SPHERE_COLOR_ORANGE
  } SphereColorType;
  
  /// Data structure for phantom sphere location table.
  typedef struct __SphereEntryType
  {
    /// Sphere diameter in milimeters.
    Types::Coordinate m_Diameter;
    
    /// Sphere center location in "RAS" coordinate. Phantom center is the coordinate space origin.
    Types::Coordinate m_CenterLocation[3];
    
    /// Grams of Copper Sulfate Penta Hydrate per liter
    double m_GramsPerLiterCSPH;
    
    /** Estimated imaging T1.
     * For the 3.0cm spheres, this value is given in the phantom manual. For the remaining
     * spheres, we derive it from an exponential fit as T1=exp(-1.83364*CSPH+7.18582).
     *
     * We store only integer values because no fractional values are given in the phantom manual.
     */
    int m_EstimatedT1;
    
    /// Sphere color.
    SphereColorType m_Color;
  } SphereEntryType;
  
  /** Table of phantom sphere locations and properties.
   * Measurements were derived manually from the following document: http://www.phantomlab.com/library/pdf/magphan_adni_manual.pdf
   * They can, therefore, be used without reference to ADNI publications.
   */
  static const Self::SphereEntryType SphereTable[Self::NumberOfSpheres];

  /// Convenience access function - get sphere radius.
  static Types::Coordinate SphereRadius( const size_t i )
  {
    return 0.5 * Self::SphereTable[i].m_Diameter;
  }

  /// Convenience access function - get sphere center as 3D vector.
  static FixedVector<3,Types::Coordinate> SphereCenter( const size_t i )
  {
    return FixedVector<3,Types::Coordinate>( Self::SphereTable[i].m_CenterLocation );
  }

  /// Create a simulated T1-weighted image of the phantom spheres.
  static UniformVolume::SmartPtr GetPhantomImage( const Types::Coordinate resolution = 1.0 /*!< Pixel size of the output image; number of pixels is determined by this and the FOV needed to cover the entire phantom. */, 
						  const bool labels = false /*!< If this is set, each sphere is drawn with its intensity equal to the index in the marker table; otherwise, estimated T1 is used. */ );
  
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMagphanEMR051_h_included_
