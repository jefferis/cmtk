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

#include "cmtkMagphanEMR051.h"

/*
 * Measurements were derived manually from the following document: http://www.phantomlab.com/library/pdf/magphan_adni_manual.pdf
 * They can, therefore, be used without reference to ADNI publications.
 */
const cmtk::Phantoms::SphereEntryType cmtk::Phantoms::MagphanEMR051SphereTable[165] =
{
  // 1x 6.0cm sphere
  { 60, { 0.0, 0.0, 0.0 },      0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  // 4x 3.0cm spheres
  { 30, {  60.0,  15.0, 0 },    0.590, cmtk::Phantoms::SPHERE_COLOR_ORANGE }, // y coordinate estimated - not marked in construction drawing
  { 30, { -60.0,  15.0, 0 },    0.430, cmtk::Phantoms::SPHERE_COLOR_RED },    // y coordinate estimated - not marked in construction drawing
  { 30, {  60.0, -15.0, 0 },    0.220, cmtk::Phantoms::SPHERE_COLOR_GREEN },  // y coordinate estimated - not marked in construction drawing
  { 30, { -60.0, -15.0, 0 },    0.295, cmtk::Phantoms::SPHERE_COLOR_YELLOW }, // y coordinate estimated - not marked in construction drawing
  // 2x 1.5cm spheres
  { 15, {  89.0, -2.9,  0.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 15, {   0.0, -2.9, 60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  // 158x 1.0cm spheres
  //   Plane 0
  //     outer ring
  { 10, {  86.4, 0.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  86.4, 0.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -86.4, 0.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -86.4, 0.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  64.7, 0.0,  64.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  64.7, 0.0, -64.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -64.7, 0.0,  64.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -64.7, 0.0, -64.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, 0.0,  86.4 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, 0.0, -86.4 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, 0.0,  86.4 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, 0.0, -86.4 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {   0.0, 0.0,  91.5 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {   0.0, 0.0, -91.5 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     middle ring
  { 10, {  30.0, 0.0,  60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, 0.0, -60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, 0.0,  60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, 0.0, -60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  60.0, 0.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  60.0, 0.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -60.0, 0.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -60.0, 0.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     inner ring
  { 10, {  30.0, 0.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, 0.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, 0.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, 0.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     single inferior mid-sagittal sphere
  { 10, {   0.0, 0.0, -60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //   Plane 1
  //     outer ring
  { 10, {  86.5, -30.0,   0.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -86.5, -30.0,   0.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {   0.0, -30.0,  86.5 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {   0.0, -30.0, -86.5 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  61.0, -30.0,  61.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  61.0, -30.0, -61.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -61.0, -30.0,  61.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, -30.0,  81.1 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, -30.0, -81.1 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, -30.0,  81.1 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, -30.0, -81.1 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  81.1, -30.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  81.1, -30.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -81.1, -30.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -81.1, -30.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     middle ring
  { 10, {  60.0, -30.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  60.0, -30.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -60.0, -30.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -60.0, -30.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, -30.0,  60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, -30.0, -60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, -30.0,  60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, -30.0, -60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     inner ring
  { 10, {  30.0, -30.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  30.0, -30.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, -30.0,  30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -30.0, -30.0, -30.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {   0.0, -30.0,  40.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {   0.0, -30.0, -40.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     single superior mid-sagittal sphere
  { 10, {   0.0, -30.0,  60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //   Plane 1b
  //     +- 15mm lateral
  { 10, {  15.0,  29.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  29.1,  65.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  29.1,  85.2 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  29.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  29.1, -65.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  29.1, -85.2 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  29.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  29.1,  65.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  29.1,  85.2 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  29.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  29.1, -65.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  29.1, -85.2 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 40mm lateral
  { 10, {  40.0,  29.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  40.0,  29.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -40.0,  29.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -40.0,  29.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 45mm lateral
  { 10, {  45.0,  29.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  45.0,  29.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0,  29.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0,  29.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  45.0,  29.1,  65.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0,  29.1, -65.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  45.0,  29.1, -73.9 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0,  29.1,  73.9 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 64.5mm lateral
  { 10, {  64.5,  29.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -64.5,  29.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 73.9mm lateral
  { 10, {  73.9,  29.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -73.9,  29.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 85.2mm lateral
  { 10, {  85.2,  29.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  85.2,  29.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -85.2,  29.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -85.2,  29.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //   Plane 2
  //     +- 15mm lateral
  { 10, {  15.0, -60.0,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0, -60.0, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0, -60.0,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0, -60.0, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0, -60.0,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0, -60.0, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0, -60.0,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0, -60.0, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0, -60.0,  67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0, -60.0, -67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0, -60.0,  67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0, -60.0, -67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 45mm lateral
  { 10, {  45.0, -60.0,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  45.0, -60.0, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0, -60.0,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0, -60.0, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 48.7mm lateral
  { 10, {  48.7, -60.0,  48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  48.7, -60.0, -48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -48.7, -60.0,  48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -48.7, -60.0, -48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 67.3mm lateral
  { 10, {  67.3, -60.0,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  67.3, -60.0, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -67.3, -60.0,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -67.3, -60.0, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //   Plane 2b (same as Plane 2 but at y=+59.1)
  //     +- 15mm lateral
  { 10, {  15.0,  59.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  59.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  59.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  59.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  59.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  59.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  59.1,  45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  59.1, -45.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  59.1,  67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  15.0,  59.1, -67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  59.1,  67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -15.0,  59.1, -67.3 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 45mm lateral
  { 10, {  45.0,  59.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  45.0,  59.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0,  59.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -45.0,  59.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 48.7mm lateral
  { 10, {  48.7,  59.1,  48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  48.7,  59.1, -48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -48.7,  59.1,  48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -48.7,  59.1, -48.7 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //     +- 67.3mm lateral
  { 10, {  67.3,  59.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, {  67.3,  59.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -67.3,  59.1,  15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  { 10, { -67.3,  59.1, -15.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE },
  //   Plane 3
  //   Plane 3b
  {  0, {   0.0, 0.0, 0.0 },  0.0, cmtk::Phantoms::SPHERE_COLOR_NONE }
};

