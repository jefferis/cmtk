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
  { 15, {   0.0, -2.9, 60.0 },  0.820, cmtk::Phantoms::SPHERE_COLOR_NONE }
  // 158x 1.0cm spheres
  //   Plane 0
  //   Plane 1
  //   Plane 1b
  //   Plane 2
  //   Plane 2b
  //   Plane 3
  //   Plane 3b
};

