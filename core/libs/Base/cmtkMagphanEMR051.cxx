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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include "cmtkMagphanEMR051.h"

/*
 * Measurements were derived manually from the following document: http://www.phantomlab.com/library/pdf/magphan_adni_manual.pdf
 * They can, therefore, be used without reference to ADNI publications.
 */
const cmtk::MagphanEMR051::SphereEntryType cmtk::MagphanEMR051::SphereTable[cmtk::MagphanEMR051::NumberOfSpheres] =
{
  // 
  // LICENSING EXCEPTION
  //   Unlike the remainder of this file, the table of phantom sphere coordinates
  //   is licensed under the CC BY 3.0 license (https://creativecommons.org/licenses/by/3.0/us/)
  //
  // 1x 6.0cm SNR sphere
  { "SNR",        60, { 0.0, 0.0, 0.0 },      0.820, 282, Self::SPHERE_COLOR_NONE }, //#000
  // 2x 1.5cm spheres
  { "15mm@90mm",  15, {  89.0, -2.9,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "15mm@60mm",  15, {   0.0, -2.9, -60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  // 4x 3.0cm spheres
  { "CNR-Orange", 30, {  60.0,  15.0, 0 },    0.590, 450, Self::SPHERE_COLOR_ORANGE }, // y coordinate estimated - not marked in construction drawing
  { "CNR-Red",    30, { -60.0,  15.0, 0 },    0.430, 600, Self::SPHERE_COLOR_RED },    // y coordinate estimated - not marked in construction drawing
  { "CNR-Yellow", 30, { -60.0, -15.0, 0 },    0.295, 750, Self::SPHERE_COLOR_YELLOW }, // y coordinate estimated - not marked in construction drawing
  { "CNR-Green",  30, {  60.0, -15.0, 0 },    0.220, 900, Self::SPHERE_COLOR_GREEN },  // y coordinate estimated - not marked in construction drawing
  // 158x 1.0cm spheres
  //   Plane 0
  //     outer ring
  { "10mm_0_01",  10, { -86.4, 0.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_02",  10, { -86.4, 0.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_03",  10, {  86.4, 0.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_04",  10, {  86.4, 0.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#010
  { "10mm_0_05",  10, { -64.7, 0.0,  64.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_06",  10, { -64.7, 0.0, -64.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_07",  10, {  64.7, 0.0,  64.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_08",  10, {  64.7, 0.0, -64.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_09",  10, { -30.0, 0.0,  86.4 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_10",  10, { -30.0, 0.0, -86.4 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_11",  10, {  30.0, 0.0,  86.4 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_12",  10, {  30.0, 0.0, -86.4 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_13",  10, {   0.0, 0.0,  91.5 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_14",  10, {   0.0, 0.0, -91.5 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#020
  //     middle ring
  { "10mm_0_15",  10, { -30.0, 0.0,  60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_16",  10, { -30.0, 0.0, -60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_17",  10, {  30.0, 0.0,  60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_18",  10, {  30.0, 0.0, -60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_19",  10, { -60.0, 0.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_20",  10, { -60.0, 0.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_21",  10, {  60.0, 0.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_22",  10, {  60.0, 0.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     inner ring
  { "10mm_0_23",  10, { -30.0, 0.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_24",  10, { -30.0, 0.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#030
  { "10mm_0_25",  10, {  30.0, 0.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_0_26",  10, {  30.0, 0.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     single inferior mid-sagittal sphere
  { "10mm_0_27",  10, {   0.0, 0.0,  60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     single right sphere
  { "10mm_0_28",  10, { -91.5, 0.0,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //   Plane 1
  //     mid-sagittal
  { "10mm_1_01",  10, {   0.0, -30.0,  40.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_02",  10, {   0.0, -30.0, -40.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_03",  10, {   0.0, -30.0,  86.5 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_04",  10, {   0.0, -30.0, -86.5 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_05",  10, {   0.0, -30.0,  60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_06",  10, {   0.0, -30.0, -60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#040
  //     +- 30mm lateral
  { "10mm_1_07",  10, {  30.0, -30.0,  81.1 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_08",  10, {  30.0, -30.0, -81.1 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_09",  10, { -30.0, -30.0,  81.1 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_10",  10, { -30.0, -30.0, -81.1 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_11",  10, {  30.0, -30.0,  60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_12",  10, {  30.0, -30.0, -60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_13",  10, { -30.0, -30.0,  60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_14",  10, { -30.0, -30.0, -60.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_15",  10, {  30.0, -30.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_16",  10, {  30.0, -30.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#050
  { "10mm_1_17",  10, { -30.0, -30.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_18",  10, { -30.0, -30.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 60mm lateral
  { "10mm_1_19",  10, {  60.0, -30.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_20",  10, {  60.0, -30.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_21",  10, { -60.0, -30.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_22",  10, { -60.0, -30.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 61mm lateral
  { "10mm_1_23",  10, {  61.0, -30.0,  61.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_24",  10, {  61.0, -30.0, -61.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_25",  10, { -61.0, -30.0,  61.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_26",  10, { -61.0, -30.0, -61.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#060
  //     +- 81.1mm lateral
  { "10mm_1_27",  10, {  81.1, -30.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_28",  10, {  81.1, -30.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_29",  10, { -81.1, -30.0,  30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_30",  10, { -81.1, -30.0, -30.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 86.5mm lateral
  { "10mm_1_31",  10, {  86.5, -30.0,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1_32",  10, { -86.5, -30.0,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //   Plane 1b
  //     +- 15mm lateral
  { "10mm_1b_01", 10, {  15.0,  29.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_02", 10, {  15.0,  29.1,  65.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_03", 10, {  15.0,  29.1,  85.2 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_04", 10, {  15.0,  29.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#070
  { "10mm_1b_05", 10, {  15.0,  29.1, -65.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_06", 10, {  15.0,  29.1, -85.2 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_07", 10, { -15.0,  29.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_08", 10, { -15.0,  29.1,  65.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_09", 10, { -15.0,  29.1,  85.2 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_10", 10, { -15.0,  29.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_11", 10, { -15.0,  29.1, -65.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_12", 10, { -15.0,  29.1, -85.2 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 40mm lateral
  { "10mm_1b_13", 10, {  40.0,  29.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_14", 10, {  40.0,  29.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#080
  { "10mm_1b_15", 10, { -40.0,  29.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_16", 10, { -40.0,  29.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 45mm lateral
  { "10mm_1b_17", 10, {  45.0,  29.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_18", 10, {  45.0,  29.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_19", 10, { -45.0,  29.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_20", 10, { -45.0,  29.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_21", 10, {  45.0,  29.1,  65.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_22", 10, { -45.0,  29.1, -65.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_23", 10, {  45.0,  29.1, -73.9 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_24", 10, { -45.0,  29.1,  73.9 },  0.820, 282, Self::SPHERE_COLOR_NONE }, //#090
  //     +- 64.5mm lateral
  { "10mm_1b_25", 10, {  64.5,  29.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_26", 10, { -64.5,  29.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 73.9mm lateral
  { "10mm_1b_27", 10, {  73.9,  29.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_28", 10, { -73.9,  29.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 85.2mm lateral
  { "10mm_1b_29", 10, {  85.2,  29.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_30", 10, {  85.2,  29.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_31", 10, { -85.2,  29.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_1b_32", 10, { -85.2,  29.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //   Plane 2
  //     +- 15mm lateral
  { "10mm_2_01",  10, {  15.0, -60.0,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_02",  10, {  15.0, -60.0, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_03",  10, { -15.0, -60.0,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_04",  10, { -15.0, -60.0, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_05",  10, {  15.0, -60.0,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_06",  10, {  15.0, -60.0, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_07",  10, { -15.0, -60.0,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_08",  10, { -15.0, -60.0, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_09",  10, {  15.0, -60.0,  67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_10",  10, {  15.0, -60.0, -67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_11",  10, { -15.0, -60.0,  67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_12",  10, { -15.0, -60.0, -67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 45mm lateral
  { "10mm_2_13",  10, {  45.0, -60.0,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_14",  10, {  45.0, -60.0, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_15",  10, { -45.0, -60.0,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_16",  10, { -45.0, -60.0, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 48.7mm lateral
  { "10mm_2_17",  10, {  48.7, -60.0,  48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_18",  10, {  48.7, -60.0, -48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_19",  10, { -48.7, -60.0,  48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_20",  10, { -48.7, -60.0, -48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 67.3mm lateral
  { "10mm_2_21",  10, {  67.3, -60.0,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_22",  10, {  67.3, -60.0, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_23",  10, { -67.3, -60.0,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2_24",  10, { -67.3, -60.0, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //   Plane 2b (same as Plane 2 but at y=+59.1)
  //     +- 15mm lateral
  { "10mm_2b_01", 10, {  15.0,  59.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_02", 10, {  15.0,  59.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_03", 10, { -15.0,  59.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_04", 10, { -15.0,  59.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_05", 10, {  15.0,  59.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_06", 10, {  15.0,  59.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_07", 10, { -15.0,  59.1,  45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_08", 10, { -15.0,  59.1, -45.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_09", 10, {  15.0,  59.1,  67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_10", 10, {  15.0,  59.1, -67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_11", 10, { -15.0,  59.1,  67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_12", 10, { -15.0,  59.1, -67.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 45mm lateral
  { "10mm_2b_13", 10, {  45.0,  59.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_14", 10, {  45.0,  59.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_15", 10, { -45.0,  59.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_16", 10, { -45.0,  59.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 48.7mm lateral
  { "10mm_2b_17", 10, {  48.7,  59.1,  48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_18", 10, {  48.7,  59.1, -48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_19", 10, { -48.7,  59.1,  48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_20", 10, { -48.7,  59.1, -48.7 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //     +- 67.3mm lateral
  { "10mm_2b_21", 10, {  67.3,  59.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_22", 10, {  67.3,  59.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_23", 10, { -67.3,  59.1,  15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_2b_24", 10, { -67.3,  59.1, -15.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //   Plane 3
  { "10mm_3_01",  10, {  28.3, -82.2,  28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_02",  10, {  28.3, -82.2, -28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_03",  10, { -28.3, -82.2,  28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_04",  10, { -28.3, -82.2, -28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_05",  10, {   0.0, -88.0,  25.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_06",  10, {   0.0, -88.0, -25.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_07",  10, {  25.0, -88.0,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_08",  10, { -25.0, -88.0,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3_09",  10, {   0.0, -89.5,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  //   Plane 3b
  { "10mm_3b_01", 10, {  28.3,  81.4,  28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_02", 10, {  28.3,  81.4, -28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_03", 10, { -28.3,  81.4,  28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_04", 10, { -28.3,  81.4, -28.3 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_05", 10, {   0.0,  87.2,  25.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_06", 10, {   0.0,  87.2, -25.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_07", 10, {  25.0,  87.2,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_08", 10, { -25.0,  87.2,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE },
  { "10mm_3b_09", 10, {   0.0,  88.6,   0.0 },  0.820, 282, Self::SPHERE_COLOR_NONE }
  // 
  // END LICENSING EXCEPTION
  //
};
