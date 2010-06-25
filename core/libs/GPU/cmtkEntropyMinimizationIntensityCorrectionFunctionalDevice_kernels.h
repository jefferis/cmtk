/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_functions_included_
#define __cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_functions_included_

/// Update output image using either additive or multiplicative bias field.
void cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( float* output, float* input, const int dims0, const int dims1, const int dims2, const int degree, const int multiply );

#endif // #ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_functions_included_
