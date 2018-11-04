/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010, 2013 SRI International
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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice.h"

#include <System/cmtkConsole.h>

namespace cmtk {

/** \addtogroup GPU */
//@{
template <unsigned int NDegreeMul>
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctionalDevice(
    const unsigned int polynomialDegreeAdd) {
  typedef EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
      FunctionalPointer;
  FunctionalPointer functional;

  switch (polynomialDegreeAdd) {
    case 0:
      functional = FunctionalPointer(
          new EntropyMinimizationIntensityCorrectionFunctionalDevice<
              0, NDegreeMul>);
      break;
    case 1:
      functional = FunctionalPointer(
          new EntropyMinimizationIntensityCorrectionFunctionalDevice<
              1, NDegreeMul>);
      break;
    case 2:
      functional = FunctionalPointer(
          new EntropyMinimizationIntensityCorrectionFunctionalDevice<
              2, NDegreeMul>);
      break;
    case 3:
      functional = FunctionalPointer(
          new EntropyMinimizationIntensityCorrectionFunctionalDevice<
              3, NDegreeMul>);
      break;
    case 4:
      functional = FunctionalPointer(
          new EntropyMinimizationIntensityCorrectionFunctionalDevice<
              4, NDegreeMul>);
      break;
    default:
      StdErr.printf(
          "ERROR: combination of polynomial degrees %u (add) and %u "
          "(mul) not supported.\n",
          polynomialDegreeAdd, NDegreeMul);
      exit(1);
  }

  return functional;
}

EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctionalDevice(
    const unsigned int polynomialDegreeAdd,
    const unsigned int polynomialDegreeMul) {
  typedef EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
      FunctionalPointer;
  FunctionalPointer functional;

  switch (polynomialDegreeMul) {
    case 0:
      functional =
          CreateEntropyMinimizationIntensityCorrectionFunctionalDevice<0>(
              polynomialDegreeAdd);
      break;
    case 1:
      functional =
          CreateEntropyMinimizationIntensityCorrectionFunctionalDevice<1>(
              polynomialDegreeAdd);
      break;
    case 2:
      functional =
          CreateEntropyMinimizationIntensityCorrectionFunctionalDevice<2>(
              polynomialDegreeAdd);
      break;
    case 3:
      functional =
          CreateEntropyMinimizationIntensityCorrectionFunctionalDevice<3>(
              polynomialDegreeAdd);
      break;
    case 4:
      functional =
          CreateEntropyMinimizationIntensityCorrectionFunctionalDevice<4>(
              polynomialDegreeAdd);
      break;
    default:
      StdErr.printf("ERROR: polynomial degree %u (mul) not supported.\n",
                    polynomialDegreeMul);
      exit(1);
  }

  return functional;
}

EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctionalDevice(
    const unsigned int polynomialDegreeAdd,
    const unsigned int polynomialDegreeMul,
    EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
        oldFunctional) {
  EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr newFunctional =
      CreateEntropyMinimizationIntensityCorrectionFunctionalDevice(
          polynomialDegreeAdd, polynomialDegreeMul);

  if (oldFunctional) {
    CoordinateVector vOld;
    oldFunctional->GetParamVector(vOld);

    CoordinateVector vNew(newFunctional->ParamVectorDim());
    vNew.SetAll(0.0);
    for (size_t degreeAdd = 0;
         degreeAdd < oldFunctional->GetNumberOfMonomialsAdd(); ++degreeAdd) {
      vNew[degreeAdd] = vOld[degreeAdd];
    }
    for (size_t degreeMul = 0;
         degreeMul < oldFunctional->GetNumberOfMonomialsMul(); ++degreeMul) {
      vNew[newFunctional->GetNumberOfMonomialsAdd() + degreeMul] =
          vOld[oldFunctional->GetNumberOfMonomialsAdd() + degreeMul];
    }
  }
  return newFunctional;
}

}  // namespace cmtk
