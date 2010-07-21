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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctional.h"

#include "System/cmtkConsole.h"

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{
template<unsigned int NDegreeMul>
EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd )
{
  typedef EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr FunctionalPointer;
  FunctionalPointer functional;

  switch ( polynomialDegreeAdd )
    {
    case 0 :
      functional = FunctionalPointer( new EntropyMinimizationIntensityCorrectionFunctional<0,NDegreeMul> );
      break;
    case 1 :
      functional = FunctionalPointer( new EntropyMinimizationIntensityCorrectionFunctional<1,NDegreeMul> );
      break;
    case 2 :
      functional = FunctionalPointer( new EntropyMinimizationIntensityCorrectionFunctional<2,NDegreeMul> );
      break;
    case 3 :
      functional = FunctionalPointer( new EntropyMinimizationIntensityCorrectionFunctional<3,NDegreeMul> );
      break;
    case 4 :
      functional = FunctionalPointer( new EntropyMinimizationIntensityCorrectionFunctional<4,NDegreeMul> );
      break;
    default:
      StdErr.printf( "ERROR: combination of polynomial degrees %d (add) and %d (mul) not supported.\n", polynomialDegreeAdd, NDegreeMul );
      exit( 1 );
    }

  return functional;
}

EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul )
{
  typedef EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr FunctionalPointer;
  FunctionalPointer functional;

  switch ( polynomialDegreeMul )
    {
    case 0 :
      functional = CreateEntropyMinimizationIntensityCorrectionFunctional<0>( polynomialDegreeAdd );
      break;
    case 1 :
      functional = CreateEntropyMinimizationIntensityCorrectionFunctional<1>( polynomialDegreeAdd );
      break;
    case 2 :
      functional = CreateEntropyMinimizationIntensityCorrectionFunctional<2>( polynomialDegreeAdd );
      break;
    case 3 :
      functional = CreateEntropyMinimizationIntensityCorrectionFunctional<3>( polynomialDegreeAdd );
      break;
    case 4 :
      functional = CreateEntropyMinimizationIntensityCorrectionFunctional<4>( polynomialDegreeAdd );
      break;
    default:
      StdErr.printf( "ERROR: polynomial degree %d (mul) not supported.\n", polynomialDegreeMul );
      exit( 1 );
    }
  
  return functional;
}

EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr
CreateEntropyMinimizationIntensityCorrectionFunctional
( const unsigned int polynomialDegreeAdd, const unsigned int polynomialDegreeMul,
  EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr oldFunctional )
{
  EntropyMinimizationIntensityCorrectionFunctionalBase::SmartPtr newFunctional = CreateEntropyMinimizationIntensityCorrectionFunctional( polynomialDegreeAdd, polynomialDegreeMul );

  if ( oldFunctional )
    {
    CoordinateVector vOld;
    oldFunctional->GetParamVector( vOld );
    
    CoordinateVector vNew( newFunctional->ParamVectorDim() );
    vNew.SetAll( 0.0 );
    for ( size_t degreeAdd = 0; degreeAdd < oldFunctional->GetNumberOfMonomialsAdd(); ++degreeAdd )
      {
      vNew[degreeAdd] = vOld[degreeAdd];
      }
    for ( size_t degreeMul = 0; degreeMul < oldFunctional->GetNumberOfMonomialsMul(); ++degreeMul )
      {
      vNew[newFunctional->GetNumberOfMonomialsAdd() + degreeMul] = vOld[oldFunctional->GetNumberOfMonomialsAdd() + degreeMul];
      }
    }
  return newFunctional;
}

} // namespace cmtk

