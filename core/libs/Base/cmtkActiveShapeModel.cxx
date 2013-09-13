/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include <Base/cmtkActiveShapeModel.h>

#include <Base/cmtkSymmetricMatrix.h>
#include <Base/cmtkEigenSystemSymmetricMatrix.h>

#include <System/cmtkConsole.h>

#include <math.h>
#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

ActiveShapeModel::ActiveShapeModel
( CoordinateVector::SmartPtr& mean, DirectionSet::SmartPtr& modes, CoordinateVector::SmartPtr& modeVariances ) :
  Mean( mean ), Modes( modes ), ModeVariances( modeVariances )
{
  NumberOfPoints = Mean->Dim;
  NumberOfModes = Modes->GetNumberOfDirections();
}

float 
ActiveShapeModel::Construct
( const Types::Coordinate *const* trainingSet, const unsigned int numberOfSamples,
  const unsigned int numberOfPoints, const unsigned int numberOfModes )
{
  if ( numberOfSamples < numberOfModes ) 
    {
    StdErr << "WARNING: number of modes of an ASM can be no higher than number of training samples.\n";
    this->Allocate( numberOfPoints, numberOfSamples );
    } 
  else
    {
    this->Allocate( numberOfPoints, numberOfModes );
    }
  

  // first, compute mean shape
  Types::Coordinate *meanPtr = Mean->Elements;
  for ( unsigned int point = 0; point < NumberOfPoints; ++point, ++meanPtr ) 
    {
    Types::Coordinate mean = trainingSet[0][point];
    for ( unsigned int sample = 1; sample < numberOfSamples; ++sample ) 
      {
      mean += trainingSet[sample][point];
      }
    (*meanPtr) = (mean / numberOfSamples );
    }
  
  // now generate covariance matrix; actually, we're using a slightly
  // modified approach following Cootes' 1995 CVIU paper. This is much
  // more efficient when the number of samples is smaller than the
  // number of data dimensions.
  SymmetricMatrix<Types::Coordinate> cc( numberOfSamples );
  
  for ( unsigned int sampleY = 0; sampleY < numberOfSamples; ++sampleY )
    {
    for ( unsigned int sampleX = 0; sampleX <= sampleY; ++sampleX ) 
      {
      Types::Coordinate ccXY = 0;
      
      const Types::Coordinate* meanPtr2 = Mean->Elements;
      for ( unsigned int point = 0; point < NumberOfPoints; ++point, ++meanPtr2 )
	{
	ccXY += ( trainingSet[sampleX][point] - (*meanPtr2) ) * ( trainingSet[sampleY][point] - (*meanPtr2) );
	}
      cc(sampleX,sampleY) = ccXY / numberOfSamples;
      }
    }
  
  // here comes the hard part: compute Eigenvectors of cc...
  // we do this in a separate routine, for clarity.
  const EigenSystemSymmetricMatrix<Types::Coordinate> eigensystem( cc );

  // determine permutation that orders eigenvectors by descending eigenvalues
  const std::vector<Types::Coordinate> eigenvalues = eigensystem.GetEigenvalues();
  std::vector<unsigned int> permutation( numberOfSamples );
  // initialize permutation array
  for ( unsigned int i = 0; i < numberOfSamples; ++i )
    permutation[i] = i;
  
  // now do a simple bubble sort
  bool sorted = false;
  while ( ! sorted ) 
    {
    sorted = true;
    for ( unsigned int i = 0; i < numberOfSamples-1; ++i )
      if ( eigenvalues[permutation[i]] < eigenvalues[permutation[i+1]] ) 
	{
	std::swap( permutation[i], permutation[i+1] );
	sorted = false;
	}
    }
  
  // now, we need to convert the eigenvectors of the simplified matrix
  // back to those of the actual covariance matrix. Again, this follows
  // Cootes et al., CVIU 1995
  for ( unsigned int mode = 0; mode < NumberOfModes; ++mode ) 
    {    
    ModeVariances->Elements[mode] = eigenvalues[permutation[mode]];

    Types::Coordinate* modePtr = (*Modes)[mode]->Elements;
    for ( unsigned int point = 0; point < NumberOfPoints; ++point, ++modePtr ) 
      {
      unsigned int fromMode = permutation[mode];
      Types::Coordinate meanValue = Mean->Elements[point];
      
      *modePtr = 0;
      for ( unsigned int sample = 0; sample < numberOfSamples; ++sample )
	*modePtr += (eigensystem.EigenvectorElement(sample,fromMode) * (trainingSet[sample][point] - meanValue) );
      }
    
    // finally, normalize mode vectors... if Geremy is right ;)
    (*(*Modes)[mode]) *= (sqrt( eigenvalues[permutation[mode]] ) / (*Modes)[mode]->EuclidNorm());
    }
  
  return 0;
}

Types::Coordinate*
ActiveShapeModel::Generate
( Types::Coordinate *const instance, const Types::Coordinate* modeWeights ) const
{
  Types::Coordinate* target = instance;
  if ( !target )
    target = Memory::ArrayC::Allocate<Types::Coordinate>( NumberOfPoints );

  memcpy( target, Mean->Elements, sizeof( *target ) * NumberOfPoints );

  if ( modeWeights )
    {
    for ( unsigned int mode = 0; mode < NumberOfModes; ++mode ) 
      {
      Types::Coordinate modeWeight = modeWeights[mode];
      
      Types::Coordinate* targetPtr = target;
      const Types::Coordinate* modePtr = (*Modes)[mode]->Elements;
      
      for ( unsigned int point = 0; point < NumberOfPoints; ++point, ++targetPtr, ++modePtr )
	(*targetPtr) += ( modeWeight * (*modePtr) );
      }
    }
  
  return target;
}

float
ActiveShapeModel::Decompose
( const CoordinateVector* input, Types::Coordinate *const weights ) const
{
  std::vector<Types::Coordinate> w( this->NumberOfModes );
  CoordinateVector deviation( *input );
  deviation -= *(this->Mean);

#define RETURN_PDF
#ifdef RETURN_PDF
  float pdf = 1.0;
  for ( size_t mode = 0; mode < this->NumberOfModes; ++mode ) 
    {
    const CoordinateVector* thisMode = (*this->Modes)[mode];
    // since Modes are orthogonal basis, we can decompose using scalar product
    w[mode] = (deviation * *thisMode) / thisMode->EuclidNorm();
    
    const Types::Coordinate variance = (*(this->ModeVariances))[mode];
    pdf *= static_cast<float>( exp( -(w[mode]*w[mode]) / (2.0 * variance) ) / sqrt( 2.0 * M_PI * variance) );
    }
#else
  float distance = 0.0;
  for ( size_t mode = 0; mode < this->NumberOfModes; ++mode ) 
    {
    const CoordinateVector* thisMode = (*this->Modes)[mode];
    // since Modes are orthogonal basis, we can decompose using scalar product
    w[mode] = (deviation * *thisMode) / thisMode->EuclidNorm();
    
    const Types::Coordinate variance = (*(this->ModeVariances))[mode];
    distance += w[mode] * w[mode] / variance;
    }
  distance = sqrt( distance );
#endif
  
  if ( weights )
    memcpy( weights, &w[0], this->NumberOfModes * sizeof( *weights ) );
  
#ifdef RETURN_PDF
  return pdf;
#else
  return distance;
#endif
}

void
ActiveShapeModel::Allocate
( const unsigned int numberOfPoints, const unsigned int numberOfModes )
{
  NumberOfModes = numberOfModes;
  NumberOfPoints = numberOfPoints;

  Modes = DirectionSet::SmartPtr( new DirectionSet( NumberOfPoints ) );
  for ( unsigned int mode = 0; mode < NumberOfModes; ++mode )
    Modes->push_back( CoordinateVector::SmartPtr( new CoordinateVector( NumberOfPoints ) ) );
  ModeVariances = CoordinateVector::SmartPtr( new CoordinateVector( NumberOfModes ) );
  Mean = CoordinateVector::SmartPtr( new CoordinateVector( NumberOfPoints ) );
}

} // namespace cmtk
