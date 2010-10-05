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

#ifndef __cmtkActiveShapeModel_h_included_
#define __cmtkActiveShapeModel_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkDirectionSet.h>

#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class for a three-dimensional active shape model.
class ActiveShapeModel
{
public:
  /// Smart pointer to active shape model.
  typedef SmartPointer<ActiveShapeModel> SmartPtr;

  /// Number of points in this model.
  unsigned int NumberOfPoints;

  /// Get number of points.
  unsigned int GetNumberOfPoints() const { return NumberOfPoints; }

  /// Point positions of the mean shape.
  CoordinateVector::SmartPtr Mean;

  /// Number of modes of variation in this model.
  unsigned int NumberOfModes;

  /// Get number of modes.
  unsigned int GetNumberOfModes() const { return NumberOfModes; }

  /// Delta vectors for the modes of variation.
  DirectionSet::SmartPtr Modes;

  /// Eigenvalue (= variance) associated with each mode
  CoordinateVector::SmartPtr ModeVariances;

  /// Default constructor.
  ActiveShapeModel() : Mean( NULL ), Modes( NULL ) {}

  /// Construct using given mean and modes.
  ActiveShapeModel( CoordinateVector::SmartPtr& mean, DirectionSet::SmartPtr& modes, CoordinateVector::SmartPtr& modeVariances );

  /** Generative constructor.
   * For a description of the parameters, see the Construct() member function.
   *\see ActiveShapeModel::Construct
   */
  ActiveShapeModel
  ( const Types::Coordinate *const* trainingSet, 
    const unsigned int numberOfSamples,
    const unsigned int numberOfPoints, 
    const unsigned int numberOfModes ) :
    Mean( NULL ), 
    Modes( NULL )
  {
    this->Construct( trainingSet, numberOfSamples, numberOfPoints, numberOfModes ); 
  }

  /** Construct model.
   *\param trainingSet The training set. This is a an array with size
   * [numberOfSamples]. Each entry is a pointer to an array of size
   * [numberOfPoints] Types::Coordinate values. Each of these latter arrays is a
   * sample vector from the training set, for example a sequential 
   * concatenation of 2-D or 3-D point coordinates. The order of these values
   * within the vectors is irrelevant, as long as the order is the same in all
   * vectors. This order also defines the order of parameters in the generated
   * model.
   *\param numberOfSamples The number of samples in the training set. This is
   * the size of the pointer array given as parameter "trainingSet".
   *\param numberOfPoints This is the number of values in each array pointed
   * to by a pointer in the "trainingSet" array.
   *\param numberOfModes Number of modes in the generated model. This can be
   * at most as many as "numberOfSamples", the number of samples in the 
   * training set. If the value of this parameter is to large, it will 
   * automatically be reduced to equal numberOfSamples.
   *\return Percentage of the total variance of the training set that is 
   * explained by the generated model.
   */
  float Construct( const Types::Coordinate *const* trainingSet, const unsigned int numberOfSamples, const unsigned int numberOfPoints, const unsigned int numberOfModes );

  /** Generate a model instance.
   */
  Types::Coordinate* Generate( Types::Coordinate *const instance, const Types::Coordinate* modeWeights ) const;

  /** Decompose a vector into mean and modes of this model.
   *@param input Input vector.
   *@param weights Weights of the modes that make up the given input vector.
   * This parameter is optional. If not given, no weights will be returned.
   *@return The value of the multivariate Gaussian PDF represented by this 
   * model atr the location of the input vector.
   */
  float Decompose( const CoordinateVector* input, Types::Coordinate *const weights = NULL ) const;

 protected:
  /** Allocate data structures.
   * If this instance already has data structures allocated, these are
   * deallocated first.
   */
  void Allocate( const unsigned int numberOfPoints, const unsigned int numberOfModes );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkActiveShapeModel_h_included_
