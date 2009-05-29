/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkFunctionalAffine2D_h_included_
#define __cmtkFunctionalAffine2D_h_included_

#include <cmtkconfig.h>

#include <cmtkFunctional.h>

#include <cmtkMacros.h>
#include <cmtkSmartPtr.h>
#include <cmtkScalarImage.h>
#include <cmtkMatrix3x3.h>

#include <cmtkHistogram.h>
#include <cmtkJointHistogram.h>

#include <cmtkScalarImageSimilarity.h>
#include <cmtkRectangle.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for the affine registration two 2D images.
 */
class FunctionalAffine2D :
  /// This is a functional.
  public Functional 
{
  /// Number of degrees of freedom.
  igsGetSetMacro(int,NumberDOFs);

  /// Similarity metric.
  igsGetSetMacro(ScalarImageSimilarity::ID,SimilarityMeasure);

  /// Flag for histogram equalization of projection data.
  igsGetSetMacro(bool,HistogramEqualization);

public:
  /// This class.
  typedef FunctionalAffine2D Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Parent class.
  typedef Functional Superclass;

  /// Contructor.
  FunctionalAffine2D( ScalarImage::SmartPtr& refImage, ScalarImage::SmartPtr& fltImage, const IntROI2D* fltROI = NULL );

  /// Contructor.
  FunctionalAffine2D( std::vector<ScalarImage::SmartPtr>& refImage, std::vector<ScalarImage::SmartPtr>& fltImage, const IntROI2D* fltROI = NULL );

  /// Virtual destructor (dummy).
  ~FunctionalAffine2D() {};
  
  /// Set reference image.
  virtual void SetRefImage( ScalarImage::SmartPtr& refImage ) 
  {
    this->RefImages.clear();
    this->RefImages.push_back( refImage );
  }

  /// Set reference image.
  virtual void SetRefImages
  ( std::vector<ScalarImage::SmartPtr>& refImages ) 
  {
    this->RefImages = refImages;
  }

  /// Set floating image (and floating image ROI).
  virtual void SetFltImage( ScalarImage::SmartPtr& fltImage ) 
  {
    this->FltImages.clear();
    this->FltImagesROI.clear();
    this->FltImages.push_back( fltImage );
    this->FltImagesROI.push_back( fltImage );
  }

  /// Set floating image (and floating image ROI).
  virtual void SetFltImages
  ( std::vector<ScalarImage::SmartPtr>& fltImages ) 
  {
    this->FltImages = fltImages;
    this->FltImagesROI = fltImages;
  }

  /// Set minimum and maximum number of bins for histogram-based similarity.
  void SetNumberOfBins( const size_t minBins, const size_t maxBins = 0 );

  /// Compute functional value.
  virtual Self::ReturnType Evaluate();

  /// Compute functional value.
  virtual  Self::ReturnType EvaluateAt( CoordinateVector& v );

  /// Compute functional value and gradient.
  virtual  Self::ReturnType EvaluateWithGradient( CoordinateVector&, CoordinateVector&, const Types::Coordinate = 1 )
  { return 0; }


  /// Set parameter vector without evaluating the functional.
  virtual void SetParamVector( CoordinateVector& v );

  /** Return current parameter vector.
   * This is useful to retrieve the intial parameter values set by this class
   * upon initialization.
   */
  virtual void GetParamVector( CoordinateVector& v );

  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const { return 8; }

  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const { return this->NumberDOFs; }

  /// Return parameter stepping.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const;

private:
  /// Reference image.
  std::vector<ScalarImage::SmartPtr> RefImages;

  /// Floating image.
  std::vector<ScalarImage::SmartPtr> FltImages;

  /** ROI from floating image.
   * If the constructor is called without ROI, this is the same object as
   * FltImage.
   */
  std::vector<ScalarImage::SmartPtr> FltImagesROI;

  /// Persistent image similarity object.
  mutable ScalarImageSimilarityMemory ImageSimilarityMemory;
  
  /// Compute appropriate similarity measure for a single image pair.
  Self::ReturnType GetSimilarity( const ScalarImage* img0,  const ScalarImage* img1 ) const;

  /// Compute appropriate similarity measure for multiple images.
  Self::ReturnType GetSimilarity( std::vector<const ScalarImage*>& imgs0,  std::vector<const ScalarImage*>& imgs1 ) const;

  /// Parameter vector.
  CoordinateVector Parameters;

  /// Transformation.
  CoordinateMatrix3x3 Transformation;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFunctionalAffine2D_h_included_
