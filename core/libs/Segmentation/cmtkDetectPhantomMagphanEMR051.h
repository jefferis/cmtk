/*
//
//  Copyright 2012, 2013 SRI International
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

#ifndef __cmtkDetectPhantomMagphanEMR051_h_included_
#define __cmtkDetectPhantomMagphanEMR051_h_included_

#include <cmtkconfig.h>

#include <System/cmtkException.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkDetectedPhantomMagphanEMR051.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Class for detecting landmark locations of the Magphan EMR051 structural imaging phantom
 * (a.k.a The ADNI Phantom).
 */
class DetectPhantomMagphanEMR051
{
public:
  /// This class.
  typedef DetectPhantomMagphanEMR051 Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const for  this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Spatial coordinate vector.
  typedef UniformVolume::SpaceVectorType SpaceVectorType;

  /// Exception thrown if a landmark sphere cannot be localized in the search region.
  class NoSphereInSearchRegion : public Exception {};
  
  /// Exception thrown if the field of view is insufficient.
  class OutsideFieldOfView : public Exception 
  {
  public:
    // Constructor takes index and predicted location of offending landmark.
    OutsideFieldOfView( const size_t idx, const UniformVolume::CoordinateVectorType& v ) : m_Idx( idx ), m_Location( v ) {}

    // Offending landmark index.
    size_t m_Idx;

    // Offending predicted location.
    UniformVolume::CoordinateVectorType m_Location;
  };
  
  /// Constructor: detect all landmark spheres.
  DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage, const bool tolerant = false /*!< If this is set, partial sphere truncation will be tolerated.*/ );

  /// Get comprehensive description of phantom as detected in image.
  DetectedPhantomMagphanEMR051::SmartPtr GetDetectedPhantom();
  
  /// Get rigid phantom-to-image transformation.
  AffineXform::SmartConstPtr GetPhantomToImageTransformationRigid() const
  {
    return this->m_PhantomToImageTransformationRigid;
  }
  
  /// Get affine phantom-to-image transformation.
  AffineXform::SmartConstPtr GetPhantomToImageTransformationAffine() const
  {
    return this->m_PhantomToImageTransformationAffine;
  }
  
  /// Get the image-space label map of detected spheres.
  UniformVolume::SmartPtr GetDetectedSpheresLabelMap();

  /// Get expected landmark locations given rigid phantom-to-image transformation.
  LandmarkList GetExpectedLandmarks( const bool includeUnreliable = false /*!< If true, include unreliable landmarks, i.e., SNR and CNR spheres. */ ) const;

  /// Get actual, detected landmark locations.
  LandmarkList GetDetectedLandmarks( const bool includeOutliers = false /*!< If true, include landmarks detected as outliers based on linear affine transformation fitting residual */ ) const;

private:
  /// Flag for correction of (linear) bias field for each sphere.
  bool m_CorrectSphereBiasField;

  /** Flag for tolerant operation - if set, we are lenient when spheres are slightly truncated.
   * This should be considered a last resort, and both phantom scans and results should be carefully inspected.
   */
  bool m_TolerateTruncation;

  /// Erode SNR sphere by this many pixels for SNR computation
  byte m_ErodePixelsSNR;
  
  /// Erode CNR spheres by this many pixels for SNR computation
  byte m_ErodePixelsCNR;
  
  /// Image of the phantom.
  UniformVolume::SmartConstPtr m_PhantomImage;

  /** Evolving exclusion mask.
   * When a sphere is detected, its volume is marked as off-limits herein so other spheres are not incorrectly detected in the same place.
   */
  UniformVolume::SmartPtr m_ExcludeMask;

  /** Temporary inclusion mask.
   * When we detect a sphere in a specific pre-determined area, this mask contains as non-zero the candidate pixels.
   */
  UniformVolume::SmartPtr m_IncludeMask;

  /// Class for landmarks and validity flags
  class LandmarkType
  {
  public:
    /// Default constructor.
    LandmarkType() : m_Location( 0 ), m_Valid( false ) {}

    /// Constructor.
    LandmarkType( const Self::SpaceVectorType& location, const bool valid = true ) : m_Location( location ), m_Valid( valid ) {}

    /// Location of the landmark
    Self::SpaceVectorType m_Location;

    /// Is this landmark valid?
    bool m_Valid;
  };

  /// The detected sphere centroid landmarks in image space.
  std::vector<Self::LandmarkType> m_Landmarks;

  /// Fitted rigid transformation from phantom to image coordinates.
  AffineXform::SmartPtr m_PhantomToImageTransformationRigid;

  /// Fitted affine transformation from phantom to image coordinates.
  AffineXform::SmartPtr m_PhantomToImageTransformationAffine;

  /// Residuals of landmark locations after linear transformation fit.
  std::vector<Types::Coordinate> m_LandmarkFitResiduals;

  /** Find a sphere of given radius.
   *\return Location of the sphere center - this is the point where the filter response is maximal.
   */
  Self::SpaceVectorType FindSphere( const TypedArray& filterResponse /*!< Response of the matched filter used for sphere finding. */ );

  /** Find a sphere in a band of given radius.
   * If the given search region is already excluded from searching based on previously identified spheres, then the
   * NoSphereInSearchRegion exception is thrown.
   *\return Location of the sphere center. This is the location of maximum filter response in the search band, minus exclusion.
   */
  Self::SpaceVectorType FindSphereAtDistance( const TypedArray& filterResponse /*!< Response of the matched filter used for sphere finding. */, 
					      const Self::SpaceVectorType& bandCenter /*!< Center of the band to search in. */, 
					      const Types::Coordinate bandRadius /*!< Radius of the band to search in. If this is zero, the search region is a sphere around bandCenter.*/, 
					      const Types::Coordinate bandWidth /*!< Width of the band to search in.*/ );
  
  /// Refine sphere position based on intensity-weighted center of mass.
  Self::SpaceVectorType RefineSphereLocation( const Self::SpaceVectorType& estimate, const Types::Coordinate radius, const int label );

  /// Get margin (in pixels) for the bipolar sphere detection matched filter.
  int GetBipolarFilterMargin() const
  {
    return 1;
  }

  /// Get safety margin (in mm) around detected spheres - no other sphere centers are permitted within this margin.
  Types::Coordinate GetSphereExcludeSafetyMargin() const
  {
    return 10.0;
  }

  /// Get margin in pixels for center-of-mass-based refinement.
  int GetRefineMarginPixels() const
  {
    return 1;
  }

  /// Get threshold for detecting outliers based on landmark fitting residuals.
  Types::Coordinate GetLandmarkFitResidualThreshold() const
  {
    return 5.0; // if we're off bymore than 5mm, that's probably not due to actual distortion
  }

  /** Compute landmark fitting residuals under given linear transformation
   *\return The maximum residual over all landmarks. This can be compared with a threshold to
   * determine whether refinement of outliers is necessary.
   */
  Types::Coordinate ComputeLandmarkFitResiduals( const AffineXform& xform /*!< Linear transformation fitted to the landmarks.*/ );

  /// Try to refine outlier (by current fitted linear transformation) landmarks.
  void RefineOutlierLandmarks( const TypedArray& filterResponse /*!< Existing filter response map. */ );

  /// Exclude outlier landmarks and re-fit linear transformation.
  void ExcludeOutlierLandmarks();

  /// Get the mean and standard deviation of intensities within a sphere of given location and radius.
  void GetSphereMeanStdDeviation( Types::DataItem& mean /*!< Reference to return mean intensity */, Types::DataItem& stdev /*!< Reference to return intensity standard deviation */, 
				  const Self::SpaceVectorType& center /*!< Center coordinate of the sphere. */, const Types::Coordinate radius /*!< Radius of the sphere */, 
				  const int erodeBy /*!< Number of pixels to erode the sphere by before computing mean and standard deviation. */,
				  const int biasFieldDegree /*!< Polynomial degree of the estimated bias field before computing mean and standard deviation (0 = no bias field correction) */ );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDetectPhantomMagphanEMR051_h_included_
