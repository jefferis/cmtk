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

#ifndef __cmtkMultiChannelRegistrationFunctional_h_included_
#define __cmtkMultiChannelRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkFunctional.h>

#include <cmtkUniformVolume.h>
#include <cmtkSmartPtr.h>

#include <cmtkUniformVolumeInterpolator.h>
#include <cmtkLinearInterpolator.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for multi-channel registration functionals. */
template<class TDataType = float, class TInterpolator = UniformVolumeInterpolator<Interpolators::Linear> >
class MultiChannelRegistrationFunctional :
  /** Inherit functional interface. */
  public Functional
{
public:
  /** This class. */
  typedef MultiChannelRegistrationFunctional<TDataType,TInterpolator> Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** Superclass. */
  typedef Functional Superclass;

  /** Real value type for data representation. */
  typedef TDataType DataType;

  /** Default constructor. */
  MultiChannelRegistrationFunctional() : m_NumberOfChannels( 0 ), m_NormalizedMI( false ) {}
  
  /** Destructor: free all converted image arrays. */
  virtual ~MultiChannelRegistrationFunctional()
  {
    this->ClearAllChannels();
  }

  /** Set flag for normalized vs. standard MI */
  void SetNormalizedMI( const bool nmi = true )
  {
    this->m_NormalizedMI = nmi;
  }

  /** Reset channels, clear all images. */
  virtual void ClearAllChannels();

  /** Add reference channel. */
  virtual void AddReferenceChannel( UniformVolume::SmartPtr& channel );

  /** Add reference channels from stl container. */
  template <class ForwardIterator>
  void AddReferenceChannels( ForwardIterator first, ForwardIterator last )
  {
    while ( first != last )
      {
      this->AddReferenceChannel( *first );
      ++first;
      }
  }

  /** Get number of reference channels. */
  size_t GetNumberOfReferenceChannels() const { return this->m_ReferenceChannels.size(); }

  /** Get a reference channel image. */
  UniformVolume::SmartPtr& GetReferenceChannel( const size_t idx ) { return this->m_ReferenceChannels[idx]; }

  /** Get constant pointer to reference channel image. */
  const UniformVolume* GetReferenceChannel( const size_t idx ) const { return this->m_ReferenceChannels[idx]; }

  /** Add floating channel. */
  virtual void AddFloatingChannel( UniformVolume::SmartPtr& channel );

  /** Add floating channels from stl container. */
  template <class ForwardIterator>
  void AddFloatingChannels( ForwardIterator first, ForwardIterator last )
  {
    while ( first != last )
      {
      this->AddFloatingChannel( *first );
      ++first;
      }
  }

  /** Get number of floating channels. */
  size_t GetNumberOfFloatingChannels() const { return this->m_FloatingChannels.size(); }

  /** Get a floating channel image. */
  UniformVolume::SmartPtr& GetFloatingChannel( const size_t idx ) { return this->m_FloatingChannels[idx]; }

  /** Get constant pointer to floating channel image. */
  const UniformVolume* GetFloatingChannel( const size_t idx ) const { return this->m_FloatingChannels[idx]; }

  /** Vector of reference images. */
  std::vector<UniformVolume::SmartPtr> m_ReferenceChannels;

  /** Vector of floating images. */
  std::vector<UniformVolume::SmartPtr> m_FloatingChannels;

protected:
  /// Interpolators for the floating image channels.
  std::vector<typename TInterpolator::SmartPtr> m_FloatingInterpolators;

  /** Total number of channels. This is the sum of the floating and reference channel vector sizes. */
  size_t m_NumberOfChannels;

  /// Grid dimensions of the reference volume.
  int m_ReferenceDims[3];

  /// Extents of the reference volume in real-world coordinates.
  Types::Coordinate m_ReferenceSize[3];

  /// Inverse pixel deltas of the reference volume.
  Types::Coordinate m_ReferenceInvDelta[3];

  /// Beginning of the rectangular crop region in the reference volume.
  int m_ReferenceCropFrom[3];

  /// End of the rectangular crop region in the reference volume.
  int m_ReferenceCropTo[3];

  /// Grid dimensions of the floating volume.
  int m_FloatingDims[3];

  /// Extents of the floating volume in real-world coordinates.
  Types::Coordinate m_FloatingSize[3];

  /// Inverse pixel deltas of the floating volume.
  Types::Coordinate m_FloatingInverseDelta[3];

  /// Starting coordinate of the floating's cropping region.
  Types::Coordinate m_FloatingCropFrom[3];

  /// End coordinate of the floating's cropping region.
  Types::Coordinate m_FloatingCropTo[3];
 
  /** Update all transformation-related data after init or refine. */
  virtual void NewReferenceChannelGeometry() {}

protected:
  /** Flag for normalized vs. standard mutual information. */
  bool m_NormalizedMI;

private:
  /** Verify size and geometry of newly added image channels against already added channels. */
  void VerifyImageSize( const UniformVolume* imgA, const UniformVolume* imgB );
};

//@}

} // namespace cmtk

#include <cmtkMultiChannelRegistrationFunctional.txx>

#endif // #ifndef __cmtkMultiChannelRegistrationFunctional_h_included_
