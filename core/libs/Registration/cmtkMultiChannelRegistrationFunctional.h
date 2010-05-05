/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkMultiChannelRegistrationFunctional_h_included_
#define __cmtkMultiChannelRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkMultiChannelRegistrationFunctionalBase.h>

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
template<class TInterpolator = UniformVolumeInterpolator<Interpolators::Linear> >
class MultiChannelRegistrationFunctional :
  /** Inherit non-template implementation and interface. */
  public MultiChannelRegistrationFunctionalBase
{
public:
  /** This class. */
  typedef MultiChannelRegistrationFunctional<TInterpolator> Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** Superclass. */
  typedef MultiChannelRegistrationFunctionalBase Superclass;

  /** Add floating channel. */
  virtual void AddFloatingChannel( UniformVolume::SmartPtr& channel );

protected:
  /// Interpolators for the floating image channels.
  std::vector<typename TInterpolator::SmartPtr> m_FloatingInterpolators;
};

//@}

} // namespace cmtk

#include <cmtkMultiChannelRegistrationFunctional.txx>

#endif // #ifndef __cmtkMultiChannelRegistrationFunctional_h_included_
