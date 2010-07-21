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

#ifndef __cmtkGroupwiseRegistrationFunctionalXformTemplate_h_included_
#define __cmtkGroupwiseRegistrationFunctionalXformTemplate_h_included_

#include <cmtkconfig.h>

#include "Registration/cmtkGroupwiseRegistrationFunctionalXformTemplateBase.h"

#include "System/cmtkSmartPtr.h"
#include "System/cmtkThreads.h"
#include "System/cmtkThreadPool.h"

#include "Base/cmtkUniformVolume.h"
#include "Base/cmtkXform.h"

#include <vector>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Trannsformation-dependent class template for groupwise registration functionals.
 * This template provides the common generic interface for all transformation-model dependent 
 * specialized templates.
 */
template<class TXform>
class GroupwiseRegistrationFunctionalXformTemplate : 
  /** Inherit from generic groupwise functional. */
  public GroupwiseRegistrationFunctionalXformTemplateBase<TXform>
{
public:
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalXformTemplateBase<TXform> Superclass;
  
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalXformTemplate<TXform> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  GroupwiseRegistrationFunctionalXformTemplate();

  /// Destructor.
  virtual ~GroupwiseRegistrationFunctionalXformTemplate() {}

protected:
  /** Interpolate given moving image to template.
   *\param idx Index of of to reformat to template. This also determines which
   *  transformation is used.
   *\param destination The reformatted pixel data is stored in this array.
   *  Sufficient memory (for as many pixels as there are in the template grid)
   *  must be allocated there.
   */
  virtual void InterpolateImage( const size_t idx, byte* const destination );

private:
  /// Thread parameters with no further data.
  typedef ThreadParameters<Self> ThreadParametersType;

  /// Thread function parameters for image interpolation.
  class InterpolateImageThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParametersType
  {
  public:
    /// Index of the image to be interpolated.
    size_t m_Idx;

    /// Pointer to storage that will hold the reformatted pixel data.
    byte* m_Destination;    

    /// Return parameter: number of reformatted pixels outside floating image field of view. This is to detect pathological transformation parameters.
    size_t m_NumberOfOutsidePixels;
  };
  
  /// Task info blocks.
  std::vector<InterpolateImageThreadParameters> m_InterpolateTaskInfo;

  /// Image interpolation thread function.
  static void InterpolateImageThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCont );

  /// Image interpolation thread function with probabilistic sampling.
  static void InterpolateImageProbabilisticThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCont );
};

//@}

} // namespace cmtk

#include "Registration/cmtkGroupwiseRegistrationFunctionalXformTemplate_Affine.h"
#include "Registration/cmtkGroupwiseRegistrationFunctionalXformTemplate_SplineWarp.h"

#endif // #ifndef __cmtkGroupwiseRegistrationFunctionalXformTemplate_h_included_
