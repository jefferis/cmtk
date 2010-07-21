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

#ifndef __cmtkMultiLevelOptimizer_h_included_
#define __cmtkMultiLevelOptimizer_h_included_

#include <cmtkconfig.h>

#include "Registration/cmtkOptimizer.h"
#include "Registration/cmtkOptimizerBase.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Multi-level optimizer.
 * This class implements a multi-level optimizer. In fact, this is a
 * meta-optimizer, which uses a separate instance of an actual low-level
 * optimizer to perform the actual optimization. This class keeps a list of
 * several different functionals with corresponding initial and final
 * optimization step sizes that are sequentially handled by the low-level
 * optimizer.
 */
class MultiLevelOptimizer : 
    public OptimizerBase
{
public:
  /** Constructor.
   */
  MultiLevelOptimizer( Optimizer::SmartPtr& optimizer )
    : m_Optimizer( optimizer ) {}

  /// Set low-level optimizer.
  void SetOptimizer( Optimizer::SmartPtr& optimizer )
  {
    this->m_Optimizer = optimizer;
  }

  /// Virtual destructor.
  virtual ~MultiLevelOptimizer() {}

  /** Add functional and step sizes.
   *\return Number of functionals defined so far.
   */
  size_t AddFunctional( Functional::SmartPtr functional, const Types::Coordinate initialStepSize, const Types::Coordinate finalStepSize )
  {
    this->m_FunctionalList.push_back( FunctionalWithStepSizes::SmartPtr( new FunctionalWithStepSizes( functional, initialStepSize, finalStepSize ) ) );
    return this->m_FunctionalList.size();
  }
 
  /** Perform the optimization.
   */
  virtual CallbackResult Optimize( CoordinateVector&, const Types::Coordinate = 1, const Types::Coordinate = 0 );

private:
  /// The actual low-level optimizer.
  Optimizer::SmartPtr m_Optimizer;

  /// Functional with step initial and final size.
  class FunctionalWithStepSizes
  {
  public:
    /// This class.
    typedef FunctionalWithStepSizes Self;

    /// Smart pointer.
    typedef SmartPointer<Self> SmartPtr;

    /// Default constructor.
    FunctionalWithStepSizes() :
      m_Functional( NULL ), m_InitialStepSize( 0.0 ), m_FinalStepSize( 0.0 ) 
    {}

    /// Constructor.
    FunctionalWithStepSizes( Functional::SmartPtr functional, const Types::Coordinate initialStepSize, const Types::Coordinate finalStepSize ) :
      m_Functional( functional ), 
      m_InitialStepSize( initialStepSize ), 
      m_FinalStepSize( finalStepSize ) 
    {}

    /// Smart pointer to functional.
    Functional::SmartPtr m_Functional;

    /// Initial step size.
    Types::Coordinate m_InitialStepSize;

    /// Final step size.
    Types::Coordinate m_FinalStepSize;
  };

  /// Type for list of functionals with step sizes.
  typedef std::list<FunctionalWithStepSizes::SmartPtr> FunctionalListType;

  /// List of functionals with step sizes.
  FunctionalListType m_FunctionalList;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMultiLevelOptimizer_h_included_
