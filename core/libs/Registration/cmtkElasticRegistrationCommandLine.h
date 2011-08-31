/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkElasticRegistrationCommandLine_h_included_
#define __cmtkElasticRegistrationCommandLine_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkElasticRegistration.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/** Class for command line controlled  voxel registration.
 *\author T. Rohlfing
 *\version $Revision$ $Date$
 */
class ElasticRegistrationCommandLine : 
  /// Inherit generic elastic registration.
  public ElasticRegistration 
{
public:
  /// This class.
  typedef ElasticRegistrationCommandLine Self;

  /// Smart pointer.
  typedef SmartPointer<ElasticRegistrationCommandLine> SmartPtr;

  /// Parent class.
  typedef ElasticRegistration Superclass;

  /** Constructor.
   *\param argc Number of command line arguments; this should be the argc
   * parameter of the main() function.
   *\param argv Array of command line arguments; this should be the argv
   * parameter of the main() function.
   */
  ElasticRegistrationCommandLine ( const int argc, const char *argv[] );

  /// Destructor.
  ~ElasticRegistrationCommandLine ();

  /** Perform registration.
   */
  virtual CallbackResult Register ();

protected:
  /** Initialize registration.
   * So far, this function has no effect other than calling the equivalent
   * inherited function.
   */
  virtual CallbackResult InitRegistration();

  /** Output registration result.
   */
  virtual void OutputResult ( const CoordinateVector* v, const CallbackResult irq = CALLBACK_OK );
  
  /** Enter resolution level.
   * An information is printed to stderr and to the protocol file if one is
   * written.
   */
  virtual void EnterResolution( CoordinateVector::SmartPtr&, Functional::SmartPtr&, const int, const int );

  /** Leave resolution level.
   * The transformation found so far is written to a file if desired.
   */
  virtual int DoneResolution( CoordinateVector::SmartPtr&, Functional::SmartPtr&, const int, const int );

private:
  /** Name of input studylist.
   * This is defined by the --inlist command line parameter (not currently
   * supported).
   */
  const char *InputStudylist;

  /** Name of output studylist.
   * This is defined by the -o or --outlist command line option.
   */
  const char *Studylist;

  /** Name of first study to be registered.
   * This is given as the first non-option command line paramter.
   */
  char *Study1;

  /** Name of second study to be registered.
   * This is given as the second non-option command line paramter.
   */
  char *Study2;
  
#ifdef CMTK_USE_SQLITE
  /// Database to update after registration completes.
  const char* m_UpdateDB;
#endif

  /// Filename for rigidity constraint map.
  const char* RigidityConstraintMapFilename;

  /** Name of protocol output file.
   * This is defined by the -p or --protocol command line option.
   */
  const char *Protocol;

  /** Name of elapsed time output file.
   * This is defined by the -t or --time command line option.
   */
  const char *Time;

  /** Select whether too create intermediate warp output files (level-xx.list).
   */
  bool m_OutputIntermediate;

  /// Write deformation to studylist archive.
  void OutputWarp ( const char* ) const;

  /// Name of output transformation file in ITK format.
  const char* m_OutputPathITK;

  /// Path for reformatted floating image.
  const char* m_ReformattedImagePath;

public:
  /// Static pointer to this object.
  static Self* StaticThis;

  /// Counter for intermediate result files.
  int IntermediateResultIndex;

  /// Write intermediate deformation file.
  void OutputIntermediate( const bool incrementCount = true );
};

//@}

} // namespace cmtk

/// Signal handler that writes intermediate result during a level.
extern "C" void cmtkElasticRegistrationCommandLineDispatchSIGUSR1( int sig );



#endif // #ifndef __cmtkElasticRegistrationCommandLine_h_included_
