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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkAffineRegistrationCommandLine_h_included_
#define __cmtkAffineRegistrationCommandLine_h_included_

#include <cmtkconfig.h>

#include <cmtkAffineRegistration.h>

namespace 
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for command line-controlled affine registration.
 *@author T. Rohlfing
 */
class AffineRegistrationCommandLine : 
  /// Inherit generic affine registration.
  public AffineRegistration 
{
public:
  /// This class.
  typedef AffineRegistrationCommandLine Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Convenience typedef.
  typedef AffineRegistration Superclass;

  /** Constructor.
   *@param argc Number of command line arguments; this should be the argc
   * parameter of the main() function.
   *@param argv Array of command line arguments; this should be the argv
   * parameter of the main() function.
   *@exception ConstructorFailed This exception is thrown if there where
   * invalid or unknown options or missing required parameters. In all these
   * cases, an information text describing the known options will have been
   * written to the standard error stream before throwing the exception.
   */
  AffineRegistrationCommandLine ( int argc, char *argv [] );

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
   * This function write the transformation that was found to a studylist
   * archive with the name provided by command line arguments. The result is 
   * also printed to stderr in parameter list form.
   */
  virtual void OutputResult ( const CoordinateVector* ) const;

  /** Enter resolution level.
   * An information is printed to stderr and to the protocol file if one is
   * written.
   */
  virtual void EnterResolution( CoordinateVector::SmartPtr&, Functional::SmartPtr&, const int, const int );

private:
  /// Number of levels for automatic parameter generation.
  unsigned int m_AutoMultiLevels;

  /** Name of output studylist.
   * This is defined by the -o or --outlist command line option.
   */
  const char *Studylist;

  /** Name of the output matrix file.
    * This is defined by the "--out-matrix" command line argument.
    */
  const char* OutMatrixName;

  /** Name of the output parameter file.
    * This is defined by the "--out-params" command line argument.
    */
  const char* OutParametersName;

  /** Name of first study to be registered.
   * This is given as the first non-option command line paramter.
   */
  const char *Study1;

  /** Name of second study to be registered.
   * This is given as the second non-option command line paramter.
   */
  const char *Study2;

  /** Name of protocol output file.
   * This is defined by the -p or --protocol command line option.
   */
  const char *Protocol;

  /** Name of elapsed time output file.
   * This is defined by the -t or --time command line option.
   */
  const char *Time;

  /** Verbosity flag.
   * This is set to 'on' by -v or --verbose, and set to 'off' by -q or --quiet.
   */
  bool Verbose;

  /** Flag for initial center-of-mass translation.
   * This defaults to 'no' and is set to 'yes' by -i or --initxlate command
   * line switch.
   */
  bool InitXlate;

  /// Output result as matrix (text) file.
  void OutputResultMatrix( const char* matrixName ) const;

  /// Output result (and statistics) as studylist archive.
  void OutputResultParameters( const char* paramsName, const CoordinateVector& v ) const;

  /// Output result (and statistics) as studylist archive.
  void OutputResultList( const char* studyList ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef _COMMANDLINEVOXELREGISTRATION_H_

