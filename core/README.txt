
The Computational Morphometry Toolkit
==========================================================================


Release Notes -- CMTK 2.1.2
===========================

This is a bug fix release. Previous versions of CMTK do not compile with the
latest version of zlib (1.2.6); this is now fixed. This release also fixes
another problem with zlib on Mac, where the system zlib is not compatible with
the MacPorts zlib, which we need for Qt support and others.

This release also adds some new options to the "registrationx",
"groupwise_init", "groupwise_affine", and "lvote" tools.

Finally, this release introduces tools and library support for fitting
parametric transformations to deformation fields as well as linear
transformations to nonrigid transformations. The next major release of CMTK
will use this functionality to introduce CMTK support for directly-manipulated
free-form deformations (DM-FFD).

For a complete list of changes and fixes, see the CHANGELOG file.


CMTK has been built and tested on the following platforms:

- Linux 32bit (Fedora 15), gcc 4.6.1
- Linux 64bit (Fedora 15), gcc 4.6.2
- Linux 32bit (Fedora 16), gcc 4.6.2
- Linux 64bit (Fedora 16), gcc 4.6.2, CUDA 3.2
- Linux, i386, Oracle Solaris Studio 12.3 C++ 5.12 2011/11/16
- MacOSX 10.6, x86_64, gcc 4.2.1, CUDA 4.1
- MacOSX 10.6, x86_64, MacPorts gcc 4.6.2, CUDA 4.1
- MacOSX 10.6, x86_64, llvm-gcc-4.2.1
- MacOSX 10.6, x86_64, clang 2.0
- Cygwin, gcc 4.5.3
- Windows XP, VisualStudio 9 (2008 Express Edition), CUDA 3.2
- Windows XP, VisualStudio 10SP1 (2010 Express Edition), CUDA 4.1


Platform-Specific Issues
========================


MacOS-X
-------

- XCode gcc compiler on MacOS cannot build CMTK with OpenMP support enabled 
  at the same time. This is due to an Apple bug and has nothing to do with 
  CMTK per se.

  http://www.nitrc.org/tracker/index.php?func=detail&aid=5451&group_id=212&atid=877

  Workaround 1: use MacPorts gcc compiler.

  Workaround 2: set "CMTK_USE_OPENMP" configuration option to "OFF"


- Code coverage tests are only supported with gcc compiler and SDK 10.6. Older
  SDKs or the clang and llvm compiler front-ends do not support code coverage.

  http://www.nitrc.org/tracker/index.php?func=detail&aid=5450&group_id=212&atid=877


- To use a pre-compiled binary distribution of CMTK on MacOS, the following
  MacPorts packages have to be installed under /opt/local:
  - qt4
  - sqlite3
  - dcmtk


SolarisStudio Compiler, Linux/Intel
-----------------------------------

- SolarisStudio C++ 12.2 crashes when compiling the "Numerics" library with
  full optimization, -O3. A bug report has been filed with, and accepted by,
  Oracle:

  http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6989625

  Workaround: build "MinSizeRel" configuration, which sets optimization level
    to O2. Note that OpenMP must be disabled, because otherwise optimization
    is bumped back to O3 by default.

  This problem is also fixed in SolarisStudio 12.3, but see next issue.


- SolarisStudio C++ 12.3 crashes when compiling the file
  "cmtkEchoPlanarUnwarpFunctional.cxx" with OpenMP support turned on.

  Workaround: set "CMTK_USE_OPENMP" configuration option to "OFF"; this will
    lead to much of CMTK's functionality executing single-threaded.


Open64 Compiler
---------------

- CMTK does not build in Release mode with the Open64 compiler due to internal
  compilar errors.

==========================================================================

This software is available from

  http://www.nitrc.org/projects/cmtk/

==========================================================================
