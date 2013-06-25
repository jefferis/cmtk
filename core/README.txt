
The Computational Morphometry Toolkit
==========================================================================

Release Notes -- CMTK 2.4.0
===========================

Contributed pipeline scripts from the N-CANDA project, which briefly appeared
in the CMTK code tree, have been moved into their own, project-specific
repository.

unwarp_image_phantom default behaviour has been changed to multi-iteration
fitting. Also now supports residual-controlled fitting.


Platform Support
================

CMTK has been built and tested on the following platforms:

- Linux 64bit (Fedora 19), gcc 4.8.1, CUDA 3.2
- Linux 64bit (Fedora 19), clang 3.3
- Linux 64bit (Fedora 17), gcc 4.7.2, CUDA 3.2
- Linux, i386, Oracle Solaris Studio 12.3 C++ 5.12 2011/11/16
- MacOSX 10.6, x86_64, gcc 4.2.1, CUDA 4.1
- MacOSX 10.6, x86_64, MacPorts gcc 4.7.2, CUDA 4.1
- MacOSX 10.6, x86_64, llvm-gcc-4.2.1
- MacOSX 10.6, x86_64, clang 2.0
- Cygwin, gcc 4.5.3
- Windows XP, VisualStudio 10SP1 (2010 Express Edition), CUDA 4.1
- Windows 7 64bit, VisualStudio 2012 Express for Desktop, CUDA 4.1


Platform-Specific Issues
========================


MacOS-X
-------

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
  compiler errors.


==========================================================================

This software is available from

  http://www.nitrc.org/projects/cmtk/

==========================================================================
