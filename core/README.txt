
The Computational Morphometry Toolkit
==========================================================================

Release Notes -- CMTK 3.3.1
===========================

This is a minor bugfix release that addresses an issue with polynomial
transformations that could cause out-of-bounds memory access.


Release Notes -- CMTK 3.3.0
===========================

This release replaces signed 32bit integers for pixel indexing with signed
64bit integers, thus removing the previous limitation of allowing no more
than 2B pixels per image.


Platform Support
================

CMTK has been built and tested on the following platforms:

- Linux 64bit (Fedora 23), gcc 5.3.1
- Cygwin, gcc 4.9.3
- Windows7 x64, VisualStudio 2013 Express for Desktop


Platform-Specific Issues
========================


MacOS-X
-------

- Code coverage tests are only supported with gcc compiler and SDK 10.6. Older
  SDKs or the clang and llvm compiler front-ends do not support code coverage.

  http://www.nitrc.org/tracker/index.php?func=detail&aid=5450&group_id=212&atid=877


- To use a pre-combiled binary version of CMTK that was compiled on MacOS with
  MacPorts compilers, the following MacPorts packages have to be installed
  under /opt/local:
  - qt4
  - sqlite3
  - dcmtk


SolarisStudio Compiler, Linux/Intel
-----------------------------------

- SolarisStudio cannot build CMTK due to a non-standard, incompatible C++ STL.


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
