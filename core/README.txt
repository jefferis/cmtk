
The Computational Morphometry Toolkit
==========================================================================


Release Notes -- CMTK 1.7.2
===========================

For a complete list of changes and fixes, see the CHANGELOG file.

CMTK has been built and tested on the following platforms:

- Linux 32bit (Fedora 13), gcc 4.4.5, CUDA 3.2
- Linux 32bit (Fedora 15), gcc 4.6.0
- Linux 64bit (Fedora 13), gcc 4.4.5, CUDA 3.2
- Linux 64bit (Fedora 14), gcc 4.5.1
- Linux, i386, Oracle Solaris Studio 12.2 C++ 5.11 2010/08/13
- MacOSX 10.4, i386, gcc 4.0.1
- MacOSX 10.6, x86_64, gcc 4.2.1, CUDA 3.2
- MacOSX 10.6, x86_64, llvm-gcc-4.2.1
- MacOSX 10.6, x86_64, clang 2.0
- Cygwin, gcc 4.3.4
- Windows XP, VisualStudio 9 (2008 Express Edition), CUDA 3.2
- Windows XP, VisualStudio 10SP1 (2010 Express Edition), CUDA 4.0
- OpenSolaris 11 Express, SunStudio 12.1, i386 (CC 5.10)
- OpenSolaris 11 Express, SunStudio 12.1, x86_64 (CC 5.10)
- OpenSolaris 11 Express, Oracle/SunStudio 12.2, x86_64


Platform-Specific Issues
========================


MacOS-X
-------

- Compilers on MacOS cannot build CMTK with shared libraries and OpenMP
  support enabled at the same time. This is due to an Apple bug and has
  nothing to do with CMTK per se. 

  http://www.nitrc.org/tracker/index.php?func=detail&aid=5451&group_id=212&atid=877

  Workaround: build CMTK with static libraries.


- Code coverage tests are only supported with gcc compiler and SDK 10.6. Older
  SDKs or the clang and llvm compiler front-ends do not support code coverage.

  http://www.nitrc.org/tracker/index.php?func=detail&aid=5450&group_id=212&atid=877


==========================================================================

This software is available from

  http://www.nitrc.org/projects/cmtk/

==========================================================================
