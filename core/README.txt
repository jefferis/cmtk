
The Computational Morphometry Toolkit
==========================================================================


Release Notes -- CMTK 1.7.0
===========================

This is a feature and bugfix release of CMTK.

CMTK now employs a toolkit-wide system of multi-level verbose output, which
replaces the former, tool-specific output. Consistency of verbose output and
streams used for this purpose has been futher improved.

This release also greatly improves CMTK's reliability on different platforms.

On the Windows platform, a long-existing bug has been fixed that broke
parallel computation. Also, CMTK now supports application-level automated
testing, driven by Cygwin-supplied "sh" shell. In the process of setting up
testing, numerous Windows-specific bugs and problems have been identified and
fixed. In particular, CMTK now compiles using Visual C++ with OpenMP parallel
processing support enabled.

On the MacOS platform, this release introduces support for Apple's Grand
Central Dispatch parallel processing framework. This also replaces, in part,
the OpenMP parallelization on MacOS, which is broken with shared library
builds (Thanks a lot, Apple!)

For a complete list of changes and fixes, see the CHANGELOG file.

CMTK has been built and tested on the following platforms:

- Linux 32bit (Fedora 13), gcc 4.4.5, CUDA 3.2
- Linux 64bit (Fedora 13), gcc 4.4.5, CUDA 3.2
- Linux 64bit (Fedora 14), gcc 4.5.1
- Linux, i386, Oracle/SunStudio C++ 5.11 (Express June 2010)
- MacOSX 10.4, i386, gcc 4.0.1
- MacOSX 10.6, x86_64, gcc 4.2.1, CUDA 3.2
- MacOSX 10.6, x86_64, llvm-gcc-4.2.1
- MacOSX 10.6, x86_64, clang 2.0
- Cygwin, gcc 4.3.4
- Windows XP, VisualStudio 9 (2008 Express Edition), CUDA 3.2
- Windows XP, VisualStudio 10SP1 (2010 Express Edition), CUDA 3.2
- OpenSolaris, SunStudio 12.1, i386 (CC 5.10)
- OpenSolaris, SunStudio 12.1, x86_64 (CC 5.10)
- OpenSolaris, Oracle/SunStudio 12.2, x86_64


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
