
The Computational Morphometry Toolkit
==========================================================================

Release Notes -- CMTK 1.6.0
===========================

This is a maintenance, consolidation, and feature release of CMTK. Note that
this release breaks backward compatibility by eliminating some obsolete and
redundant command line tools:

1. Groupwise registration tools have been consolidated. The new
   "groupwise_affine" tool replaces "congeal" and "groupwise_rmi", and the new
   "groupwise_warp" tool replaces "congeal_warp" and "groupwise_rmi_warp".

2. The old "convert" tool has been removed. The "convertx" tool should now be
   used instead.

3. The old "probe_xform" tool has been removed. Use "gregxform" instead.


The most notable improvements are:

1. Volume injection and reconstruction tools now preserve the physical coordinate
   space of the input images.

2. The "mk_phantom_3d" tool can now determine the location of an MRS voxel
   from GE DICOM files and draw the voxel in the correct image location.


For a complete list of changes and fixes, see the CHANGELOG file.


CMTK has been built and tested on the following platforms:

- Linux 32bit (Fedora 13), gcc 4.4.5, CUDA 3.2
- Linux 64bit (Fedora 13), gcc 4.4.5, CUDA 3.2
- Linux 64bit (Fedora 14), gcc 4.5.1
- Linux, i386, Oracle/SunStudio C++ 5.11 (Express June 2010)
- MacOSX 10.4, i386, gcc 4.0.1
- MacOSX 10.6, x86_64, gcc 4.2.1, CUDA 3.2
- Cygwin, gcc 4.5.0
- Windows XP, VisualStudio 9 (2008 Express Edition), CUDA 3.2
- Windows XP, VisualStudio 10 (2010 Express Edition)
- OpenSolaris, SunStudio 12.1, i386 (CC 5.10)
- OpenSolaris, SunStudio 12.1, x86_64 (CC 5.10)
- OpenSolaris, Oracle/SunStudio 12.2, x86_64


Note that as of this release, pre-built binaries for the MacOS/PPC platform
are no longer provided.

==========================================================================

This software is available from

  http://www.nitrc.org/projects/cmtk/

==========================================================================
