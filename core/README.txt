
The Computational Morphometry Toolkit
==========================================================================


Release Notes -- CMTK 3.2.2
===========================

This is a minor bug fix release that corrects errors in the vtkxform tool, the
handling of inverse B-spline transformations in cmtk::XformList, and a
configuration issue with Greg Jefferis' munger.pl script.


Release Notes -- CMTK 3.2.1
===========================

ADDED FEATURES:

1. Compute "dimming" of lowest-intensity spheres in ADNI phantom scans

2. Global scale factor support in "vol2csv" tool


BUGFIXES:

1. Treat missing location information in Siemens DICOM CSA data structure as
   zero. This fixes incorrectly generated image-to-space transformation
   matrices.


Release Notes -- CMTK 3.2.0
===========================

ADDED FEATURES:

1. Paths for non-standard installation directories

CMTK_BINARY_DIR and CMTK_LIRBARY_DIR as used in a number of shell scripts
included with CMTK can now be overridden by setting them before executing one
such script. This should allow CMTK to 	work when it was installed/unpacked
into a path other than the installation path defined during build.


2. Additional information in XML sidecar files

The following additionl information is now written to XML sidecar files
created by the "dcm2image" tool:  

 - time units for TR, TE, TI, and dwell time.

 - sign of the phase encode direction from Siemens EPI files

 - slice times (relative to first slice in each volume) from Siemens DICOM
   files


3. New command line tools:

"geomatch" can compare grid dimensions, pixel sizes, and/or transformation
matrices of two or more images to ensure that they match.


4. New features in existing tools:

- The "sequence" tool can now optionally create histogram files in CSV format.

- Distance map-based erosion of multi-label fields, both in the library and
  exposed via the new "--erode-distance-multilabel"operation of "convertx"

- Exact isotropic resampling in UniformVolume class, exposed at command line
  via new "--resample-exact" operation of the "convertx" tool.

- The "epiunwarp" tool now has optional sixth positional parameter to export
  reverse deformation field to a given path.

- symmetric affine registration now also supports in-plane transformations.

- new operations in "convertx" tool, "--labels" and "--grey", set and reset
  the "intent_code" field in NIfTI images to NIFTI_INTENT_LABEL and 0, respectively.

- the "xform2dfield" tool can now write deformation fields in NIfTI format and
  fit for use by FSL's "applywarp" tool.

- freq_dim, phase_dim, slice_dim, as well as slice times and acquisition order
  fields are now properly set in NIFTI headers for files created by
  "dcm2image," so long as these data were available in the input DICOM files.

- The CNR centroid fallback (for phantoms with broken-off SNR sphere) can now
  be forced in the "detect_adni_phantom" tool, for cases where the SNR sphere
  is floating in a position that does not trigger the automatic fallback.


USER-VISIBLE CHANGES:

1. The "--average" operator of the "imagemath" tool no longer returns missing
values for a pixel where a single one of the input  images has missing
(padded) data.

2. Where appropriate, image processing operators now set the data class of
their outputs as "LABEL". This includes a variety of binarization/thresholding
techniques, label combination algorithms, and segmentations.

The "reformatx" tool now propagates floating image data class (label or grey)
to output, as represented by the "intent_code" header field when written in
NIfTI format.


Release Notes -- CMTK 3.1.0
===========================

ADDED FEATURES:

1. Polynomial Transformations

A new class of coordinate transformations based on polynomials (up to 4th
order) has been added. This is mostly useful for fitting transformations to
landmark correspondences derived from scans of the ADNI phantom.


2. Basic Support for Philips DWI DICOMs

Diffusion-related information (b-value, b-vector) are now extracted also from
recent Philips DICOM files. Since Philips DWI series appear to contain an
average diffusion weighted image (with b-value = 0 and b-vector = 0), the
dcm2image tool no longer writes the bVector tag to XML sidecar files if the
vector is zero but the b value is not. (This is so the necessary information
for tensor reconstruction can be easily extracted). In our experience, neither
Siemens nor GE write zero b-vectors with non-zero b-values, so behaviour
should not change for files from these manufacturers.


USER-VISIBLE CHANGES:

1. Global Scale Correction for Jacobian Maps

The "reformatx" tool no longer corrects Jacobian maps for global
scale by default. This seems like a more natural way to do this,
avoiding surprises.  To force the old behaviour, add
"--jacobian-correct-global" switch to the command  line.


2. "fview" Tool Command Line

The "fview" GUI no longer expects the fixed and moving image paths on the
command line (these are commonly and conveniently extracted from the
transformation). New, optional arguments allow providing or overriding the
deduced paths when necessary.


DEVELOPER CHANGES:

1. MacOS Configurations

Inconsistencies in pre-defined configurations for MacOS X builds have been
resolved and unnecessary external dependencies cleaned up.


Release Notes -- CMTK 3.0.0
===========================


CHANGES THAT BREAK BACKWARD COMPATIBILITY:
------------------------------------------

1. Affine Transformations:

Affine transformation matrices were broken when shear and non-isotropic scale
were used. Existing transformation files will continue to be read and generate
the same matrices as before, but these matrices will not strictly be
containing the specified scale and shear components. Newly-created
transformation files will use a different meaning of the "shear" coefficients.

https://www.nitrc.org/tracker/index.php?func=detail&aid=7179&group_id=212&atid=877


2. NIFTI Image Import/Export:

Handling of NIFTI qform and sform has changed in a way that breaks regression
tests, but should otherwise be benign. This should also make CMTK largely
NIFTI-compliant in the it puts image-to-physical coordinates into the header's
qform fields, rather than sform as before.

https://www.nitrc.org/tracker/index.php?func=detail&aid=7169&group_id=212&atid=877



OTHER USER-VISIBLE CHANGES:
---------------------------

1. Contributed pipeline scripts from the N-CANDA project, which briefly appeared
in the CMTK code tree, have been moved into their own, project-specific
repository.

2. unwarp_image_phantom default behaviour has been changed to multi-iteration
fitting. Also now supports residual-controlled fitting.

3. mk_analyze_hdr and mk_nifti_hdr default behaviour has changed - data type now
defaults to "byte" and orientation (for Analyze) now defaults to "axial",
rather than being "UNKNOWN"

4. dcm2image default behaviour has changed - potentially identifiable metadata is
no longer embedded in the "Description" field of NIFTI or Analyze images
created from DICOM files.

5. xform2scalar now puts padding pixels where application of transformation
sequence failed (e.g., due to failed numerical inversion of nonrigid
transformation).

6. fit_spline_xform was broken due to two bugs in the spline fitting code. These
have been fixed, but as a result, the tool now generates different output (as
it should, since the previous output was plain invalid).



Platform Support
================

CMTK has been built and tested on the following platforms:

- Linux 64bit (Fedora 20), gcc 4.8.2, CUDA 6.0 RC
- Linux 64bit (Fedora 20), clang 3.3
- Linux 64bit (Fedora 19), gcc 4.8.2
- Linux, i386, Oracle Solaris Studio 12.3 C++ 5.12 2011/11/16
- MacOSX 10.6, x86_64, gcc 4.2.1, CUDA 4.1
- MacOSX 10.6, x86_64, MacPorts gcc 4.8.2, CUDA 4.1
- MacOSX 10.6, x86_64, llvm-gcc-4.2.1
- MacOSX 10.6, x86_64, clang 2.0
- Cygwin, gcc 4.8.2
- Windows XP, VisualStudio 10SP1 (2010 Express Edition), CUDA 4.1
- Windows 7 x64, VisualStudio 2012 Express for Desktop, CUDA 5.5
- Windows 7 x64, VisualStudio 2013 Express for Desktop


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
