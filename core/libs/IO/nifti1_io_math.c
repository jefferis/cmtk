#include "nifti1_io_math.h"   /* typedefs, prototypes, macros, etc. */

/*****===================================================================*****/
/*****     Sample functions to deal with NIFTI-1 and ANALYZE files       *****/
/*****...................................................................*****/
/*****            This code is released to the public domain.            *****/
/*****...................................................................*****/
/*****  Author: Robert W Cox, SSCC/DIRP/NIMH/NIH/DHHS/USA/EARTH          *****/
/*****  Date:   August 2003                                              *****/
/*****...................................................................*****/
/*****  Neither the National Institutes of Health (NIH), nor any of its  *****/
/*****  employees imply any warranty of usefulness of this software for  *****/
/*****  any purpose, and do not assume any liability for damages,        *****/
/*****  incidental or otherwise, caused by any use of this document.     *****/
/*****===================================================================*****/

/** \file nifti1_io.c
    \brief main collection of nifti1 i/o routines
           - written by Bob Cox, SSCC NIMH
           - revised by Mark Jenkinson, FMRIB
           - revised by Rick Reynolds, SSCC, NIMH
           - revised by Kate Fissell, University of Pittsburgh

        The library history can be viewed via "nifti_tool -nifti_hist".
    <br>The library version can be viewed via "nifti_tool -nifti_ver".
 */

/*! global history and version strings, for printing */
static char const * const gni_history[] =
{
  "----------------------------------------------------------------------\n"
  "history (of nifti library changes):\n"
  "\n",
  "0.0  August, 2003 [rwcox]\n"
  "     (Robert W Cox of the National Institutes of Health, SSCC/DIRP/NIMH)\n"
  "   - initial version\n"
  "\n",
  "0.1  July/August, 2004 [Mark Jenkinson]\n"
  "     (FMRIB Centre, University of Oxford, UK)\n"
  "   - Mainly adding low-level IO and changing things to allow gzipped\n"
  "     files to be read and written\n"
  "   - Full backwards compatability should have been maintained\n"
  "\n",
  "0.2  16 Nov 2004 [rickr]\n"
  "     (Rick Reynolds of the National Institutes of Health, SSCC/DIRP/NIMH)\n"
  "   - included Mark's changes in the AFNI distribution (including znzlib/)\n"
  "     (HAVE_ZLIB is commented out for the standard distribution)\n"
  "   - modified nifti_validfilename() and nifti_makebasename()\n"
  "   - added nifti_find_file_extension()\n"
  "\n",
  "0.3  3 Dec 2004 [rickr]\n"
  "   - note: header extensions are not yet checked for\n"
  "   - added formatted history as global string, for printing\n"
  "   - added nifti_disp_lib_hist(), to display the nifti library history\n"
  "   - added nifti_disp_lib_version(), to display the nifti library history\n",
  "   - re-wrote nifti_findhdrname()\n"
  "       o used nifti_find_file_extension()\n"
  "       o changed order of file tests (default is .nii, depends on input)\n"
  "       o free hdrname on failure\n"
  "   - made similar changes to nifti_findimgname()\n"
  "   - check for NULL return from nifti_findhdrname() calls\n",
  "   - removed most of ERREX() macros\n"
  "   - modified nifti_image_read()\n"
  "       o added debug info and error checking (on gni_debug > 0, only)\n"
  "       o fail if workingname is NULL\n"
  "       o check for failure to open header file\n"
  "       o free workingname on failure\n"
  "       o check for failure of nifti_image_load()\n"
  "       o check for failure of nifti_convert_nhdr2nim()\n",
  "   - changed nifti_image_load() to int, and check nifti_read_buffer return\n"
  "   - changed nifti_read_buffer() to fail on short read, and to count float\n"
  "     fixes (to print on debug)\n"
  "   - changed nifti_image_infodump to print to stderr\n"
  "   - updated function header comments, or moved comments above header\n"
  "   - removed const keyword\n"
  "   - added LNI_FERR() macro for error reporting on input files\n"
  "\n",
  "0.4  10 Dec 2004 [rickr]  - added header extensions\n"
  "   - in nifti1_io.h:\n"
  "       o added num_ext and ext_list to the definition of nifti_image\n"
  "       o made many functions static (more to follow)\n"
  "       o added LNI_MAX_NIA_EXT_LEN, for max nifti_type 3 extension length\n",
  "   - added __DATE__ to version output in nifti_disp_lib_version()\n"
  "   - added nifti_disp_matrix_orient() to print orientation information\n"
  "   - added '.nia' as a valid file extension in nifti_find_file_extension()\n"
  "   - added much more debug output\n"
  "   - in nifti_image_read(), in the case of an ASCII header, check for\n"
  "     extensions after the end of the header\n",
  "   - added nifti_read_extensions() function\n"
  "   - added nifti_read_next_extension() function\n"
  "   - added nifti_add_exten_to_list() function\n"
  "   - added nifti_check_extension() function\n"
  "   - added nifti_write_extensions() function\n"
  "   - added nifti_extension_size() function\n"
  "   - in nifti_set_iname_offest():\n"
  "       o adjust offset by the extension size and the extender size\n",
  "       o fixed the 'ceiling modulo 16' computation\n"
  "   - in nifti_image_write_hdr_img2(): \n"
  "       o added extension writing\n"
  "       o check for NULL return from nifti_findimgname()\n"
  "   - include number of extensions in nifti_image_to_ascii() output\n"
  "   - in nifti_image_from_ascii():\n"
  "       o return bytes_read as a parameter, computed from the final spos\n"
  "       o extract num_ext from ASCII header\n"
  "\n",
  "0.5  14 Dec 2004 [rickr]  - added sub-brick reading functions\n"
  "   - added nifti_brick_list type to nifti1_io.h, along with new prototypes\n"
  "   - added main nifti_image_read_bricks() function, with description\n"
  "   - added nifti_image_load_bricks() - library function (requires nim)\n"
  "   - added valid_nifti_brick_list() - library function\n"
  "   - added free_NBL() - library function\n",
  "   - added update_nifti_image_for_brick_list() for dimension update\n"
  "   - added nifti_load_NBL_bricks(), nifti_alloc_NBL_mem(),\n"
  "           nifti_copynsort() and force_positive() (static functions)\n"
  "   - in nifti_image_read(), check for failed load only if read_data is set\n"
  "   - broke most of nifti_image_load() into nifti_image_load_prep()\n"
  "\n",
  "0.6  15 Dec 2004 [rickr]  - added sub-brick writing functionality\n"
  "   - in nifti1_io.h, removed znzlib directory from include - all nifti\n"
  "       library files are now under the nifti directory\n"
  "   - nifti_read_extensions(): print no offset warning for nifti_type 3\n"
  "   - nifti_write_all_data():\n"
  "       o pass nifti_brick_list * NBL, for optional writing\n"
  "       o if NBL, write each sub-brick, sequentially\n",
  "   - nifti_set_iname_offset(): case 1 must have sizeof() cast to int\n"
  "   - pass NBL to nifti_image_write_hdr_img2(), and allow NBL or data\n"
  "   - added nifti_image_write_bricks() wrapper for ...write_hdr_img2()\n"
  "   - included compression abilities\n"
  "\n",
  "0.7  16 Dec 2004 [rickr] - minor changes to extension reading\n"
  "\n",
  "0.8  21 Dec 2004 [rickr] - restrict extension reading, and minor changes\n"
  "   - in nifti_image_read(), compute bytes for extensions (see remaining)\n"
  "   - in nifti_read_extensions(), pass 'remain' as space for extensions,\n"
  "        pass it to nifti_read_next_ext(), and update for each one read \n"
  "   - in nifti_check_extension(), require (size <= remain)\n",
  "   - in update_nifti_image_brick_list(), update nvox\n"
  "   - in nifti_image_load_bricks(), make explicit check for nbricks <= 0\n"
  "   - in int_force_positive(), check for (!list)\n"
  "   - in swap_nifti_header(), swap sizeof_hdr, and reorder to struct order\n"
  "   - change get_filesize functions to signed ( < 0 is no file or error )\n",
  "   - in nifti_validfilename(), lose redundant (len < 0) check\n"
  "   - make print_hex_vals() static\n"
  "   - in disp_nifti_1_header, restrict string field widths\n"
  "\n",
  "0.9  23 Dec 2004 [rickr] - minor changes\n"
  "   - broke ASCII header reading out of nifti_image_read(), into new\n"
  "        functions has_ascii_header() and read_ascii_image()\n",
  "   - check image_read failure and znzseek failure\n"
  "   - altered some debug output\n"
  "   - nifti_write_all_data() now returns an int\n"
  "\n",
  "0.10 29 Dec 2004 [rickr]\n"
  "   - renamed nifti_valid_extension() to nifti_check_extension()\n"
  "   - added functions nifti_makehdrname() and nifti_makeimgname()\n"
  "   - added function valid_nifti_extensions()\n"
  "   - in nifti_write_extensions(), check for validity before writing\n",
  "   - rewrote nifti_image_write_hdr_img2():\n"
  "       o set write_data and leave_open flags from write_opts\n"
  "       o add debug print statements\n"
  "       o use nifti_write_ascii_image() for the ascii case\n"
  "       o rewrote the logic of all cases to be easier to follow\n",
  "   - broke out code as nifti_write_ascii_image() function\n"
  "   - added debug to top-level write functions, and free the znzFile\n"
  "   - removed unused internal function nifti_image_open()\n"
  "\n",
  "0.11 30 Dec 2004 [rickr] - small mods\n"
  "   - moved static function prototypes from header to C file\n"
  "   - free extensions in nifti_image_free()\n"
  "\n",
  "1.0  07 Jan 2005 [rickr] - INITIAL RELEASE VERSION\n"
  "   - added function nifti_set_filenames()\n"
  "   - added function nifti_read_header()\n"
  "   - added static function nhdr_looks_good()\n"
  "   - added static function need_nhdr_swap()\n"
  "   - exported nifti_add_exten_to_list symbol\n",
  "   - fixed #bytes written in nifti_write_extensions()\n"
  "   - only modify offset if it is too small (nifti_set_iname_offset)\n"
  "   - added nifti_type 3 to nifti_makehdrname and nifti_makeimgname\n"
  "   - added function nifti_set_filenames()\n"
  "\n",
  "1.1  07 Jan 2005 [rickr]\n"
  "   - in nifti_read_header(), swap if needed\n"
  "\n",
  "1.2  07 Feb 2005 [kate fissell c/o rickr] \n"
  "   - nifti1.h: added doxygen comments for main struct and #define groups\n"
  "   - nifti1_io.h: added doxygen comments for file and nifti_image struct\n"
  "   - nifti1_io.h: added doxygen comments for file and some functions\n"
  "   - nifti1_io.c: changed nifti_copy_nim_info to use memcpy\n"
  "\n",
  "1.3  09 Feb 2005 [rickr]\n"
  "   - nifti1.h: added doxygen comments for extension structs\n"
  "   - nifti1_io.h: put most #defines in #ifdef _NIFTI1_IO_C_ block\n"
  "   - added a doxygen-style description to every exported function\n"
  "   - added doxygen-style comments within some functions\n"
  "   - re-exported many znzFile functions that I had made static\n"
  "   - re-added nifti_image_open (sorry, Mark)\n"
  "   - every exported function now has 'nifti' in the name (19 functions)\n",
  "   - made sure every alloc() has a failure test\n"
  "   - added nifti_copy_extensions function, for use in nifti_copy_nim_info\n"
  "   - nifti_is_gzfile: added initial strlen test\n"
  "   - nifti_set_filenames: added set_byte_order parameter option\n"
  "     (it seems appropriate to set the BO when new files are associated)\n"
  "   - disp_nifti_1_header: prints to stdout (a.o.t. stderr), with fflush\n"
  "\n",
  "1.4  23 Feb 2005 [rickr] - sourceforge merge\n"
  "   - merged into the nifti_io CVS directory structure at sourceforge.net\n"
  "   - merged in 4 changes by Mark, and re-added his const keywords\n"
  "   - cast some pointers to (void *) for -pedantic compile option\n"
  "   - added nifti_free_extensions()\n"
  "\n",
  "1.5  02 Mar 2005 [rickr] - started nifti global options\n"
  "   - gni_debug is now g_opts.debug\n"
  "   - added validity check parameter to nifti_read_header\n"
  "   - need_nhdr_swap no longer does test swaps on the stack\n"
  "\n",
  "1.6  05 April 2005 [rickr] - validation and collapsed_image_read\n"
  "   - added nifti_read_collapsed_image(), an interface for reading partial\n"
  "     datasets, specifying a subset of array indices\n"
  "   - for read_collapsed_image, added static functions: rci_read_data(),\n"
  "     rci_alloc_mem(), and make_pivot_list()\n",
  "   - added nifti_nim_is_valid() to check for consistency (more to do)\n"
  "   - added nifti_nim_has_valid_dims() to do many dimensions tests\n"
  "\n",
  "1.7  08 April 2005 [rickr]\n"
  "   - added nifti_update_dims_from_array() - to update dimensions\n"
  "   - modified nifti_makehdrname() and nifti_makeimgname():\n"
  "       if prefix has a valid extension, use it (else make one up)\n"
  "   - added nifti_get_intlist - for making an array of ints\n"
  "   - fixed init of NBL->bsize in nifti_alloc_NBL_mem()  {thanks, Bob}\n"
  "\n",
  "1.8  14 April 2005 [rickr]\n"
  "   - added nifti_set_type_from_names(), for nifti_set_filenames()\n"
  "     (only updates type if number of files does not match it)\n"
  "   - added is_valid_nifti_type(), just to be sure\n"
  "   - updated description of nifti_read_collapsed_image() for *data change\n"
  "     (if *data is already set, assume memory exists for results)\n"
  "   - modified rci_alloc_mem() to allocate only if *data is NULL\n"
  "\n",
  "1.9  19 April 2005 [rickr]\n"
  "   - added extension codes NIFTI_ECODE_COMMENT and NIFTI_ECODE_XCEDE\n"
  "   - added nifti_type codes NIFTI_MAX_ECODE and NIFTI_MAX_FTYPE\n"
  "   - added nifti_add_extension() {exported}\n"
  "   - added nifti_fill_extension() as a static function\n"
  "   - added nifti_is_valid_ecode() {exported}\n",
  "   - nifti_type values are now NIFTI_FTYPE_* file codes\n"
  "   - in nifti_read_extensions(), decrement 'remain' by extender size, 4\n"
  "   - in nifti_set_iname_offset(), case 1, update if offset differs\n"
  "   - only output '-d writing nifti file' if debug > 1\n"
  "\n",
  "1.10 10 May 2005 [rickr]\n"
  "   - files are read using ZLIB only if they end in '.gz'\n"
  "\n",
  "1.11 12 August 2005 [kate fissell]\n"
  "   - Kate's 0.2 release packaging, for sourceforge\n"
  "\n",
  "1.12 17 August 2005 [rickr] - comment (doxygen) updates\n"
  "   - updated comments for most functions (2 updates from Cinly Ooi)\n"
  "   - added nifti_type_and_names_match()\n"
  "\n",
  "1.12a 24 August 2005 [rickr] - remove all tabs from Clibs/*/*.[ch]\n",
  "1.12b 25 August 2005 [rickr] - changes by Hans Johnson\n",
  "1.13  25 August 2005 [rickr]\n",
  "   - finished changes by Hans for Insight\n"
  "   - added const in all appropraite parameter locations (30-40)\n"
  "     (any pointer referencing data that will not change)\n"
  "   - shortened all string constants below 509 character limit\n"
  "1.14  28 October 2005 [HJohnson]\n",
  "   - use nifti_set_filenames() in nifti_convert_nhdr2nim()\n"
  "1.15  02 November 2005 [rickr]\n",
  "   - added skip_blank_ext to nifti_global_options\n"
  "   - added nifti_set_skip_blank_ext(), to set option\n"
  "   - if skip_blank_ext and no extensions, do not read/write extender\n"
  "1.16 18 November 2005 [rickr]\n",
  "   - removed any test or access of dim[i], i>dim[0]\n"
  "   - do not set pixdim for collapsed dims to 1.0, leave them as they are\n"
  "   - added magic and dim[i] tests in nifti_hdr_looks_good()\n"
  "   - added 2 size_t casts\n"
  "1.17 22 November 2005 [rickr]\n",
  "   - in hdr->nim, for i > dim[0], pass 0 or 1, else set to 1\n"
  "1.18 02 March 2006 [rickr]\n",
  "   - in nifti_alloc_NBL_mem(), fixed nt=0 case from 1.17 change\n"
  "1.19 23 May 2006 [HJohnson,rickr]\n",
  "   - nifti_write_ascii_image(): free(hstr)\n"
  "   - nifti_copy_extensions(): clear num_ext and ext_list\n"
  "1.20 27 Jun 2006 [rickr]\n",
  "   - nifti_findhdrname(): fixed assign of efirst to match stated logic\n"
  "     (problem found by Atle Bjørnerud)\n"
  "1.21 05 Sep 2006 [rickr] update for nifticlib-0.4 release\n",
  "   - was reminded to actually add nifti_set_skip_blank_ext()\n"
  "   - init g_opts.skip_blank_ext to 0\n"
  "1.22 01 Jun 2007 nifticlib-0.5 release\n",
  "1.23 05 Jun 2007 nifti_add_exten_to_list: revert on failure, free old list\n"
  "1.24 07 Jun 2007 nifti_copy_extensions: use esize-8 for data size\n"
  "1.25 12 Jun 2007 [rickr] EMPTY_IMAGE creation\n",
  "   - added nifti_make_new_header() - to create from dims/dtype\n"
  "   - added nifti_make_new_nim() - to create from dims/dtype/fill\n"
  "   - added nifti_is_valid_datatype(), and more debug info\n",
  "1.26 27 Jul 2007 [rickr] handle single volumes > 2^31 bytes (but < 2^32)\n",
  "1.27 28 Jul 2007 [rickr] nim->nvox, NBL-bsize are now type size_t\n"
  "1.28 30 Jul 2007 [rickr] size_t updates\n",
  "1.29 08 Aug 2007 [rickr] for list, valid_nifti_brick_list requires 3 dims\n"
  "1.30 08 Nov 2007 [Yaroslav/rickr]\n"
  "   - fix ARM struct alignment problem in byte-swapping routines\n",
  "1.31 29 Nov 2007 [rickr] for nifticlib-1.0.0\n"
  "   - added nifti_datatype_to/from_string routines\n"
  "   - added DT_RGBA32/NIFTI_TYPE_RGBA32 datatype macros (2304)\n"
  "   - added NIFTI_ECODE_FREESURFER (14)\n",
  "1.32 08 Dec 2007 [rickr]\n"
  "   - nifti_hdr_looks_good() allows ANALYZE headers (req. by V. Luccio)\n"
  "   - added nifti_datatype_is_valid()\n",
  "1.33 05 Feb 2008 [hansj,rickr] - block nia.gz use\n"
  "1.34 13 Jun 2008 [rickr] - added nifti_compiled_with_zlib()\n"
  "1.35 03 Aug 2008 [rickr]\n",
  "   - deal with swapping, so that CPU type does not affect output\n"
  "     (motivated by C Burns)\n"
  "   - added nifti_analyze75 structure and nifti_swap_as_analyze()\n"
  "   - previous swap_nifti_header is saved as old_swap_nifti_header\n"
  "   - also swap UNUSED fields in nifti_1_header struct\n",
  "1.36 07 Oct 2008 [rickr]\n",
  "   - added nifti_NBL_matches_nim() check for write_bricks()\n"
  "1.37 10 Mar 2009 [rickr]\n",
  "   - H Johnson cast updates (06 Feb)\n"
  "   - added NIFTI_ECODE_PYPICKLE for PyNIfTI (06 Feb)\n"
  "   - added NIFTI_ECODEs 18-28 for the LONI MiND group\n"
  "1.38 28 Apr 2009 [rickr]\n",
  "   - uppercase extensions are now valid (requested by M. Coursolle)\n"
  "   - nifti_set_allow_upper_fext controls this option (req by C. Ooi)\n"
  "1.39 23 Jun 2009 [rickr]: added 4 checks of alloc() returns\n",
  "1.40 16 Mar 2010 [rickr]: added NIFTI_ECODE_VOXBO for D. Kimberg\n",
  "1.41 28 Apr 2010 [rickr]: added NIFTI_ECODE_CARET for J. Harwell\n",
  "1.42 06 Jul 2010 [rickr]: trouble with large (gz) files\n",
  "   - noted/investigated by M Hanke and Y Halchenko\n"
  "   - fixed znzread/write, noting example by M Adler\n"
  "   - changed nifti_swap_* routines/calls to take size_t (6)\n"
  "1.43 07 Jul 2010 [rickr]: fixed znzR/W to again return nmembers\n",
  "1.44 19 Jul 2013 [rickr]: ITK compatibility updates from H Johnson\n",
  "----------------------------------------------------------------------\n"
};
static const char gni_version[] = "nifti library version 1.44 (19 July, 2013)";

/*---------------------------------------------------------------------------*/
/*! Given the quaternion parameters (etc.), compute a transformation matrix.

   See comments in nifti1.h for details.
     - qb,qc,qd = quaternion parameters
     - qx,qy,qz = offset parameters
     - dx,dy,dz = grid stepsizes (non-negative inputs are set to 1.0)
     - qfac     = sign of dz step (< 0 is negative; >= 0 is positive)

   <pre>
   If qx=qy=qz=0, dx=dy=dz=1, then the output is a rotation matrix.
   For qfac >= 0, the rotation is proper.
   For qfac <  0, the rotation is improper.
   </pre>

   \see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h
   \see nifti_mat44_to_quatern, nifti_make_orthog_mat44,
       nifti_mat44_to_orientation

*//*-------------------------------------------------------------------------*/
mat44 nifti_quatern_to_mat44( float qb, float qc, float qd,
                              float qx, float qy, float qz,
                              float dx, float dy, float dz, float qfac )
{
   mat44 R ;
   double a,b=qb,c=qc,d=qd , xd,yd,zd ;

   /* last row is always [ 0 0 0 1 ] */

   R.m[3][0]=R.m[3][1]=R.m[3][2] = 0.0f ; R.m[3][3]= 1.0f ;

   /* compute a parameter from b,c,d */

   a = 1.0l - (b*b + c*c + d*d) ;
   if( a < 1.e-7l ){                   /* special case */
     a = 1.0l / sqrt(b*b+c*c+d*d) ;
     b *= a ; c *= a ; d *= a ;        /* normalize (b,c,d) vector */
     a = 0.0l ;                        /* a = 0 ==> 180 degree rotation */
   } else{
     a = sqrt(a) ;                     /* angle = 2*arccos(a) */
   }

   /* load rotation matrix, including scaling factors for voxel sizes */

   xd = (dx > 0.0) ? dx : 1.0l ;       /* make sure are positive */
   yd = (dy > 0.0) ? dy : 1.0l ;
   zd = (dz > 0.0) ? dz : 1.0l ;

   if( qfac < 0.0 ) zd = -zd ;         /* left handedness? */

   R.m[0][0] = (float)( (a*a+b*b-c*c-d*d) * xd) ;
   R.m[0][1] = 2.0l * (b*c-a*d        ) * yd ;
   R.m[0][2] = 2.0l * (b*d+a*c        ) * zd ;
   R.m[1][0] = 2.0l * (b*c+a*d        ) * xd ;
   R.m[1][1] = (float)( (a*a+c*c-b*b-d*d) * yd) ;
   R.m[1][2] = 2.0l * (c*d-a*b        ) * zd ;
   R.m[2][0] = 2.0l * (b*d-a*c        ) * xd ;
   R.m[2][1] = 2.0l * (c*d+a*b        ) * yd ;
   R.m[2][2] = (float)( (a*a+d*d-c*c-b*b) * zd) ;

   /* load offsets */

   R.m[0][3] = qx ; R.m[1][3] = qy ; R.m[2][3] = qz ;

   return R ;
}

/*---------------------------------------------------------------------------*/
/*! Given the 3x4 upper corner of the matrix R, compute the quaternion
   parameters that fit it.

     - Any NULL pointer on input won't get assigned (e.g., if you don't want
       dx,dy,dz, just pass NULL in for those pointers).
     - If the 3 input matrix columns are NOT orthogonal, they will be
       orthogonalized prior to calculating the parameters, using
       the polar decomposition to find the orthogonal matrix closest
       to the column-normalized input matrix.
     - However, if the 3 input matrix columns are NOT orthogonal, then
       the matrix produced by nifti_quatern_to_mat44 WILL have orthogonal
       columns, so it won't be the same as the matrix input here.
       This "feature" is because the NIFTI 'qform' transform is
       deliberately not fully general -- it is intended to model a volume
       with perpendicular axes.
     - If the 3 input matrix columns are not even linearly independent,
       you'll just have to take your luck, won't you?

   \see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h

   \see nifti_quatern_to_mat44, nifti_make_orthog_mat44,
       nifti_mat44_to_orientation
*//*-------------------------------------------------------------------------*/
void nifti_mat44_to_quatern( mat44 R ,
                             float *qb, float *qc, float *qd,
                             float *qx, float *qy, float *qz,
                             float *dx, float *dy, float *dz, float *qfac )
{
   double r11,r12,r13 , r21,r22,r23 , r31,r32,r33 ;
   double xd,yd,zd , a,b,c,d ;
   mat33 P,Q ;

   /* offset outputs are read write out of input matrix  */

   ASSIF(qx,R.m[0][3]) ; ASSIF(qy,R.m[1][3]) ; ASSIF(qz,R.m[2][3]) ;

   /* load 3x3 matrix into local variables */

   r11 = R.m[0][0] ; r12 = R.m[0][1] ; r13 = R.m[0][2] ;
   r21 = R.m[1][0] ; r22 = R.m[1][1] ; r23 = R.m[1][2] ;
   r31 = R.m[2][0] ; r32 = R.m[2][1] ; r33 = R.m[2][2] ;

   /* compute lengths of each column; these determine grid spacings  */

   xd = sqrt( r11*r11 + r21*r21 + r31*r31 ) ;
   yd = sqrt( r12*r12 + r22*r22 + r32*r32 ) ;
   zd = sqrt( r13*r13 + r23*r23 + r33*r33 ) ;

   /* if a column length is zero, patch the trouble */

   if( xd == 0.0l ){ r11 = 1.0l ; r21 = r31 = 0.0l ; xd = 1.0l ; }
   if( yd == 0.0l ){ r22 = 1.0l ; r12 = r32 = 0.0l ; yd = 1.0l ; }
   if( zd == 0.0l ){ r33 = 1.0l ; r13 = r23 = 0.0l ; zd = 1.0l ; }

   /* assign the output lengths */

   ASSIF(dx,(float)xd) ; ASSIF(dy,(float)yd) ; ASSIF(dz,(float)zd) ;

   /* normalize the columns */

   r11 /= xd ; r21 /= xd ; r31 /= xd ;
   r12 /= yd ; r22 /= yd ; r32 /= yd ;
   r13 /= zd ; r23 /= zd ; r33 /= zd ;

   /* At this point, the matrix has normal columns, but we have to allow
      for the fact that the hideous user may not have given us a matrix
      with orthogonal columns.

      So, now find the orthogonal matrix closest to the current matrix.

      One reason for using the polar decomposition to get this
      orthogonal matrix, rather than just directly orthogonalizing
      the columns, is so that inputting the inverse matrix to R
      will result in the inverse orthogonal matrix at this point.
      If we just orthogonalized the columns, this wouldn't necessarily hold. */

   Q.m[0][0] = (float)r11 ; Q.m[0][1] = (float)r12 ; Q.m[0][2] = (float)r13 ; /* load Q */
   Q.m[1][0] = (float)r21 ; Q.m[1][1] = (float)r22 ; Q.m[1][2] = (float)r23 ;
   Q.m[2][0] = (float)r31 ; Q.m[2][1] = (float)r32 ; Q.m[2][2] = (float)r33 ;

   P = nifti_mat33_polar(Q) ;  /* P is orthog matrix closest to Q */

   r11 = P.m[0][0] ; r12 = P.m[0][1] ; r13 = P.m[0][2] ; /* unload */
   r21 = P.m[1][0] ; r22 = P.m[1][1] ; r23 = P.m[1][2] ;
   r31 = P.m[2][0] ; r32 = P.m[2][1] ; r33 = P.m[2][2] ;

   /*                            [ r11 r12 r13 ]               */
   /* at this point, the matrix  [ r21 r22 r23 ] is orthogonal */
   /*                            [ r31 r32 r33 ]               */

   /* compute the determinant to determine if it is proper */

   zd = r11*r22*r33-r11*r32*r23-r21*r12*r33
       +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;  /* should be -1 or 1 */

   if( zd > 0 ){             /* proper */
     ASSIF(qfac,1.0f) ;
   } else {                  /* improper ==> flip 3rd column */
     ASSIF(qfac,-1.0f) ;
     r13 = -r13 ; r23 = -r23 ; r33 = -r33 ;
   }

   /* now, compute quaternion parameters */

   a = r11 + r22 + r33 + 1.0l ;

   if( a > 0.5l ){                /* simplest case */
     a = 0.5l * sqrt(a) ;
     b = 0.25l * (r32-r23) / a ;
     c = 0.25l * (r13-r31) / a ;
     d = 0.25l * (r21-r12) / a ;
   } else {                       /* trickier case */
     xd = 1.0 + r11 - (r22+r33) ;  /* 4*b*b */
     yd = 1.0 + r22 - (r11+r33) ;  /* 4*c*c */
     zd = 1.0 + r33 - (r11+r22) ;  /* 4*d*d */
     if( xd > 1.0 ){
       b = 0.5l * sqrt(xd) ;
       c = 0.25l* (r12+r21) / b ;
       d = 0.25l* (r13+r31) / b ;
       a = 0.25l* (r32-r23) / b ;
     } else if( yd > 1.0 ){
       c = 0.5l * sqrt(yd) ;
       b = 0.25l* (r12+r21) / c ;
       d = 0.25l* (r23+r32) / c ;
       a = 0.25l* (r13-r31) / c ;
     } else {
       d = 0.5l * sqrt(zd) ;
       b = 0.25l* (r13+r31) / d ;
       c = 0.25l* (r23+r32) / d ;
       a = 0.25l* (r21-r12) / d ;
     }
     if( a < 0.0l ){ b=-b ; c=-c ; d=-d; a=-a; }
   }

   ASSIF(qb,(float)b) ; ASSIF(qc,(float)c) ; ASSIF(qd,(float)d) ;
   return ;
}

/*---------------------------------------------------------------------------*/
/*! Compute the inverse of a bordered 4x4 matrix.

     <pre>
   - Some numerical code fragments were generated by Maple 8.
   - If a singular matrix is input, the output matrix will be all zero.
   - You can check for this by examining the [3][3] element, which will
     be 1.0 for the normal case and 0.0 for the bad case.

     The input matrix should have the form:
        [ r11 r12 r13 v1 ]
        [ r21 r22 r23 v2 ]
        [ r31 r32 r33 v3 ]
        [  0   0   0   1 ]
     </pre>
*//*-------------------------------------------------------------------------*/
mat44 nifti_mat44_inverse( mat44 R )
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33,v1,v2,v3 , deti ;
   mat44 Q ;
                                                       /*  INPUT MATRIX IS:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 v1 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 v2 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 v3 ] */
   v1  = R.m[0][3]; v2  = R.m[1][3]; v3  = R.m[2][3];  /* [  0   0   0   1 ] */

   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;

   if( deti != 0.0l ) deti = 1.0l / deti ;

   Q.m[0][0] = (float)( deti*( r22*r33-r32*r23) ) ;
   Q.m[0][1] = (float)( deti*(-r12*r33+r32*r13) ) ;
   Q.m[0][2] = (float)( deti*( r12*r23-r22*r13) ) ;
   Q.m[0][3] = (float)( deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3
                     -r22*v1*r33-r32*r13*v2+r32*v1*r23) ) ;

   Q.m[1][0] = (float)( deti*(-r21*r33+r31*r23) ) ;
   Q.m[1][1] = (float)( deti*( r11*r33-r31*r13) ) ;
   Q.m[1][2] = (float)( deti*(-r11*r23+r21*r13) ) ;
   Q.m[1][3] = (float)( deti*( r11*r23*v3-r11*v2*r33-r21*r13*v3
                     +r21*v1*r33+r31*r13*v2-r31*v1*r23) ) ;

   Q.m[2][0] = (float)( deti*( r21*r32-r31*r22) ) ;
   Q.m[2][1] = (float)( deti*(-r11*r32+r31*r12) ) ;
   Q.m[2][2] = (float)( deti*( r11*r22-r21*r12) ) ;
   Q.m[2][3] = (float)( deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3
                     -r21*r32*v1-r31*r12*v2+r31*r22*v1) ) ;

   Q.m[3][0] = Q.m[3][1] = Q.m[3][2] = 0.0l ;
   Q.m[3][3] = (deti == 0.0l) ? 0.0l : 1.0l ; /* failure flag if deti == 0 */

   return Q ;
}

/*---------------------------------------------------------------------------*/
/*! Input 9 floats and make an orthgonal mat44 out of them.

   Each row is normalized, then nifti_mat33_polar() is used to orthogonalize
   them.  If row #3 (r31,r32,r33) is input as zero, then it will be taken to
   be the cross product of rows #1 and #2.

   This function can be used to create a rotation matrix for transforming
   an oblique volume to anatomical coordinates.  For this application:
    - row #1 (r11,r12,r13) is the direction vector along the image i-axis
    - row #2 (r21,r22,r23) is the direction vector along the image j-axis
    - row #3 (r31,r32,r33) is the direction vector along the slice direction
      (if available; otherwise enter it as 0's)

   The first 2 rows can be taken from the DICOM attribute (0020,0037)
   "Image Orientation (Patient)".

   After forming the rotation matrix, the complete affine transformation from
   (i,j,k) grid indexes to (x,y,z) spatial coordinates can be computed by
   multiplying each column by the appropriate grid spacing:
    - column #1 (R.m[0][0],R.m[1][0],R.m[2][0]) by delta-x
    - column #2 (R.m[0][1],R.m[1][1],R.m[2][1]) by delta-y
    - column #3 (R.m[0][2],R.m[1][2],R.m[2][2]) by delta-z

   and by then placing the center (x,y,z) coordinates of voxel (0,0,0) into
   the column #4 (R.m[0][3],R.m[1][3],R.m[2][3]).

   \sa nifti_quatern_to_mat44, nifti_mat44_to_quatern,
       nifti_mat44_to_orientation
*//*-------------------------------------------------------------------------*/
mat44 nifti_make_orthog_mat44( float r11, float r12, float r13 ,
                               float r21, float r22, float r23 ,
                               float r31, float r32, float r33  )
{
   mat44 R ;
   mat33 Q , P ;
   double val ;

   R.m[3][0] = R.m[3][1] = R.m[3][2] = 0.0l ; R.m[3][3] = 1.0l ;

   Q.m[0][0] = r11 ; Q.m[0][1] = r12 ; Q.m[0][2] = r13 ; /* load Q */
   Q.m[1][0] = r21 ; Q.m[1][1] = r22 ; Q.m[1][2] = r23 ;
   Q.m[2][0] = r31 ; Q.m[2][1] = r32 ; Q.m[2][2] = r33 ;

   /* normalize row 1 */

   val = Q.m[0][0]*Q.m[0][0] + Q.m[0][1]*Q.m[0][1] + Q.m[0][2]*Q.m[0][2] ;
   if( val > 0.0l ){
     val = 1.0l / sqrt(val) ;
     Q.m[0][0] *= (float)val ; Q.m[0][1] *= (float)val ; Q.m[0][2] *= (float)val ;
   } else {
     Q.m[0][0] = 1.0l ; Q.m[0][1] = 0.0l ; Q.m[0][2] = 0.0l ;
   }

   /* normalize row 2 */

   val = Q.m[1][0]*Q.m[1][0] + Q.m[1][1]*Q.m[1][1] + Q.m[1][2]*Q.m[1][2] ;
   if( val > 0.0l ){
     val = 1.0l / sqrt(val) ;
     Q.m[1][0] *= (float)val ; Q.m[1][1] *= (float)val ; Q.m[1][2] *= (float)val ;
   } else {
     Q.m[1][0] = 0.0l ; Q.m[1][1] = 1.0l ; Q.m[1][2] = 0.0l ;
   }

   /* normalize row 3 */

   val = Q.m[2][0]*Q.m[2][0] + Q.m[2][1]*Q.m[2][1] + Q.m[2][2]*Q.m[2][2] ;
   if( val > 0.0l ){
     val = 1.0l / sqrt(val) ;
     Q.m[2][0] *= (float)val ; Q.m[2][1] *= (float)val ; Q.m[2][2] *= (float)val ;
   } else {
     Q.m[2][0] = Q.m[0][1]*Q.m[1][2] - Q.m[0][2]*Q.m[1][1] ;  /* cross */
     Q.m[2][1] = Q.m[0][2]*Q.m[1][0] - Q.m[0][0]*Q.m[1][2] ;  /* product */
     Q.m[2][2] = Q.m[0][0]*Q.m[1][1] - Q.m[0][1]*Q.m[1][0] ;
   }

   P = nifti_mat33_polar(Q) ;  /* P is orthog matrix closest to Q */

   R.m[0][0] = P.m[0][0] ; R.m[0][1] = P.m[0][1] ; R.m[0][2] = P.m[0][2] ;
   R.m[1][0] = P.m[1][0] ; R.m[1][1] = P.m[1][1] ; R.m[1][2] = P.m[1][2] ;
   R.m[2][0] = P.m[2][0] ; R.m[2][1] = P.m[2][1] ; R.m[2][2] = P.m[2][2] ;

   R.m[0][3] = R.m[1][3] = R.m[2][3] = 0.0f ; return R ;
}

/*----------------------------------------------------------------------*/
/*! compute the inverse of a 3x3 matrix
*//*--------------------------------------------------------------------*/
mat33 nifti_mat33_inverse( mat33 R )   /* inverse of 3x3 matrix */
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33 , deti ;
   mat33 Q ;
                                                       /*  INPUT MATRIX:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;

   if( deti != 0.0l ) deti = 1.0l / deti ;

   Q.m[0][0] = (float)( deti*( r22*r33-r32*r23) ) ;
   Q.m[0][1] = (float)( deti*(-r12*r33+r32*r13) ) ;
   Q.m[0][2] = (float)( deti*( r12*r23-r22*r13) ) ;

   Q.m[1][0] = (float)( deti*(-r21*r33+r31*r23) ) ;
   Q.m[1][1] = (float)( deti*( r11*r33-r31*r13) ) ;
   Q.m[1][2] = (float)( deti*(-r11*r23+r21*r13) ) ;

   Q.m[2][0] = (float)( deti*( r21*r32-r31*r22) ) ;
   Q.m[2][1] = (float)( deti*(-r11*r32+r31*r12) ) ;
   Q.m[2][2] = (float)( deti*( r11*r22-r21*r12) ) ;

   return Q ;
}

/*----------------------------------------------------------------------*/
/*! compute the determinant of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float nifti_mat33_determ( mat33 R )   /* determinant of 3x3 matrix */
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33 ;
                                                       /*  INPUT MATRIX:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

   return (float)(r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13) ;
}

/*----------------------------------------------------------------------*/
/*! compute the max row norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float nifti_mat33_rownorm( mat33 A )  /* max row norm of 3x3 matrix */
{
   float r1,r2,r3 ;

   r1 = (float)( fabs(A.m[0][0])+fabs(A.m[0][1])+fabs(A.m[0][2]) ) ;
   r2 = (float)( fabs(A.m[1][0])+fabs(A.m[1][1])+fabs(A.m[1][2]) ) ;
   r3 = (float)( fabs(A.m[2][0])+fabs(A.m[2][1])+fabs(A.m[2][2]) ) ;
   if( r1 < r2 ) r1 = r2 ;
   if( r1 < r3 ) r1 = r3 ;
   return r1 ;
}

/*----------------------------------------------------------------------*/
/*! compute the max column norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float nifti_mat33_colnorm( mat33 A )  /* max column norm of 3x3 matrix */
{
   float r1,r2,r3 ;

   r1 = (float)( fabs(A.m[0][0])+fabs(A.m[1][0])+fabs(A.m[2][0]) ) ;
   r2 = (float)( fabs(A.m[0][1])+fabs(A.m[1][1])+fabs(A.m[2][1]) ) ;
   r3 = (float)( fabs(A.m[0][2])+fabs(A.m[1][2])+fabs(A.m[2][2]) ) ;
   if( r1 < r2 ) r1 = r2 ;
   if( r1 < r3 ) r1 = r3 ;
   return r1 ;
}

/*----------------------------------------------------------------------*/
/*! multiply 2 3x3 matrices
*//*--------------------------------------------------------------------*/
mat33 nifti_mat33_mul( mat33 A , mat33 B )  /* multiply 2 3x3 matrices */
{
   mat33 C ; int i,j ;
   for( i=0 ; i < 3 ; i++ )
    for( j=0 ; j < 3 ; j++ )
      C.m[i][j] =  A.m[i][0] * B.m[0][j]
                 + A.m[i][1] * B.m[1][j]
                 + A.m[i][2] * B.m[2][j] ;
   return C ;
}

/*---------------------------------------------------------------------------*/
/*! polar decomposition of a 3x3 matrix

   This finds the closest orthogonal matrix to input A
   (in both the Frobenius and L2 norms).

   Algorithm is that from NJ Higham, SIAM J Sci Stat Comput, 7:1160-1174.
*//*-------------------------------------------------------------------------*/
mat33 nifti_mat33_polar( mat33 A )
{
   mat33 X , Y , Z ;
   float alp,bet,gam,gmi , dif=1.0f ;
   int k=0 ;

   X = A ;

   /* force matrix to be nonsingular */

   gam = nifti_mat33_determ(X) ;
   while( gam == 0.0 ){        /* perturb matrix */
     gam = (float)( 0.00001 * ( 0.001 + nifti_mat33_rownorm(X) ) ) ;
     X.m[0][0] += gam ; X.m[1][1] += gam ; X.m[2][2] += gam ;
     gam = nifti_mat33_determ(X) ;
   }

   while(1){
     Y = nifti_mat33_inverse(X) ;
     if( dif > 0.3 ){     /* far from convergence */
       alp = (float)( sqrt( nifti_mat33_rownorm(X) * nifti_mat33_colnorm(X) ) ) ;
       bet = (float)( sqrt( nifti_mat33_rownorm(Y) * nifti_mat33_colnorm(Y) ) ) ;
       gam = (float)( sqrt( bet / alp ) ) ;
       gmi = (float)( 1.0 / gam ) ;
     } else {
       gam = gmi = 1.0f ;  /* close to convergence */
     }
     Z.m[0][0] = (float)( 0.5 * ( gam*X.m[0][0] + gmi*Y.m[0][0] ) ) ;
     Z.m[0][1] = (float)( 0.5 * ( gam*X.m[0][1] + gmi*Y.m[1][0] ) ) ;
     Z.m[0][2] = (float)( 0.5 * ( gam*X.m[0][2] + gmi*Y.m[2][0] ) ) ;
     Z.m[1][0] = (float)( 0.5 * ( gam*X.m[1][0] + gmi*Y.m[0][1] ) ) ;
     Z.m[1][1] = (float)( 0.5 * ( gam*X.m[1][1] + gmi*Y.m[1][1] ) ) ;
     Z.m[1][2] = (float)( 0.5 * ( gam*X.m[1][2] + gmi*Y.m[2][1] ) ) ;
     Z.m[2][0] = (float)( 0.5 * ( gam*X.m[2][0] + gmi*Y.m[0][2] ) ) ;
     Z.m[2][1] = (float)( 0.5 * ( gam*X.m[2][1] + gmi*Y.m[1][2] ) ) ;
     Z.m[2][2] = (float)( 0.5 * ( gam*X.m[2][2] + gmi*Y.m[2][2] ) ) ;

     dif = (float)( fabs(Z.m[0][0]-X.m[0][0])+fabs(Z.m[0][1]-X.m[0][1])
          +fabs(Z.m[0][2]-X.m[0][2])+fabs(Z.m[1][0]-X.m[1][0])
          +fabs(Z.m[1][1]-X.m[1][1])+fabs(Z.m[1][2]-X.m[1][2])
          +fabs(Z.m[2][0]-X.m[2][0])+fabs(Z.m[2][1]-X.m[2][1])
          +fabs(Z.m[2][2]-X.m[2][2])                          );

     k = k+1 ;
     if( k > 100 || dif < 3.e-6 ) break ;  /* convergence or exhaustion */
     X = Z ;
   }

   return Z ;
}

