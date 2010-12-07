#ifndef cmtk_zlib_mangle_h
#define cmtk_zlib_mangle_h

/*

This header file mangles all symbols exported from the zlib library.
It is included in all files while building the zlib library. Due to
namespace pollution, no zlib headers should be included in .h files in
VTK.

The following command was used to obtain the symbol list:

nm libcmtkzlib.so |grep " [TRD] "

This is the way to recreate the whole list:

nm libcmtkzlib.so |grep " [TRD] " | awk '{ print "#define "$3" cmtk_zlib_"$3 }'

REMOVE the "_init" and "_fini" entries.

*/

#define adler32 cmtk_zlib_adler32
#define adler32_combine cmtk_zlib_adler32_combine
#define compress cmtk_zlib_compress
#define compress2 cmtk_zlib_compress2
#define compressBound cmtk_zlib_compressBound
#define crc32 cmtk_zlib_crc32
#define crc32_combine cmtk_zlib_crc32_combine
#define get_crc_table cmtk_zlib_get_crc_table
#define deflate cmtk_zlib_deflate
#define deflateBound cmtk_zlib_deflateBound
#define deflateCopy cmtk_zlib_deflateCopy
#define deflateEnd cmtk_zlib_deflateEnd
#define deflateInit2_ cmtk_zlib_deflateInit2_
#define deflateInit_ cmtk_zlib_deflateInit_
#define deflateParams cmtk_zlib_deflateParams
#define deflatePrime cmtk_zlib_deflatePrime
#define deflateReset cmtk_zlib_deflateReset
#define deflateSetDictionary cmtk_zlib_deflateSetDictionary
#define deflateSetHeader cmtk_zlib_deflateSetHeader
#define deflateTune cmtk_zlib_deflateTune
#define deflate_copyright cmtk_zlib_deflate_copyright
#define gzclearerr cmtk_zlib_gzclearerr
#define gzclose cmtk_zlib_gzclose
#define gzdirect cmtk_zlib_gzdirect
#define gzdopen cmtk_zlib_gzdopen
#define gzeof cmtk_zlib_gzeof
#define gzerror cmtk_zlib_gzerror
#define gzflush cmtk_zlib_gzflush
#define gzgetc cmtk_zlib_gzgetc
#define gzgets cmtk_zlib_gzgets
#define gzopen cmtk_zlib_gzopen
#define gzprintf cmtk_zlib_gzprintf
#define gzputc cmtk_zlib_gzputc
#define gzputs cmtk_zlib_gzputs
#define gzread cmtk_zlib_gzread
#define gzrewind cmtk_zlib_gzrewind
#define gzseek cmtk_zlib_gzseek
#define gzsetparams cmtk_zlib_gzsetparams
#define gztell cmtk_zlib_gztell
#define gzungetc cmtk_zlib_gzungetc
#define gzwrite cmtk_zlib_gzwrite
#define inflate_fast cmtk_zlib_inflate_fast
#define inflate cmtk_zlib_inflate
#define inflateCopy cmtk_zlib_inflateCopy
#define inflateEnd cmtk_zlib_inflateEnd
#define inflateGetHeader cmtk_zlib_inflateGetHeader
#define inflateInit2_ cmtk_zlib_inflateInit2_
#define inflateInit_ cmtk_zlib_inflateInit_
#define inflatePrime cmtk_zlib_inflatePrime
#define inflateReset cmtk_zlib_inflateReset
#define inflateSetDictionary cmtk_zlib_inflateSetDictionary
#define inflateSync cmtk_zlib_inflateSync
#define inflateSyncPoint cmtk_zlib_inflateSyncPoint
#define inflate_copyright cmtk_zlib_inflate_copyright
#define inflate_table cmtk_zlib_inflate_table
#define _dist_code cmtk_zlib__dist_code
#define _length_code cmtk_zlib__length_code
#define _tr_align cmtk_zlib__tr_align
#define _tr_flush_block cmtk_zlib__tr_flush_block
#define _tr_init cmtk_zlib__tr_init
#define _tr_stored_block cmtk_zlib__tr_stored_block
#define _tr_tally cmtk_zlib__tr_tally
#define uncompress cmtk_zlib_uncompress
#define zError cmtk_zlib_zError
#define z_errmsg cmtk_zlib_z_errmsg
#define zcalloc cmtk_zlib_zcalloc
#define zcfree cmtk_zlib_zcfree
#define zlibCompileFlags cmtk_zlib_zlibCompileFlags
#define zlibVersion cmtk_zlib_zlibVersion

#endif
