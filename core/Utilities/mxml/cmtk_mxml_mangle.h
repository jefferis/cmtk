#ifndef cmtk_mxml_mangle_h
#define cmtk_mxml_mangle_h

/*

This header file mangles all symbols exported from the mxml library.
It is included in all files while building the mxml library. Due to
namespace pollution, no mxml headers should be included in .h files in
VTK.

The following command was used to obtain the symbol list:

nm libcmtkmxml.so |grep " [TRD] "

This is the way to recreate the whole list:

nm libcmtkmxml.so |grep " [TRD] " | awk '{ print "#define "$3" cmtk_mxml_"$3 }'

REMOVE the "_init" and "_fini" entries.

*/

#define mxmlElementDeleteAttr cmtk_mxmlElementDeleteAttr
#define mxmlElementGetAttr cmtk_mxmlElementGetAttr
#define mxmlElementSetAttr cmtk_mxmlElementSetAttr
#define mxmlElementSetAttrf cmtk_mxmlElementSetAttrf
#define _mxml_entity_cb cmtk__mxml_entity_cb
#define mxmlEntityAddCallback cmtk_mxmlEntityAddCallback
#define mxmlEntityGetName cmtk_mxmlEntityGetName
#define mxmlEntityGetValue cmtk_mxmlEntityGetValue
#define mxmlEntityRemoveCallback cmtk_mxmlEntityRemoveCallback
#define mxmlIndexDelete cmtk_mxmlIndexDelete
#define mxmlIndexEnum cmtk_mxmlIndexEnum
#define mxmlIndexFind cmtk_mxmlIndexFind
#define mxmlIndexNew cmtk_mxmlIndexNew
#define mxmlIndexReset cmtk_mxmlIndexReset
#define _mxml_global cmtk__mxml_global
#define mxml_error cmtk_mxml_error
#define mxml_ignore_cb cmtk_mxml_ignore_cb
#define mxml_integer_cb cmtk_mxml_integer_cb
#define mxml_opaque_cb cmtk_mxml_opaque_cb
#define mxml_real_cb cmtk_mxml_real_cb
#define mxmlSetCDATA cmtk_mxmlSetCDATA
#define mxmlSetCustom cmtk_mxmlSetCustom
#define mxmlSetElement cmtk_mxmlSetElement
#define mxmlSetInteger cmtk_mxmlSetInteger
#define mxmlSetOpaque cmtk_mxmlSetOpaque
#define mxmlSetReal cmtk_mxmlSetReal
#define mxmlSetText cmtk_mxmlSetText
#define mxmlSetTextf cmtk_mxmlSetTextf
#define mxmlLoadFd cmtk_mxmlLoadFd
#define mxmlLoadFile cmtk_mxmlLoadFile
#define mxmlLoadString cmtk_mxmlLoadString
#define mxmlSAXLoadFd cmtk_mxmlSAXLoadFd
#define mxmlSAXLoadFile cmtk_mxmlSAXLoadFile
#define mxmlSAXLoadString cmtk_mxmlSAXLoadString
#define mxmlSaveAllocString cmtk_mxmlSaveAllocString
#define mxmlSaveFd cmtk_mxmlSaveFd
#define mxmlSaveFile cmtk_mxmlSaveFile
#define mxmlSaveString cmtk_mxmlSaveString
#define mxmlSetCustomHandlers cmtk_mxmlSetCustomHandlers
#define mxmlSetErrorCallback cmtk_mxmlSetErrorCallback
#define mxmlSetWrapMargin cmtk_mxmlSetWrapMargin
#define mxmlAdd cmtk_mxmlAdd
#define mxmlDelete cmtk_mxmlDelete
#define mxmlNewCDATA cmtk_mxmlNewCDATA
#define mxmlNewCustom cmtk_mxmlNewCustom
#define mxmlNewElement cmtk_mxmlNewElement
#define mxmlNewInteger cmtk_mxmlNewInteger
#define mxmlNewOpaque cmtk_mxmlNewOpaque
#define mxmlNewReal cmtk_mxmlNewReal
#define mxmlNewText cmtk_mxmlNewText
#define mxmlNewTextf cmtk_mxmlNewTextf
#define mxmlNewXML cmtk_mxmlNewXML
#define mxmlRelease cmtk_mxmlRelease
#define mxmlRemove cmtk_mxmlRemove
#define mxmlRetain cmtk_mxmlRetain
#define mxmlFindElement cmtk_mxmlFindElement
#define mxmlWalkNext cmtk_mxmlWalkNext
#define mxmlWalkPrev cmtk_mxmlWalkPrev
#define _mxml_strdupf cmtk__mxml_strdupf
#define _mxml_vstrdupf cmtk__mxml_vstrdupf

#endif
