#ifndef cmtk_mxml_mangle_h
#define cmtk_mxml_mangle_h

/*

This header file mangles all symbols exported from the mxml library.
It is included in all files while building the mxml library. Due to
namespace pollution, no mxml headers should be included in .h files in
CMTK.

The following command was used to obtain the symbol list:

nm libcmtkMxml.so |grep " [TRD] "

This is the way to recreate the whole list:

nm libcmtkMxml.so |grep " [TRD] " | awk '{ print "#define "$3" cmtk_"$3 }'

REMOVE the "_init" and "_fini" entries.

*/

#define cmtk__mxml_entity_cb cmtk_cmtk__mxml_entity_cb
#define cmtk__mxml_global cmtk_cmtk__mxml_global
#define cmtk__mxml_strdupf cmtk_cmtk__mxml_strdupf
#define cmtk__mxml_vstrdupf cmtk_cmtk__mxml_vstrdupf
#define cmtk_mxmlAdd cmtk_cmtk_mxmlAdd
#define cmtk_mxmlDelete cmtk_cmtk_mxmlDelete
#define cmtk_mxmlElementDeleteAttr cmtk_cmtk_mxmlElementDeleteAttr
#define cmtk_mxmlElementGetAttr cmtk_cmtk_mxmlElementGetAttr
#define cmtk_mxmlElementSetAttr cmtk_cmtk_mxmlElementSetAttr
#define cmtk_mxmlElementSetAttrf cmtk_cmtk_mxmlElementSetAttrf
#define cmtk_mxmlEntityAddCallback cmtk_cmtk_mxmlEntityAddCallback
#define cmtk_mxmlEntityGetName cmtk_cmtk_mxmlEntityGetName
#define cmtk_mxmlEntityGetValue cmtk_cmtk_mxmlEntityGetValue
#define cmtk_mxmlEntityRemoveCallback cmtk_cmtk_mxmlEntityRemoveCallback
#define cmtk_mxmlFindElement cmtk_cmtk_mxmlFindElement
#define cmtk_mxmlIndexDelete cmtk_cmtk_mxmlIndexDelete
#define cmtk_mxmlIndexEnum cmtk_cmtk_mxmlIndexEnum
#define cmtk_mxmlIndexFind cmtk_cmtk_mxmlIndexFind
#define cmtk_mxmlIndexNew cmtk_cmtk_mxmlIndexNew
#define cmtk_mxmlIndexReset cmtk_cmtk_mxmlIndexReset
#define cmtk_mxmlLoadFd cmtk_cmtk_mxmlLoadFd
#define cmtk_mxmlLoadFile cmtk_cmtk_mxmlLoadFile
#define cmtk_mxmlLoadString cmtk_cmtk_mxmlLoadString
#define cmtk_mxmlNewCDATA cmtk_cmtk_mxmlNewCDATA
#define cmtk_mxmlNewCustom cmtk_cmtk_mxmlNewCustom
#define cmtk_mxmlNewElement cmtk_cmtk_mxmlNewElement
#define cmtk_mxmlNewInteger cmtk_cmtk_mxmlNewInteger
#define cmtk_mxmlNewOpaque cmtk_cmtk_mxmlNewOpaque
#define cmtk_mxmlNewReal cmtk_cmtk_mxmlNewReal
#define cmtk_mxmlNewText cmtk_cmtk_mxmlNewText
#define cmtk_mxmlNewTextf cmtk_cmtk_mxmlNewTextf
#define cmtk_mxmlNewXML cmtk_cmtk_mxmlNewXML
#define cmtk_mxmlRelease cmtk_cmtk_mxmlRelease
#define cmtk_mxmlRemove cmtk_cmtk_mxmlRemove
#define cmtk_mxmlRetain cmtk_cmtk_mxmlRetain
#define cmtk_mxmlSAXLoadFd cmtk_cmtk_mxmlSAXLoadFd
#define cmtk_mxmlSAXLoadFile cmtk_cmtk_mxmlSAXLoadFile
#define cmtk_mxmlSAXLoadString cmtk_cmtk_mxmlSAXLoadString
#define cmtk_mxmlSaveAllocString cmtk_cmtk_mxmlSaveAllocString
#define cmtk_mxmlSaveFd cmtk_cmtk_mxmlSaveFd
#define cmtk_mxmlSaveFile cmtk_cmtk_mxmlSaveFile
#define cmtk_mxmlSaveString cmtk_cmtk_mxmlSaveString
#define cmtk_mxmlSetCDATA cmtk_cmtk_mxmlSetCDATA
#define cmtk_mxmlSetCustom cmtk_cmtk_mxmlSetCustom
#define cmtk_mxmlSetCustomHandlers cmtk_cmtk_mxmlSetCustomHandlers
#define cmtk_mxmlSetElement cmtk_cmtk_mxmlSetElement
#define cmtk_mxmlSetErrorCallback cmtk_cmtk_mxmlSetErrorCallback
#define cmtk_mxmlSetInteger cmtk_cmtk_mxmlSetInteger
#define cmtk_mxmlSetOpaque cmtk_cmtk_mxmlSetOpaque
#define cmtk_mxmlSetReal cmtk_cmtk_mxmlSetReal
#define cmtk_mxmlSetText cmtk_cmtk_mxmlSetText
#define cmtk_mxmlSetTextf cmtk_cmtk_mxmlSetTextf
#define cmtk_mxmlSetWrapMargin cmtk_cmtk_mxmlSetWrapMargin
#define cmtk_mxmlWalkNext cmtk_cmtk_mxmlWalkNext
#define cmtk_mxmlWalkPrev cmtk_cmtk_mxmlWalkPrev
#define cmtk_mxml_error cmtk_cmtk_mxml_error
#define cmtk_mxml_ignore_cb cmtk_cmtk_mxml_ignore_cb
#define cmtk_mxml_integer_cb cmtk_cmtk_mxml_integer_cb
#define cmtk_mxml_opaque_cb cmtk_cmtk_mxml_opaque_cb
#define cmtk_mxml_real_cb cmtk_cmtk_mxml_real_cb
#define mxmlFindPath cmtk_mxmlFindPath
#define mxmlGetCDATA cmtk_mxmlGetCDATA
#define mxmlGetCustom cmtk_mxmlGetCustom
#define mxmlGetElement cmtk_mxmlGetElement
#define mxmlGetFirstChild cmtk_mxmlGetFirstChild
#define mxmlGetInteger cmtk_mxmlGetInteger
#define mxmlGetLastChild cmtk_mxmlGetLastChild
#define mxmlGetNextSibling cmtk_mxmlGetNextSibling
#define mxmlGetOpaque cmtk_mxmlGetOpaque
#define mxmlGetParent cmtk_mxmlGetParent
#define mxmlGetPrevSibling cmtk_mxmlGetPrevSibling
#define mxmlGetReal cmtk_mxmlGetReal
#define mxmlGetRefCount cmtk_mxmlGetRefCount
#define mxmlGetText cmtk_mxmlGetText
#define mxmlGetType cmtk_mxmlGetType
#define mxmlGetUserData cmtk_mxmlGetUserData
#define mxmlIndexGetCount cmtk_mxmlIndexGetCount
#define mxmlSetUserData cmtk_mxmlSetUserData

#endif
