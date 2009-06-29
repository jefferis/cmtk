#!/bin/sh

DCMTK_ROOT=$1

files="COPYRIGHT"

DCM_SOURCE_ROOT="/Users/mphasak/Downloads/dcmtk-3.5.4"
CMTK_ROOT="/Users/mphasak/Development/cmtk/trunk/core"
CMTK_DCMTK_ROOT="${CMTK_ROOT}/Utilities/dcmtk"
CMTK_DCMTK_INCLUDE="${CMTK_ROOT}/Utilities/dcmtk/dcmtk"

CONFIG_DIR="config"
CONFIG_MODULES="cfunix cfwin32 osconfig"

DCMDATA_DIR="dcmdata"
DCMDATA_MODULES="dcbytstr dcdeftag dcdict dcerror dchashdi dclist dcpcache dcstack dctagkey dctypes dcvrui dcdatset dcdicent dcelem dcfilefo dcitem dcobject dcsequen dctag dctk dcvr dcxfer"

DCMIMGLE_DIR="dcmimgle"
DCMIMGLE_MODULES="didocu diobjcou diutils"

DCMJPEG_DIR="dcmjpeg"
DCMJPEG_MODULES="djdecode djutils"

OFSTD_DIR="ofstd"
OFSTD_MODULES="ofalgo ofcast ofcond ofconsol ofglobal oflist ofstack ofstdinc ofstream ofstring ofthread oftypes"

copy_src_tgt()
{
  local module_dir=$1
  local module_list=$*

  src_h_path=${DCM_SOURCE_ROOT}/${module_dir}/include/dcmtk/${module_dir}
  tgt_h_path=${CMTK_DCMTK_INCLUDE}/${module_dir}
  mkdir -p ${tgt_h_path}
  
  src_c_path=${DCM_SOURCE_ROOT}/${module_dir}/libsrc
  #tgt_c_path=${CMTK_DCMTK_ROOT}/src/${module_dir}
  tgt_c_path=${CMTK_DCMTK_ROOT}
  mkdir -p ${tgt_c_path}

  for module in ${module_list}; do
    cp ${src_h_path}/${module}.h ${tgt_h_path}
    cp ${src_c_path}/${module}.cc ${tgt_c_path}
  done

}

copy_src_tgt ${CONFIG_DIR} ${CONFIG_MODULES};
copy_src_tgt ${DCMDATA_DIR} ${DCMDATA_MODULES};
copy_src_tgt ${DCMIMGLE_DIR} ${DCMIMGLE_MODULES};
copy_src_tgt ${DCMJPEG_DIR} ${DCMJPEG_MODULES};
copy_src_tgt ${OFSTD_DIR} ${OFSTD_MODULES};
cp ${DCM_SOURCE_ROOT}/config/include/dcmtk/config/cfunix.h.in ${CMTK_DCMTK_INCLUDE}/config
cp ${DCM_SOURCE_ROOT}/dcmdata/libsrc/dicom.dic ${CMTK_DCMTK_ROOT}


