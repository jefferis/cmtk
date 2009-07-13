#!/bin/sh

DCMTK_ROOT=$1

files="COPYRIGHT"

CMTK_ROOT="/Users/mphasak/Development/cmtk/trunk/core"
CMTK_DCMTK_ROOT="${CMTK_ROOT}/Utilities/dcmtk"
CMTK_DCMTK_INCLUDE="${CMTK_ROOT}/Utilities/dcmtk/dcmtk"

CONFIG_DIR="config"
CONFIG_MODULES="cfwin32 osconfig"

DCMDATA_DIR="dcmdata"
DCMDATA_MODULES="cmdlnarg dcbytstr dcchrstr dccodec dcdatset dcdebug dcdefine dcdeftag dcdicdir dcdicent dcdict dcdirrec dcelem dcerror dcfilefo dchashdi dcistrma dcistrmf dcitem dclist dcmetinf dcobject dcofsetl dcostrma dcostrmf dcostrmz dcovlay dcpcache dcpixel dcpixseq dcpxitem dcsequen dcstack dcswap dctag dctagkey dctk dctypes dcuid dcvm dcvr dcvrae dcvras dcvrat dcvrcs dcvrda dcvrds dcvrdt dcvrfd dcvrfl dcvris dcvrlo dcvrlt dcvrobow dcvrof dcvrpn dcvrpobw dcvrsh dcvrsl dcvrss dcvrst dcvrtm dcvrui dcvrul dcvrulup dcvrus dcvrut dcxfer"
DCMIMGLE_DIR="dcmimgle"
DCMIMGLE_MODULES="didocu diobjcou diutils"

DCMJPEG_DIR="dcmjpeg"
DCMJPEG_MODULES="djdecode djutils djdecbas djdecext djdecsps djdecpro djdecsv1 djdeclol djcparam djcodecd"

OFSTD_DIR="ofstd"
OFSTD_MODULES="ofalgo ofcast ofcond ofconsol ofglobal oflist ofstack ofstdinc ofstream ofstring ofthread oftypes ofdate ofdatime oftime"

copy_src_tgt()
{
  local module_dir=$1
  local module_list=$*

  src_h_path=${DCMTK_ROOT}/${module_dir}/include/dcmtk/${module_dir}
  tgt_h_path=${CMTK_DCMTK_INCLUDE}/${module_dir}
  mkdir -p ${tgt_h_path}
  
  src_c_path=${DCMTK_ROOT}/${module_dir}/libsrc
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
cp ${DCMTK_ROOT}/dcmdata/libsrc/dicom.dic ${CMTK_DCMTK_ROOT}


