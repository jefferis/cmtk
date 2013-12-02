/** \file nifti1_io.h
    \brief Data structures for using nifti1_io API.
           - Written by Bob Cox, SSCC NIMH
           - Revisions by Rick Reynolds, SSCC NIMH
 */
#ifndef _NIFTI_IO_HEADER_
#define _NIFTI_IO_HEADER_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*=================*/
#ifdef  __cplusplus
extern "C" {
#endif
/*=================*/

/*****===================================================================*****/
/*****         File nifti1_io.h == Declarations for nifti1_io.c          *****/
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

/*
   Modified by: Mark Jenkinson (FMRIB Centre, University of Oxford, UK)
   Date: July/August 2004

      Mainly adding low-level IO and changing things to allow gzipped files
      to be read and written
      Full backwards compatability should have been maintained

   Modified by: Rick Reynolds (SSCC/DIRP/NIMH, National Institutes of Health)
   Date: December 2004

      Modified and added many routines for I/O.
*/

#undef  ASSIF                                 /* assign v to *p, if possible */
#define ASSIF(p,v) if( (p)!=NULL ) *(p) = (v)

/********************** Some sample data structures **************************/

typedef struct {                   /** 4x4 matrix struct **/
  double m[4][4] ;
} mat44 ;

typedef struct {                   /** 3x3 matrix struct **/
  double m[3][3] ;
} mat33 ;


/*****************************************************************************/
/*--------------- Prototypes of functions defined in this file --------------*/

mat44 nifti_mat44_inverse( mat44 R ) ;

mat33 nifti_mat33_inverse( mat33 R ) ;
mat33 nifti_mat33_polar  ( mat33 A ) ;
double nifti_mat33_rownorm( mat33 A ) ;
double nifti_mat33_colnorm( mat33 A ) ;
double nifti_mat33_determ ( mat33 R ) ;
mat33 nifti_mat33_mul    ( mat33 A , mat33 B ) ;

void nifti_mat44_to_quatern( mat44 R ,
                             float *qb, float *qc, float *qd,
                             float *qx, float *qy, float *qz,
                             float *dx, float *dy, float *dz, float *qfac ) ;

mat44 nifti_quatern_to_mat44( const float qb, const float qc, const float qd,
                              const float qx, const float qy, const float qz,
                              const float dx, const float dy, const float dz, const float qfac );

mat44 nifti_make_orthog_mat44( const double r11, const double r12, const double r13 ,
                               const double r21, const double r22, const double r23 ,
                               const double r31, const double r32, const double r33  ) ;


/*=================*/
#ifdef  __cplusplus
}
#endif
/*=================*/

#endif /* _NIFTI_IO_HEADER_ */
