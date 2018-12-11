/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2009, 2013 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//

*/
/********************************************************************
AP Library version 1.2
Copyright (c) 2003-2007, Sergey Bochkanov (ALGLIB project).
See www.alglib.net or alglib.sources.ru for details.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************/

#include "ap.h"

/********************************************************************
Optimized ABLAS interface
********************************************************************/
#ifdef AP_WIN32
#include <windows.h>
extern "C"
{
typedef ap::real_value_type  (*_ddot1)(const ap::real_value_type*, const ap::real_value_type*, long);
typedef void    (*_dmove1)(const ap::real_value_type*, const ap::real_value_type*, long);
typedef void    (*_dmoves1)(const ap::real_value_type*, const ap::real_value_type*, long, ap::real_value_type);
typedef void    (*_dmoveneg1)(const ap::real_value_type*, const ap::real_value_type*, long);
typedef void    (*_dadd1)(const ap::real_value_type*, const ap::real_value_type*, long);
typedef void    (*_dadds1)(const ap::real_value_type*, const ap::real_value_type*, long, ap::real_value_type);
typedef void    (*_dsub1)(const ap::real_value_type*, const ap::real_value_type*, long);
typedef void    (*_dmuls1)(const ap::real_value_type*, long, ap::real_value_type);
}
HINSTANCE ABLAS = LoadLibrary("ablas.dll");

static _ddot1     ddot1     = ABLAS==NULL ? NULL :     (_ddot1)  GetProcAddress(ABLAS, "ASMDotProduct1");
static _dmove1    dmove1    = ABLAS==NULL ? NULL :    (_dmove1)  GetProcAddress(ABLAS, "ASMMove1");
static _dmoves1   dmoves1   = ABLAS==NULL ? NULL :   (_dmoves1)  GetProcAddress(ABLAS, "ASMMoveS1");
static _dmoveneg1 dmoveneg1 = ABLAS==NULL ? NULL : (_dmoveneg1)  GetProcAddress(ABLAS, "ASMMoveNeg1");
static _dadd1     dadd1     = ABLAS==NULL ? NULL :     (_dadd1)  GetProcAddress(ABLAS, "ASMAdd1");
static _dadds1    dadds1    = ABLAS==NULL ? NULL :    (_dadds1)  GetProcAddress(ABLAS, "ASMAddS1");
static _dsub1     dsub1     = ABLAS==NULL ? NULL :     (_dsub1)  GetProcAddress(ABLAS, "ASMSub1");
static _dmuls1    dmuls1    = ABLAS==NULL ? NULL :     (_dmuls1) GetProcAddress(ABLAS, "ASMMulS1");
#endif

const ap::real_value_type ap::machineepsilon = 5E-16;
const ap::real_value_type ap::maxrealnumber  = 1E300;
const ap::real_value_type ap::minrealnumber  = 1E-300;

/********************************************************************
ap::complex operations
********************************************************************/
bool ap::operator==(const ap::complex& lhs, const ap::complex& rhs)
{ return lhs.x==rhs.x && lhs.y==rhs.y; }

bool ap::operator!=(const ap::complex& lhs, const ap::complex& rhs)
{ return lhs.x!=rhs.x || lhs.y!=rhs.y; }

const ap::complex ap::operator+(const ap::complex& lhs)
{ return lhs; }

const ap::complex ap::operator-(const ap::complex& lhs)
{ return ap::complex(-lhs.x, -lhs.y); }

const ap::complex ap::operator+(const ap::complex& lhs, const ap::complex& rhs)
{ ap::complex r = lhs; r += rhs; return r; }

const ap::complex ap::operator+(const ap::complex& lhs, const ap::real_value_type& rhs)
{ ap::complex r = lhs; r += rhs; return r; }

const ap::complex ap::operator+(const ap::real_value_type& lhs, const ap::complex& rhs)
{ ap::complex r = rhs; r += lhs; return r; }

const ap::complex ap::operator-(const ap::complex& lhs, const ap::complex& rhs)
{ ap::complex r = lhs; r -= rhs; return r; }

const ap::complex ap::operator-(const ap::complex& lhs, const ap::real_value_type& rhs)
{ ap::complex r = lhs; r -= rhs; return r; }

const ap::complex ap::operator-(const ap::real_value_type& lhs, const ap::complex& rhs)
{ ap::complex r = lhs; r -= rhs; return r; }

const ap::complex ap::operator*(const ap::complex& lhs, const ap::complex& rhs)
{ return ap::complex(lhs.x*rhs.x - lhs.y*rhs.y,  lhs.x*rhs.y + lhs.y*rhs.x); }

const ap::complex ap::operator*(const ap::complex& lhs, const ap::real_value_type& rhs)
{ return ap::complex(lhs.x*rhs,  lhs.y*rhs); }

const ap::complex ap::operator*(const ap::real_value_type& lhs, const ap::complex& rhs)
{ return ap::complex(lhs*rhs.x,  lhs*rhs.y); }

const ap::complex ap::operator/(const ap::complex& lhs, const ap::complex& rhs)
{
    ap::complex result;
    ap::real_value_type e;
    ap::real_value_type f;
    if( fabs(rhs.y)<fabs(rhs.x) )
    {
        e = rhs.y/rhs.x;
        f = rhs.x+rhs.y*e;
        result.x = (lhs.x+lhs.y*e)/f;
        result.y = (lhs.y-lhs.x*e)/f;
    }
    else
    {
        e = rhs.x/rhs.y;
        f = rhs.y+rhs.x*e;
        result.x = (lhs.y+lhs.x*e)/f;
        result.y = (-lhs.x+lhs.y*e)/f;
    }
    return result;
}

const ap::complex ap::operator/(const ap::real_value_type& lhs, const ap::complex& rhs)
{
    ap::complex result;
    ap::real_value_type e;
    ap::real_value_type f;
    if( fabs(rhs.y)<fabs(rhs.x) )
    {
        e = rhs.y/rhs.x;
        f = rhs.x+rhs.y*e;
        result.x = lhs/f;
        result.y = -lhs*e/f;
    }
    else
    {
        e = rhs.x/rhs.y;
        f = rhs.y+rhs.x*e;
        result.x = lhs*e/f;
        result.y = -lhs/f;
    }
    return result;
}

const ap::complex ap::operator/(const ap::complex& lhs, const ap::real_value_type& rhs)
{ return ap::complex(lhs.x/rhs, lhs.y/rhs); }

const ap::complex ap::conj(const ap::complex &z)
{ return ap::complex(z.x, -z.y); }

const ap::complex ap::csqr(const ap::complex &z)
{ return ap::complex(z.x*z.x-z.y*z.y, 2*z.x*z.y); }

/********************************************************************
BLAS functions
********************************************************************/
ap::real_value_type ap::vdotproduct(const ap::real_value_type *v1, const ap::real_value_type *v2, int N)
{
#ifdef AP_WIN32
    if( ddot1!=NULL )
        return ddot1(v1, v2, N);
#endif
    return ap::_vdotproduct<ap::real_value_type>(v1, v2, N);
}

ap::complex ap::vdotproduct(const ap::complex *v1, const ap::complex *v2, int N)
{
    return ap::_vdotproduct<ap::complex>(v1, v2, N);
}

void ap::vmove(double *vdst, const double* vsrc, int N)
{
#ifdef AP_WIN32
    if( dmove1!=NULL )
    {
        dmove1(vdst, vsrc, N);
        return;
    }
#endif
    ap::_vmove<double>(vdst, vsrc, N);
}

void ap::vmove(float *vdst, const float* vsrc, int N)
{
//#ifdef AP_WIN32
//    if( dmove1!=NULL )
//    {
//        dmove1(vdst, vsrc, N);
//        return;
//    }
//#endif
    ap::_vmove<float>(vdst, vsrc, N);
}

void ap::vmove(ap::complex *vdst, const ap::complex* vsrc, int N)
{
    ap::_vmove<ap::complex>(vdst, vsrc, N);
}

void ap::vmove(double *vdst, const double *vsrc, int N, double alpha)
{
#ifdef AP_WIN32
    if( dmoves1!=NULL )
    {
        dmoves1(vdst, vsrc, N, alpha);
        return;
    }
#endif
    ap::_vmove<double,double>(vdst, vsrc, N, alpha);
}

void ap::vmove(float *vdst, const float *vsrc, int N, float alpha)
{
//#ifdef AP_WIN32
//    if( dmoves1!=NULL )
//    {
//        dmoves1(vdst, vsrc, N, alpha);
//        return;
//    }
//#endif
    ap::_vmove<float,float>(vdst, vsrc, N, alpha);
}

void ap::vmove(ap::complex *vdst, const ap::complex *vsrc, int N, ap::real_value_type alpha)
{
    ap::_vmove<ap::complex,ap::real_value_type>(vdst, vsrc, N, alpha);
}

void ap::vmove(ap::complex *vdst, const ap::complex *vsrc, int N, ap::complex alpha)
{
    ap::_vmove<ap::complex,ap::complex>(vdst, vsrc, N, alpha);
}

void ap::vadd(ap::real_value_type *vdst, const ap::real_value_type *vsrc, int N)
{
#ifdef AP_WIN32
    if( dadd1!=NULL )
    {
        dadd1(vdst, vsrc, N);
        return;
    }
#endif
    ap::_vadd<ap::real_value_type>(vdst, vsrc, N);
}

void ap::vadd(ap::complex *vdst, const ap::complex *vsrc, int N)
{
    ap::_vadd<ap::complex>(vdst, vsrc, N);
}

void ap::vadd(ap::real_value_type *vdst, const ap::real_value_type *vsrc, int N, ap::real_value_type alpha)
{
#ifdef AP_WIN32
    if( dadds1!=NULL )
    {
        dadds1(vdst, vsrc, N, alpha);
        return;
    }
#endif
    ap::_vadd<ap::real_value_type,ap::real_value_type>(vdst, vsrc, N, alpha);
}

void ap::vadd(ap::complex *vdst, const ap::complex *vsrc, int N, ap::real_value_type alpha)
{
  ap::_vadd<ap::complex,ap::real_value_type>(vdst, vsrc, N, alpha);
}

void ap::vadd(ap::complex *vdst, const ap::complex *vsrc, int N, ap::complex alpha)
{
    ap::_vadd<ap::complex,ap::complex>(vdst, vsrc, N, alpha);
}

void ap::vsub(ap::real_value_type *vdst, const ap::real_value_type *vsrc, int N)
{
#ifdef AP_WIN32
    if( dsub1!=NULL )
    {
        dsub1(vdst, vsrc, N);
        return;
    }
#endif
    ap::_vsub<ap::real_value_type>(vdst, vsrc, N);
}

void ap::vsub(ap::complex *vdst, const ap::complex *vsrc, int N)
{
    ap::_vsub<ap::complex>(vdst, vsrc, N);
}

void ap::vsub(ap::real_value_type *vdst, const ap::real_value_type *vsrc, int N, ap::real_value_type alpha)
{
#ifdef AP_WIN32
    if( dadds1!=NULL )
    {
        dadds1(vdst, vsrc, N, -alpha);
        return;
    }
#endif
    ap::_vsub<ap::real_value_type,ap::real_value_type>(vdst, vsrc, N, alpha);
}

void ap::vsub(ap::complex *vdst, const ap::complex *vsrc, int N, ap::real_value_type alpha)
{
    ap::_vsub<ap::complex,ap::real_value_type>(vdst, vsrc, N, alpha);
}

void ap::vsub(ap::complex *vdst, const ap::complex *vsrc, int N, ap::complex alpha)
{
    ap::_vsub<ap::complex,ap::complex>(vdst, vsrc, N, alpha);
}

void ap::vmul(ap::real_value_type *vdst, int N, ap::real_value_type alpha)
{
#ifdef AP_WIN32
    if( dmuls1!=NULL )
    {
        dmuls1(vdst, N, alpha);
        return;
    }
#endif
    ap::_vmul<ap::real_value_type,ap::real_value_type>(vdst, N, alpha);
}

void ap::vmul(ap::complex *vdst, int N, ap::real_value_type alpha)
{
    ap::_vmul<ap::complex,ap::real_value_type>(vdst, N, alpha);
}

void ap::vmul(ap::complex *vdst, int N, ap::complex alpha)
{
    ap::_vmul<ap::complex,ap::complex>(vdst, N, alpha);
}

/********************************************************************
standard functions
********************************************************************/
int ap::sign(ap::real_value_type x)
{
    if( x>0 ) return  1;
    if( x<0 ) return -1;
    return 0;
}

int ap::round(ap::real_value_type x)
{ return int(floor(x+0.5)); }

int ap::trunc(ap::real_value_type x)
{ return int(x>0 ? floor(x) : ceil(x)); }

int ap::ifloor(ap::real_value_type x)
{ return int(floor(x)); }

ap::real_value_type ap::pi()
{ return 3.14159265358979323846; }

ap::real_value_type ap::sqr(ap::real_value_type x)
{ return x*x; }

int ap::maxint(int m1, int m2)
{
    return m1>m2 ? m1 : m2;
}

int ap::minint(int m1, int m2)
{
    return m1>m2 ? m2 : m1;
}

ap::real_value_type ap::maxreal(ap::real_value_type m1, ap::real_value_type m2)
{
    return m1>m2 ? m1 : m2;
}

ap::real_value_type ap::minreal(ap::real_value_type m1, ap::real_value_type m2)
{
    return m1>m2 ? m2 : m1;
}

/********************************************************************
Service routines:
********************************************************************/
void* ap::amalloc(size_t size, size_t alignment)
{
    if( alignment<=1 )
    {
        //
        // no alignment, just call malloc
        //
        void *block = malloc(sizeof(void*)+size);
        void **p = (void**)block;
        *p = block;
        return (void*)((char*)block+sizeof(void*));
    }
    else
    {
        //
        // align.
        //
        void *block = malloc(alignment-1+sizeof(void*)+size);
        char *result = (char*)block+sizeof(void*);
        //if( ((unsigned int)(result))%alignment!=0 )
        //    result += alignment - ((unsigned int)(result))%alignment;
        if( (result-(char*)0)%alignment!=0 )
            result += alignment - (result-(char*)0)%alignment;
        *((void**)(result-sizeof(void*))) = block;
        return result;
    }
}

void ap::afree(void *block)
{
    void *p = *((void**)((char*)block-sizeof(void*)));
    free(p);
}

int ap::vlen(int n1, int n2)
{
    return n2-n1+1;
}

