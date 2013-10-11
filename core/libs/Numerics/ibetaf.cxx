/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010, 2013 SRI International
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
/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
      pseudocode.

See subroutines comments for additional copyrights.

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
*************************************************************************/


#include "ibetaf.h"

namespace 
alglib
{

ap::real_value_type incompletebetafe(ap::real_value_type a,
     ap::real_value_type b,
     ap::real_value_type x,
     ap::real_value_type big,
     ap::real_value_type biginv);
ap::real_value_type incompletebetafe2(ap::real_value_type a,
     ap::real_value_type b,
     ap::real_value_type x,
     ap::real_value_type big,
     ap::real_value_type biginv);
ap::real_value_type incompletebetaps(ap::real_value_type a, ap::real_value_type b, ap::real_value_type x, ap::real_value_type maxgam);

/*************************************************************************
Incomplete beta integral

Returns incomplete beta integral of the arguments, evaluated
from zero to x.  The function is defined as

                 x
    -            -
   | (a+b)      | |  a-1     b-1
 -----------    |   t   (1-t)   dt.
  -     -     | |
 | (a) | (b)   -
                0

The domain of definition is 0 <= x <= 1.  In this
implementation a and b are restricted to positive values.
The integral from x to 1 may be obtained by the symmetry
relation

   1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).

The integral is evaluated by a continued fraction expansion
or, when b*x is small, by a power series.

ACCURACY:

Tested at uniformly distributed random points (a,b,x) with a and b
in "domain" and x between 0 and 1.
                                       Relative error
arithmetic   domain     # trials      peak         rms
   IEEE      0,5         10000       6.9e-15     4.5e-16
   IEEE      0,85       250000       2.2e-13     1.7e-14
   IEEE      0,1000      30000       5.3e-12     6.3e-13
   IEEE      0,10000    250000       9.3e-11     7.1e-12
   IEEE      0,100000    10000       8.7e-10     4.8e-11
Outputs smaller than the IEEE gradual underflow threshold
were excluded from these statistics.

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
ap::real_value_type incompletebeta(ap::real_value_type a, ap::real_value_type b, ap::real_value_type x)
{
    ap::real_value_type result;
    ap::real_value_type t;
    ap::real_value_type xc;
    ap::real_value_type w;
    ap::real_value_type y;
    int flag;
    ap::real_value_type sg;
    ap::real_value_type big;
    ap::real_value_type biginv;
    ap::real_value_type maxgam;
    ap::real_value_type minlog;
    ap::real_value_type maxlog;

    big = 4.503599627370496e15;
    biginv = 2.22044604925031308085e-16;
    maxgam = 171.624376956302725;
    minlog = log(ap::minrealnumber);
    maxlog = log(ap::maxrealnumber);

#ifndef NO_AP_ASSERT
    ap::ap_error::make_assertion(a>0&&b>0, "Domain error in IncompleteBeta");
    ap::ap_error::make_assertion(x>=0&&x<=1, "Domain error in IncompleteBeta");
#endif

    if( x==0 )
    {
        result = 0;
        return result;
    }
    if( x==1 )
    {
        result = 1;
        return result;
    }
    flag = 0;
    if( b*x<=1.0&&x<=0.95 )
    {
        result = incompletebetaps(a, b, x, maxgam);
        return result;
    }
    w = 1.0-x;
    if( x>a/(a+b) )
    {
        flag = 1;
        t = a;
        a = b;
        b = t;
        xc = x;
        x = w;
    }
    else
    {
        xc = w;
    }
    if( flag==1&&b*x<=1.0&&x<=0.95 )
    {
        t = incompletebetaps(a, b, x, maxgam);
        if( t<=ap::machineepsilon )
        {
            result = 1.0-ap::machineepsilon;
        }
        else
        {
            result = 1.0-t;
        }
        return result;
    }
    y = x*(a+b-2.0)-(a-1.0);
    if( y<0.0 )
    {
        w = incompletebetafe(a, b, x, big, biginv);
    }
    else
    {
        w = incompletebetafe2(a, b, x, big, biginv)/xc;
    }
    y = a*log(x);
    t = b*log(xc);
    if( a+b<maxgam&&fabs(y)<maxlog&&fabs(t)<maxlog )
    {
        t = pow(xc, b);
        t = t*pow(x, a);
        t = t/a;
        t = t*w;
        t = t*(gamma(a+b)/(gamma(a)*gamma(b)));
        if( flag==1 )
        {
            if( t<=ap::machineepsilon )
            {
                result = 1.0-ap::machineepsilon;
            }
            else
            {
                result = 1.0-t;
            }
        }
        else
        {
            result = t;
        }
        return result;
    }
    y = y+t+lngamma(a+b, sg)-lngamma(a, sg)-lngamma(b, sg);
    y = y+log(w/a);
    if( y<minlog )
    {
        t = 0.0;
    }
    else
    {
        t = exp(y);
    }
    if( flag==1 )
    {
        if( t<=ap::machineepsilon )
        {
            t = 1.0-ap::machineepsilon;
        }
        else
        {
            t = 1.0-t;
        }
    }
    result = t;
    return result;
}


/*************************************************************************
Continued fraction expansion #1 for incomplete beta integral

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
ap::real_value_type incompletebetafe(ap::real_value_type a,
     ap::real_value_type b,
     ap::real_value_type x,
     ap::real_value_type big,
     ap::real_value_type biginv)
{
    ap::real_value_type result;
    ap::real_value_type xk;
    ap::real_value_type pk;
    ap::real_value_type pkm1;
    ap::real_value_type pkm2;
    ap::real_value_type qk;
    ap::real_value_type qkm1;
    ap::real_value_type qkm2;
    ap::real_value_type k1;
    ap::real_value_type k2;
    ap::real_value_type k3;
    ap::real_value_type k4;
    ap::real_value_type k5;
    ap::real_value_type k6;
    ap::real_value_type k7;
    ap::real_value_type k8;
    ap::real_value_type r;
    ap::real_value_type t;
    ap::real_value_type ans;
    ap::real_value_type thresh;
    int n;

    k1 = a;
    k2 = a+b;
    k3 = a;
    k4 = a+1.0;
    k5 = 1.0;
    k6 = b-1.0;
    k7 = k4;
    k8 = a+2.0;
    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0*ap::machineepsilon;
    do
    {
        xk = -x*k1*k2/(k3*k4);
        pk = pkm1+pkm2*xk;
        qk = qkm1+qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        xk = x*k5*k6/(k7*k8);
        pk = pkm1+pkm2*xk;
        qk = qkm1+qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( qk!=0 )
        {
            r = pk/qk;
        }
        if( r!=0 )
        {
            t = fabs((ans-r)/r);
            ans = r;
        }
        else
        {
            t = 1.0;
        }
        if( t<thresh )
        {
            break;
        }
        k1 = k1+1.0;
        k2 = k2+1.0;
        k3 = k3+2.0;
        k4 = k4+2.0;
        k5 = k5+1.0;
        k6 = k6-1.0;
        k7 = k7+2.0;
        k8 = k8+2.0;
        if( fabs(qk)+fabs(pk)>big )
        {
            pkm2 = pkm2*biginv;
            pkm1 = pkm1*biginv;
            qkm2 = qkm2*biginv;
            qkm1 = qkm1*biginv;
        }
        if( fabs(qk)<biginv||fabs(pk)<biginv )
        {
            pkm2 = pkm2*big;
            pkm1 = pkm1*big;
            qkm2 = qkm2*big;
            qkm1 = qkm1*big;
        }
        n = n+1;
    }
    while(n!=300);
    result = ans;
    return result;
}


/*************************************************************************
Continued fraction expansion #2
for incomplete beta integral

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
ap::real_value_type incompletebetafe2(ap::real_value_type a,
     ap::real_value_type b,
     ap::real_value_type x,
     ap::real_value_type big,
     ap::real_value_type biginv)
{
    ap::real_value_type result;
    ap::real_value_type xk;
    ap::real_value_type pk;
    ap::real_value_type pkm1;
    ap::real_value_type pkm2;
    ap::real_value_type qk;
    ap::real_value_type qkm1;
    ap::real_value_type qkm2;
    ap::real_value_type k1;
    ap::real_value_type k2;
    ap::real_value_type k3;
    ap::real_value_type k4;
    ap::real_value_type k5;
    ap::real_value_type k6;
    ap::real_value_type k7;
    ap::real_value_type k8;
    ap::real_value_type r;
    ap::real_value_type t;
    ap::real_value_type ans;
    ap::real_value_type z;
    ap::real_value_type thresh;
    int n;

    k1 = a;
    k2 = b-1.0;
    k3 = a;
    k4 = a+1.0;
    k5 = 1.0;
    k6 = a+b;
    k7 = a+1.0;
    k8 = a+2.0;
    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    z = x/(1.0-x);
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0*ap::machineepsilon;
    do
    {
        xk = -z*k1*k2/(k3*k4);
        pk = pkm1+pkm2*xk;
        qk = qkm1+qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        xk = z*k5*k6/(k7*k8);
        pk = pkm1+pkm2*xk;
        qk = qkm1+qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( qk!=0 )
        {
            r = pk/qk;
        }
        if( r!=0 )
        {
            t = fabs((ans-r)/r);
            ans = r;
        }
        else
        {
            t = 1.0;
        }
        if( t<thresh )
        {
            break;
        }
        k1 = k1+1.0;
        k2 = k2-1.0;
        k3 = k3+2.0;
        k4 = k4+2.0;
        k5 = k5+1.0;
        k6 = k6+1.0;
        k7 = k7+2.0;
        k8 = k8+2.0;
        if( fabs(qk)+fabs(pk)>big )
        {
            pkm2 = pkm2*biginv;
            pkm1 = pkm1*biginv;
            qkm2 = qkm2*biginv;
            qkm1 = qkm1*biginv;
        }
        if( fabs(qk)<biginv||fabs(pk)<biginv )
        {
            pkm2 = pkm2*big;
            pkm1 = pkm1*big;
            qkm2 = qkm2*big;
            qkm1 = qkm1*big;
        }
        n = n+1;
    }
    while(n!=300);
    result = ans;
    return result;
}


/*************************************************************************
Power series for incomplete beta integral.
Use when b*x is small and x not too close to 1.

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
ap::real_value_type incompletebetaps(ap::real_value_type a, ap::real_value_type b, ap::real_value_type x, ap::real_value_type maxgam)
{
    ap::real_value_type result;
    ap::real_value_type s;
    ap::real_value_type t;
    ap::real_value_type u;
    ap::real_value_type v;
    ap::real_value_type n;
    ap::real_value_type t1;
    ap::real_value_type z;
    ap::real_value_type ai;
    ap::real_value_type sg;

    ai = 1.0/a;
    u = (1.0-b)*x;
    v = u/(a+1.0);
    t1 = v;
    t = u;
    n = 2.0;
    s = 0.0;
    z = ap::machineepsilon*ai;
    while(fabs(v)>z)
    {
        u = (n-b)*x/n;
        t = t*u;
        v = t/(a+n);
        s = s+v;
        n = n+1.0;
    }
    s = s+t1;
    s = s+ai;
    u = a*log(x);
    if( a+b<maxgam&&fabs(u)<log(ap::maxrealnumber) )
    {
        t = gamma(a+b)/(gamma(a)*gamma(b));
        s = s*t*pow(x, a);
    }
    else
    {
        t = lngamma(a+b, sg)-lngamma(a, sg)-lngamma(b, sg)+u+log(s);
        if( t<log(ap::minrealnumber) )
        {
            s = 0.0;
        }
        else
        {
            s = exp(t);
        }
    }
    result = s;
    return result;
}

} // namespace alglib
