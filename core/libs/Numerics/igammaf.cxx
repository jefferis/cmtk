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

#include "igammaf.h"

/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
ap::real_value_type incompletegamma(ap::real_value_type a, ap::real_value_type x)
{
    ap::real_value_type result;
    ap::real_value_type igammaepsilon;
    ap::real_value_type ans;
    ap::real_value_type ax;
    ap::real_value_type c;
    ap::real_value_type r;
    ap::real_value_type tmp;

    igammaepsilon = 0.000000000000001;
    if( x<=0||a<=0 )
    {
        result = 0;
        return result;
    }
    if( x>1&&x>a )
    {
        result = 1-incompletegammac(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax<-709.78271289338399 )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    r = a;
    c = 1;
    ans = 1;
    do
    {
        r = r+1;
        c = c*x/r;
        ans = ans+c;
    }
    while(c/ans>igammaepsilon);
    result = ans*ax/a;
    return result;
}


/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
ap::real_value_type incompletegammac(ap::real_value_type a, ap::real_value_type x)
{
    ap::real_value_type result;
    ap::real_value_type igammaepsilon;
    ap::real_value_type igammabignumber;
    ap::real_value_type igammabignumberinv;
    ap::real_value_type ans;
    ap::real_value_type ax;
    ap::real_value_type c;
    ap::real_value_type yc;
    ap::real_value_type r;
    ap::real_value_type t;
    ap::real_value_type y;
    ap::real_value_type z;
    ap::real_value_type pk;
    ap::real_value_type pkm1;
    ap::real_value_type pkm2;
    ap::real_value_type qk;
    ap::real_value_type qkm1;
    ap::real_value_type qkm2;
    ap::real_value_type tmp;

    igammaepsilon = 0.000000000000001;
    igammabignumber = 4503599627370496.0;
    igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
    if( x<=0||a<=0 )
    {
        result = 1;
        return result;
    }
    if( x<1||x<a )
    {
        result = 1-incompletegamma(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax<-709.78271289338399 )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    y = 1-a;
    z = x+y+1;
    c = 0;
    pkm2 = 1;
    qkm2 = x;
    pkm1 = x+1;
    qkm1 = z*x;
    ans = pkm1/qkm1;
    do
    {
        c = c+1;
        y = y+1;
        z = z+2;
        yc = y*c;
        pk = pkm1*z-pkm2*yc;
        qk = qkm1*z-qkm2*yc;
        if( qk!=0 )
        {
            r = pk/qk;
            t = fabs((ans-r)/r);
            ans = r;
        }
        else
        {
            t = 1;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( fabs(pk)>igammabignumber )
        {
            pkm2 = pkm2*igammabignumberinv;
            pkm1 = pkm1*igammabignumberinv;
            qkm2 = qkm2*igammabignumberinv;
            qkm1 = qkm1*igammabignumberinv;
        }
    }
    while(t>igammaepsilon);
    result = ans*ax;
    return result;
}





