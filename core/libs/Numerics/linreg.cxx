/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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
//  $Revision: 2373 $
//
//  $LastChangedDate: 2010-09-29 11:17:31 -0700 (Wed, 29 Sep 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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

#include "linreg.h"

static const int lrvnum = 5;

static void lrinternal(const ap::real_2d_array& xy,
     const ap::real_1d_array& s,
     int npoints,
     int nvars,
     int& info,
     ap::linearmodel& lm,
     ap::lrreport& ar);

/*************************************************************************
Linear regression

Subroutine builds model:

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)

and model found in ALGLIB format, covariation matrix, training set  errors
(rms,  average,  average  relative)   and  leave-one-out  cross-validation
estimate of the generalization error. CV  estimate calculated  using  fast
algorithm with O(NPoints*NVars) complexity.

When  covariation  matrix  is  calculated  standard deviations of function
values are assumed to be equal to RMS error on the training set.

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuild(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     ap::linearmodel& lm,
     ap::lrreport& ar)
{
    ap::real_1d_array s;
    int i;
    ap::real_value_type sigma2;

    if( npoints<=nvars+1||nvars<1 )
    {
        info = -1;
        return;
    }
    s.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        s(i) = 1;
    }
    lrbuilds(xy, s, npoints, nvars, info, lm, ar);
    if( info<0 )
    {
        return;
    }
    sigma2 = ap::sqr(ar.rmserror)*npoints/(npoints-nvars-1);
    for(i = 0; i <= nvars; i++)
    {
        ap::vmul(&ar.c(i, 0), ap::vlen(0,nvars), sigma2);
    }
}


/*************************************************************************
Linear regression

Variant of LRBuild which uses vector of standatd deviations (errors in
function values).

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    S           -   standard deviations (errors in function values)
                    array[0..NPoints-1], S[i]>0.
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    * -2, if S[I]<=0
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuilds(const ap::real_2d_array& xy,
     const ap::real_1d_array& s,
     int npoints,
     int nvars,
     int& info,
     ap::linearmodel& lm,
     ap::lrreport& ar)
{
    ap::real_2d_array xyi;
    ap::real_1d_array x;
    ap::real_1d_array means;
    ap::real_1d_array sigmas;
    int i;
    int j;
    ap::real_value_type v;
    int offs;
    ap::real_value_type mean;
    ap::real_value_type variance;
    ap::real_value_type skewness;
    ap::real_value_type kurtosis;

    
    //
    // Test parameters
    //
    if( npoints<=nvars+1||nvars<1 )
    {
        info = -1;
        return;
    }
    
    //
    // Copy data, add one more column (constant term)
    //
    xyi.setbounds(0, npoints-1, 0, nvars+1);
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&xyi(i, 0), &xy(i, 0), ap::vlen(0,nvars-1));
        xyi(i,nvars) = 1;
        xyi(i,nvars+1) = xy(i,nvars);
    }
    
    //
    // Standartization
    //
    x.setbounds(0, npoints-1);
    means.setbounds(0, nvars-1);
    sigmas.setbounds(0, nvars-1);
    for(j = 0; j <= nvars-1; j++)
    {
        ap::vmove(x.getvector(0, npoints-1), xy.getcolumn(j, 0, npoints-1));
        calculatemoments(x, npoints, mean, variance, skewness, kurtosis);
        means(j) = mean;
        sigmas(j) = sqrt(variance);
        if( sigmas(j)==0 )
        {
            sigmas(j) = 1;
        }
        for(i = 0; i <= npoints-1; i++)
        {
            xyi(i,j) = (xyi(i,j)-means(j))/sigmas(j);
        }
    }
    
    //
    // Internal processing
    //
    lrinternal(xyi, s, npoints, nvars+1, info, lm, ar);
    if( info<0 )
    {
        return;
    }
    
    //
    // Un-standartization
    //
    offs = ap::round(lm.w(3));
    for(j = 0; j <= nvars-1; j++)
    {
        
        //
        // Constant term is updated (and its covariance too,
        // since it gets some variance from J-th component)
        //
        lm.w(offs+nvars) = lm.w(offs+nvars)-lm.w(offs+j)*means(j)/sigmas(j);
        v = means(j)/sigmas(j);
        ap::vsub(&ar.c(nvars, 0), &ar.c(j, 0), ap::vlen(0,nvars), v);
        ap::vsub(ar.c.getcolumn(nvars, 0, nvars), ar.c.getcolumn(j, 0, nvars), v);
        
        //
        // J-th term is updated
        //
        lm.w(offs+j) = lm.w(offs+j)/sigmas(j);
        v = 1/sigmas(j);
        ap::vmul(&ar.c(j, 0), ap::vlen(0,nvars), v);
        ap::vmul(ar.c.getcolumn(j, 0, nvars), v);
    }
}


/*************************************************************************
Like LRBuildS, but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildzs(const ap::real_2d_array& xy,
     const ap::real_1d_array& s,
     int npoints,
     int nvars,
     int& info,
     ap::linearmodel& lm,
     ap::lrreport& ar)
{
    ap::real_2d_array xyi;
    ap::real_1d_array x;
    ap::real_1d_array c;
    int i;
    int j;
    ap::real_value_type v;
    int offs;
    ap::real_value_type mean;
    ap::real_value_type variance;
    ap::real_value_type skewness;
    ap::real_value_type kurtosis;

    
    //
    // Test parameters
    //
    if( npoints<=nvars+1||nvars<1 )
    {
        info = -1;
        return;
    }
    
    //
    // Copy data, add one more column (constant term)
    //
    xyi.setbounds(0, npoints-1, 0, nvars+1);
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&xyi(i, 0), &xy(i, 0), ap::vlen(0,nvars-1));
        xyi(i,nvars) = 0;
        xyi(i,nvars+1) = xy(i,nvars);
    }
    
    //
    // Standartization: unusual scaling
    //
    x.setbounds(0, npoints-1);
    c.setbounds(0, nvars-1);
    for(j = 0; j <= nvars-1; j++)
    {
        ap::vmove(x.getvector(0, npoints-1), xy.getcolumn(j, 0, npoints-1));
        calculatemoments(x, npoints, mean, variance, skewness, kurtosis);
        if( fabs(mean)>sqrt(variance) )
        {
            
            //
            // variation is relatively small, it is better to
            // bring mean value to 1
            //
            c(j) = mean;
        }
        else
        {
            
            //
            // variation is large, it is better to bring variance to 1
            //
            if( variance==0 )
            {
                variance = 1;
            }
            c(j) = sqrt(variance);
        }
        for(i = 0; i <= npoints-1; i++)
        {
            xyi(i,j) = xyi(i,j)/c(j);
        }
    }
    
    //
    // Internal processing
    //
    lrinternal(xyi, s, npoints, nvars+1, info, lm, ar);
    if( info<0 )
    {
        return;
    }
    
    //
    // Un-standartization
    //
    offs = ap::round(lm.w(3));
    for(j = 0; j <= nvars-1; j++)
    {
        
        //
        // J-th term is updated
        //
        lm.w(offs+j) = lm.w(offs+j)/c(j);
        v = 1/c(j);
        ap::vmul(&ar.c(j, 0), ap::vlen(0,nvars), v);
        ap::vmul(ar.c.getcolumn(j, 0, nvars), v);
    }
}


/*************************************************************************
Like LRBuild but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildz(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     ap::linearmodel& lm,
     ap::lrreport& ar)
{
    ap::real_1d_array s;
    int i;
    ap::real_value_type sigma2;

    if( npoints<=nvars+1||nvars<1 )
    {
        info = -1;
        return;
    }
    s.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        s(i) = 1;
    }
    lrbuildzs(xy, s, npoints, nvars, info, lm, ar);
    if( info<0 )
    {
        return;
    }
    sigma2 = ap::sqr(ar.rmserror)*npoints/(npoints-nvars-1);
    for(i = 0; i <= nvars; i++)
    {
        ap::vmul(&ar.c(i, 0), ap::vlen(0,nvars), sigma2);
    }
}


/*************************************************************************
Unpacks coefficients of linear model.

INPUT PARAMETERS:
    LM          -   linear model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables (one less than number
                    of coefficients)

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrunpack(const ap::linearmodel& lm, ap::real_1d_array& v, int& nvars)
{
    int offs;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==lrvnum, "LINREG: Incorrect LINREG version!");
    nvars = ap::round(lm.w(2));
    offs = ap::round(lm.w(3));
    v.setbounds(0, nvars);
    ap::vmove(&v(0), &lm.w(offs), ap::vlen(0,nvars));
}


/*************************************************************************
"Packs" coefficients and creates linear model in ALGLIB format (LRUnpack
reversed).

INPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables

OUTPUT PAREMETERS:
    LM          -   linear model.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrpack(const ap::real_1d_array& v, int nvars, ap::linearmodel& lm)
{
    int offs;

    lm.w.setbounds(0, 4+nvars);
    offs = 4;
    lm.w(0) = 4+nvars+1;
    lm.w(1) = lrvnum;
    lm.w(2) = nvars;
    lm.w(3) = offs;
    ap::vmove(&lm.w(offs), &v(0), ap::vlen(offs,offs+nvars));
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   linear model
    X       -   input vector,  array[0..NVars-1].

Result:
    value of linear model regression estimate

  -- ALGLIB --
     Copyright 03.09.2008 by Bochkanov Sergey
*************************************************************************/
ap::real_value_type lrprocess(const ap::linearmodel& lm, const ap::real_1d_array& x)
{
    ap::real_value_type result;
    ap::real_value_type v;
    int offs;
    int nvars;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==lrvnum, "LINREG: Incorrect LINREG version!");
    nvars = ap::round(lm.w(2));
    offs = ap::round(lm.w(3));
    v = ap::vdotproduct(&x(0), &lm.w(offs), ap::vlen(0,nvars-1));
    result = v+lm.w(offs+nvars);
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
ap::real_value_type lrrmserror(const ap::linearmodel& lm,
     const ap::real_2d_array& xy,
     int npoints)
{
    ap::real_value_type result;
    int i;
    ap::real_value_type v;
    int offs;
    int nvars;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==lrvnum, "LINREG: Incorrect LINREG version!");
    nvars = ap::round(lm.w(2));
    offs = ap::round(lm.w(3));
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        v = ap::vdotproduct(&xy(i, 0), &lm.w(offs), ap::vlen(0,nvars-1));
        v = v+lm.w(offs+nvars);
        result = result+ap::sqr(v-xy(i,nvars));
    }
    result = sqrt(result/npoints);
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
ap::real_value_type lravgerror(const ap::linearmodel& lm,
     const ap::real_2d_array& xy,
     int npoints)
{
    ap::real_value_type result;
    int i;
    ap::real_value_type v;
    int offs;
    int nvars;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==lrvnum, "LINREG: Incorrect LINREG version!");
    nvars = ap::round(lm.w(2));
    offs = ap::round(lm.w(3));
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        v = ap::vdotproduct(&xy(i, 0), &lm.w(offs), ap::vlen(0,nvars-1));
        v = v+lm.w(offs+nvars);
        result = result+fabs(v-xy(i,nvars));
    }
    result = result/npoints;
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
ap::real_value_type lravgrelerror(const ap::linearmodel& lm,
     const ap::real_2d_array& xy,
     int npoints)
{
    ap::real_value_type result;
    int i;
    int k;
    ap::real_value_type v;
    int offs;
    int nvars;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==lrvnum, "LINREG: Incorrect LINREG version!");
    nvars = ap::round(lm.w(2));
    offs = ap::round(lm.w(3));
    result = 0;
    k = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        if( xy(i,nvars)!=0 )
        {
            v = ap::vdotproduct(&xy(i, 0), &lm.w(offs), ap::vlen(0,nvars-1));
            v = v+lm.w(offs+nvars);
            result = result+fabs((v-xy(i,nvars))/xy(i,nvars));
            k = k+1;
        }
    }
    if( k!=0 )
    {
        result = result/k;
    }
    return result;
}


/*************************************************************************
Copying of LinearModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void lrcopy(const ap::linearmodel& lm1, ap::linearmodel& lm2)
{
    int k;

    k = ap::round(lm1.w(0));
    lm2.w.setbounds(0, k-1);
    ap::vmove(&lm2.w(0), &lm1.w(0), ap::vlen(0,k-1));
}


/*************************************************************************
Serialization of LinearModel strucure

INPUT PARAMETERS:
    LM      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores model,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void lrserialize(const ap::linearmodel& lm, ap::real_1d_array& ra, int& rlen)
{

    rlen = ap::round(lm.w(0))+1;
    ra.setbounds(0, rlen-1);
    ra(0) = lrvnum;
    ap::vmove(&ra(1), &lm.w(0), ap::vlen(1,rlen-1));
}


/*************************************************************************
Unserialization of DecisionForest strucure

INPUT PARAMETERS:
    RA      -   real array which stores decision forest

OUTPUT PARAMETERS:
    LM      -   unserialized structure

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void lrunserialize(const ap::real_1d_array& ra, ap::linearmodel& lm)
{

    ap::ap_error::make_assertion(ap::round(ra(0))==lrvnum, "LRUnserialize: incorrect array!");
    lm.w.setbounds(0, ap::round(ra(1))-1);
    ap::vmove(&lm.w(0), &ra(1), ap::vlen(0,ap::round(ra(1))-1));
}


/*************************************************************************
Internal linear regression subroutine
*************************************************************************/
static void lrinternal(const ap::real_2d_array& xy,
     const ap::real_1d_array& s,
     int npoints,
     int nvars,
     int& info,
     ap::linearmodel& lm,
     ap::lrreport& ar)
{
    ap::real_2d_array a;
    ap::real_2d_array u;
    ap::real_2d_array vt;
    ap::real_2d_array vm;
    ap::real_2d_array xym;
    ap::real_1d_array b;
    ap::real_1d_array sv;
    ap::real_1d_array t;
    ap::real_1d_array svi;
    ap::real_1d_array work;
    int i;
    int j;
    int k;
    int ncv;
    int na;
    int nacv;
    ap::real_value_type r;
    ap::real_value_type p;
    ap::real_value_type epstol;
    ap::lrreport ar2;
    int offs;
    ap::linearmodel tlm;

    epstol = 1000;
    
    //
    // Check for errors in data
    //
    if( npoints<nvars||nvars<1 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= npoints-1; i++)
    {
        if( s(i)<=0 )
        {
            info = -2;
            return;
        }
    }
    info = 1;
    
    //
    // Create design matrix
    //
    a.setbounds(0, npoints-1, 0, nvars-1);
    b.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        r = 1/s(i);
        ap::vmove(&a(i, 0), &xy(i, 0), ap::vlen(0,nvars-1), r);
        b(i) = xy(i,nvars)/s(i);
    }
    
    //
    // Allocate W:
    // W[0]     array size
    // W[1]     version number, 0
    // W[2]     NVars (minus 1, to be compatible with external representation)
    // W[3]     coefficients offset
    //
    lm.w.setbounds(0, 4+nvars-1);
    offs = 4;
    lm.w(0) = 4+nvars;
    lm.w(1) = lrvnum;
    lm.w(2) = nvars-1;
    lm.w(3) = offs;
    
    //
    // Solve problem using SVD:
    //
    // 0. check for degeneracy (different types)
    // 1. A = U*diag(sv)*V'
    // 2. T = b'*U
    // 3. w = SUM((T[i]/sv[i])*V[..,i])
    // 4. cov(wi,wj) = SUM(Vji*Vjk/sv[i]^2,K=1..M)
    //
    // see $15.4 of "Numerical Recipes in C" for more information
    //
    t.setbounds(0, nvars-1);
    svi.setbounds(0, nvars-1);
    ar.c.setbounds(0, nvars-1, 0, nvars-1);
    vm.setbounds(0, nvars-1, 0, nvars-1);
    if( !rmatrixsvd(a, npoints, nvars, 1, 1, 2, sv, u, vt) )
    {
        info = -4;
        return;
    }
    if( sv(0)<=0 )
    {
        
        //
        // Degenerate case: zero design matrix.
        //
        for(i = offs; i <= offs+nvars-1; i++)
        {
            lm.w(i) = 0;
        }
        ar.rmserror = lrrmserror(lm, xy, npoints);
        ar.avgerror = lravgerror(lm, xy, npoints);
        ar.avgrelerror = lravgrelerror(lm, xy, npoints);
        ar.cvrmserror = ar.rmserror;
        ar.cvavgerror = ar.avgerror;
        ar.cvavgrelerror = ar.avgrelerror;
        ar.ncvdefects = 0;
        ar.cvdefects.setbounds(0, nvars-1);
        ar.c.setbounds(0, nvars-1, 0, nvars-1);
        for(i = 0; i <= nvars-1; i++)
        {
            for(j = 0; j <= nvars-1; j++)
            {
                ar.c(i,j) = 0;
            }
        }
        return;
    }
    if( sv(nvars-1)<=epstol*ap::machineepsilon*sv(0) )
    {
        
        //
        // Degenerate case, non-zero design matrix.
        //
        // We can leave it and solve task in SVD least squares fashion.
        // Solution and covariance matrix will be obtained correctly,
        // but CV error estimates - will not. It is better to reduce
        // it to non-degenerate task and to obtain correct CV estimates.
        //
        for(k = nvars; k >= 1; k--)
        {
            if( sv(k-1)>epstol*ap::machineepsilon*sv(0) )
            {
                
                //
                // Reduce
                //
                xym.setbounds(0, npoints-1, 0, k);
                for(i = 0; i <= npoints-1; i++)
                {
                    for(j = 0; j <= k-1; j++)
                    {
                        r = ap::vdotproduct(&xy(i, 0), &vt(j, 0), ap::vlen(0,nvars-1));
                        xym(i,j) = r;
                    }
                    xym(i,k) = xy(i,nvars);
                }
                
                //
                // Solve
                //
                lrinternal(xym, s, npoints, k, info, tlm, ar2);
                if( info!=1 )
                {
                    return;
                }
                
                //
                // Convert back to un-reduced format
                //
                for(j = 0; j <= nvars-1; j++)
                {
                    lm.w(offs+j) = 0;
                }
                for(j = 0; j <= k-1; j++)
                {
                    r = tlm.w(offs+j);
                    ap::vadd(&lm.w(offs), &vt(j, 0), ap::vlen(offs,offs+nvars-1), r);
                }
                ar.rmserror = ar2.rmserror;
                ar.avgerror = ar2.avgerror;
                ar.avgrelerror = ar2.avgrelerror;
                ar.cvrmserror = ar2.cvrmserror;
                ar.cvavgerror = ar2.cvavgerror;
                ar.cvavgrelerror = ar2.cvavgrelerror;
                ar.ncvdefects = ar2.ncvdefects;
                ar.cvdefects.setbounds(0, nvars-1);
                for(j = 0; j <= ar.ncvdefects-1; j++)
                {
                    ar.cvdefects(j) = ar2.cvdefects(j);
                }
                ar.c.setbounds(0, nvars-1, 0, nvars-1);
                work.setbounds(1, nvars);
                matrixmatrixmultiply(ar2.c, 0, k-1, 0, k-1, false, vt, 0, k-1, 0, nvars-1, false, 1.0, vm, 0, k-1, 0, nvars-1, 0.0, work);
                matrixmatrixmultiply(vt, 0, k-1, 0, nvars-1, true, vm, 0, k-1, 0, nvars-1, false, 1.0, ar.c, 0, nvars-1, 0, nvars-1, 0.0, work);
                return;
            }
        }
        info = -255;
        return;
    }
    for(i = 0; i <= nvars-1; i++)
    {
        if( sv(i)>epstol*ap::machineepsilon*sv(0) )
        {
            svi(i) = 1/sv(i);
        }
        else
        {
            svi(i) = 0;
        }
    }
    for(i = 0; i <= nvars-1; i++)
    {
        t(i) = 0;
    }
    for(i = 0; i <= npoints-1; i++)
    {
        r = b(i);
        ap::vadd(&t(0), &u(i, 0), ap::vlen(0,nvars-1), r);
    }
    for(i = 0; i <= nvars-1; i++)
    {
        lm.w(offs+i) = 0;
    }
    for(i = 0; i <= nvars-1; i++)
    {
        r = t(i)*svi(i);
        ap::vadd(&lm.w(offs), &vt(i, 0), ap::vlen(offs,offs+nvars-1), r);
    }
    for(j = 0; j <= nvars-1; j++)
    {
        r = svi(j);
        ap::vmove(vm.getcolumn(j, 0, nvars-1), vt.getrow(j, 0, nvars-1), r);
    }
    for(i = 0; i <= nvars-1; i++)
    {
        for(j = i; j <= nvars-1; j++)
        {
            r = ap::vdotproduct(&vm(i, 0), &vm(j, 0), ap::vlen(0,nvars-1));
            ar.c(i,j) = r;
            ar.c(j,i) = r;
        }
    }
    
    //
    // Leave-1-out cross-validation error.
    //
    // NOTATIONS:
    // A            design matrix
    // A*x = b      original linear least squares task
    // U*S*V'       SVD of A
    // ai           i-th row of the A
    // bi           i-th element of the b
    // xf           solution of the original LLS task
    //
    // Cross-validation error of i-th element from a sample is
    // calculated using following formula:
    //
    //     ERRi = ai*xf - (ai*xf-bi*(ui*ui'))/(1-ui*ui')     (1)
    //
    // This formula can be derived from normal equations of the
    // original task
    //
    //     (A'*A)x = A'*b                                    (2)
    //
    // by applying modification (zeroing out i-th row of A) to (2):
    //
    //     (A-ai)'*(A-ai) = (A-ai)'*b
    //
    // and using Sherman-Morrison formula for updating matrix inverse
    //
    // NOTE 1: b is not zeroed out since it is much simpler and
    // does not influence final result.
    //
    // NOTE 2: some design matrices A have such ui that 1-ui*ui'=0.
    // Formula (1) can't be applied for such cases and they are skipped
    // from CV calculation (which distorts resulting CV estimate).
    // But from the properties of U we can conclude that there can
    // be no more than NVars such vectors. Usually
    // NVars << NPoints, so in a normal case it only slightly
    // influences result.
    //
    ncv = 0;
    na = 0;
    nacv = 0;
    ar.rmserror = 0;
    ar.avgerror = 0;
    ar.avgrelerror = 0;
    ar.cvrmserror = 0;
    ar.cvavgerror = 0;
    ar.cvavgrelerror = 0;
    ar.ncvdefects = 0;
    ar.cvdefects.setbounds(0, nvars-1);
    for(i = 0; i <= npoints-1; i++)
    {
        
        //
        // Error on a training set
        //
        r = ap::vdotproduct(&xy(i, 0), &lm.w(offs), ap::vlen(0,nvars-1));
        ar.rmserror = ar.rmserror+ap::sqr(r-xy(i,nvars));
        ar.avgerror = ar.avgerror+fabs(r-xy(i,nvars));
        if( xy(i,nvars)!=0 )
        {
            ar.avgrelerror = ar.avgrelerror+fabs((r-xy(i,nvars))/xy(i,nvars));
            na = na+1;
        }
        
        //
        // Error using fast leave-one-out cross-validation
        //
        p = ap::vdotproduct(&u(i, 0), &u(i, 0), ap::vlen(0,nvars-1));
        if( p>1-epstol*ap::machineepsilon )
        {
            ar.cvdefects(ar.ncvdefects) = i;
            ar.ncvdefects = ar.ncvdefects+1;
            continue;
        }
        r = s(i)*(r/s(i)-b(i)*p)/(1-p);
        ar.cvrmserror = ar.cvrmserror+ap::sqr(r-xy(i,nvars));
        ar.cvavgerror = ar.cvavgerror+fabs(r-xy(i,nvars));
        if( xy(i,nvars)!=0 )
        {
            ar.cvavgrelerror = ar.cvavgrelerror+fabs((r-xy(i,nvars))/xy(i,nvars));
            nacv = nacv+1;
        }
        ncv = ncv+1;
    }
    if( ncv==0 )
    {
        
        //
        // Something strange: ALL ui are degenerate.
        // Unexpected...
        //
        info = -255;
        return;
    }
    ar.rmserror = sqrt(ar.rmserror/npoints);
    ar.avgerror = ar.avgerror/npoints;
    if( na!=0 )
    {
        ar.avgrelerror = ar.avgrelerror/na;
    }
    ar.cvrmserror = sqrt(ar.cvrmserror/ncv);
    ar.cvavgerror = ar.cvavgerror/ncv;
    if( nacv!=0 )
    {
        ar.cvavgrelerror = ar.cvavgrelerror/nacv;
    }
}



