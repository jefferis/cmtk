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
/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee. All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
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


#include "hessenberg.h"


/*************************************************************************
Obsolete 1-based subroutine.
See RMatrixHessenberg for 0-based replacement.
*************************************************************************/
void toupperhessenberg(ap::real_2d_array& a, int n, ap::real_1d_array& tau)
{
    int i;
    int ip1;
    int nmi;
//    ap::real_value_type aii;
    ap::real_value_type v;
    ap::real_1d_array t;
    ap::real_1d_array work;

    ap::ap_error::make_assertion(n>=0, "ToUpperHessenberg: incorrect N!");
    
    //
    // Quick return if possible
    //
    if( n<=1 )
    {
        return;
    }
    tau.setbounds(1, n-1);
    t.setbounds(1, n);
    work.setbounds(1, n);
    for(i = 1; i <= n-1; i++)
    {
        
        //
        // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
        //
        ip1 = i+1;
        nmi = n-i;
        ap::vmove(t.getvector(1, nmi), a.getcolumn(i, ip1, n));
        generatereflection(t, nmi, v);
        ap::vmove(a.getcolumn(i, ip1, n), t.getvector(1, nmi));
        tau(i) = v;
        t(1) = 1;
        
        //
        // Apply H(i) to A(1:ihi,i+1:ihi) from the right
        //
        applyreflectionfromtheright(a, v, t, 1, n, i+1, n, work);
        
        //
        // Apply H(i) to A(i+1:ihi,i+1:n) from the left
        //
        applyreflectionfromtheleft(a, v, t, i+1, n, i+1, n, work);
    }
}


/*************************************************************************
Obsolete 1-based subroutine.
See RMatrixHessenbergUnpackQ for 0-based replacement.
*************************************************************************/
void unpackqfromupperhessenberg(const ap::real_2d_array& a,
     int n,
     const ap::real_1d_array& tau,
     ap::real_2d_array& q)
{
    int i;
    int j;
    ap::real_1d_array v;
    ap::real_1d_array work;
    int ip1;
    int nmi;

    if( n==0 )
    {
        return;
    }
    
    //
    // init
    //
    q.setbounds(1, n, 1, n);
    v.setbounds(1, n);
    work.setbounds(1, n);
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    
    //
    // unpack Q
    //
    for(i = 1; i <= n-1; i++)
    {
        
        //
        // Apply H(i)
        //
        ip1 = i+1;
        nmi = n-i;
        ap::vmove(v.getvector(1, nmi), a.getcolumn(i, ip1, n));
        v(1) = 1;
        applyreflectionfromtheright(q, tau(i), v, 1, n, i+1, n, work);
    }
}


/*************************************************************************
Obsolete 1-based subroutine.
See RMatrixHessenbergUnpackH for 0-based replacement.
*************************************************************************/
void unpackhfromupperhessenberg(const ap::real_2d_array& a,
     int n,
     const ap::real_1d_array&, //tau
     ap::real_2d_array& h)
{
    int i;
    int j;
    ap::real_1d_array v;
    ap::real_1d_array work;
//    int ip1;
//    int nmi;

    if( n==0 )
    {
        return;
    }
    h.setbounds(1, n, 1, n);
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= i-2; j++)
        {
            h(i,j) = 0;
        }
        j = ap::maxint(1, i-1);
        ap::vmove(&h(i, j), &a(i, j), ap::vlen(j,n));
    }
}



