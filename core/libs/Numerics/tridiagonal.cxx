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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/
/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

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


#include "tridiagonal.h"

/*************************************************************************
Reduction of a symmetric matrix which is given by its higher or lower
triangular part to a tridiagonal matrix using orthogonal similarity
transformation: Q'*A*Q=T.

Input parameters:
    A       -   matrix to be transformed
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then matrix A is given
                by its upper triangle, and the lower triangle is not used
                and not modified by the algorithm, and vice versa
                if IsUpper = False.

Output parameters:
    A       -   matrices T and Q in  compact form (see lower)
    Tau     -   array of factors which are forming matrices H(i)
                array with elements [0..N-2].
    D       -   main diagonal of symmetric matrix T.
                array with elements [0..N-1].
    E       -   secondary diagonal of symmetric matrix T.
                array with elements [0..N-2].


  If IsUpper=True, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-2) . . . H(2) H(0).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
  A(0:i-1,i+1), and tau in TAU(i).

  If IsUpper=False, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(0) H(2) . . . H(n-2).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  d   e   v1  v2  v3 )              (  d                  )
    (      d   e   v2  v3 )              (  e   d              )
    (          d   e   v3 )              (  v0  e   d          )
    (              d   e  )              (  v0  v1  e   d      )
    (                  d  )              (  v0  v1  v2  e   d  )

  where d and e denote diagonal and off-diagonal elements of T, and vi
  denotes an element of the vector defining H(i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************/
void smatrixtd(ap::real_2d_array& a,
     int n,
     bool isupper,
     ap::real_1d_array& tau,
     ap::real_1d_array& d,
     ap::real_1d_array& e)
{
    int i;
    ap::real_value_type alpha;
    ap::real_value_type taui;
    ap::real_value_type v;
    ap::real_1d_array t;
    ap::real_1d_array t2;
    ap::real_1d_array t3;

    if( n<=0 )
    {
        return;
    }
    t.setbounds(1, n);
    t2.setbounds(1, n);
    t3.setbounds(1, n);
    if( n>1 )
    {
        tau.setbounds(0, n-2);
    }
    d.setbounds(0, n-1);
    if( n>1 )
    {
        e.setbounds(0, n-2);
    }
    if( isupper )
    {
        
        //
        // Reduce the upper triangle of A
        //
        for(i = n-2; i >= 0; i--)
        {
            
            //
            // Generate elementary reflector H() = E - tau * v * v'
            //
            if( i>=1 )
            {
                ap::vmove(t.getvector(2, i+1), a.getcolumn(i+1, 0, i-1));
            }
            t(1) = a(i,i+1);
            generatereflection(t, i+1, taui);
            if( i>=1 )
            {
                ap::vmove(a.getcolumn(i+1, 0, i-1), t.getvector(2, i+1));
            }
            a(i,i+1) = t(1);
            e(i) = a(i,i+1);
            if( taui!=0 )
            {
                
                //
                // Apply H from both sides to A
                //
                a(i,i+1) = 1;
                
                //
                // Compute  x := tau * A * v  storing x in TAU
                //
                ap::vmove(t.getvector(1, i+1), a.getcolumn(i+1, 0, i));
                symmetricmatrixvectormultiply(a, isupper, 0, i, t, taui, t3);
                ap::vmove(&tau(0), &t3(1), ap::vlen(0,i));
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                v = ap::vdotproduct(tau.getvector(0, i), a.getcolumn(i+1, 0, i));
                alpha = -0.5*taui*v;
                ap::vadd(tau.getvector(0, i), a.getcolumn(i+1, 0, i), alpha);
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                ap::vmove(t.getvector(1, i+1), a.getcolumn(i+1, 0, i));
                ap::vmove(&t3(1), &tau(0), ap::vlen(1,i+1));
                symmetricrank2update(a, isupper, 0, i, t, t3, t2, ap::real_value_type(-1));
                a(i,i+1) = e(i);
            }
            d(i+1) = a(i+1,i+1);
            tau(i) = taui;
        }
        d(0) = a(0,0);
    }
    else
    {
        
        //
        // Reduce the lower triangle of A
        //
        for(i = 0; i <= n-2; i++)
        {
            
            //
            // Generate elementary reflector H = E - tau * v * v'
            //
            ap::vmove(t.getvector(1, n-i-1), a.getcolumn(i, i+1, n-1));
            generatereflection(t, n-i-1, taui);
            ap::vmove(a.getcolumn(i, i+1, n-1), t.getvector(1, n-i-1));
            e(i) = a(i+1,i);
            if( taui!=0 )
            {
                
                //
                // Apply H from both sides to A
                //
                a(i+1,i) = 1;
                
                //
                // Compute  x := tau * A * v  storing y in TAU
                //
                ap::vmove(t.getvector(1, n-i-1), a.getcolumn(i, i+1, n-1));
                symmetricmatrixvectormultiply(a, isupper, i+1, n-1, t, taui, t2);
                ap::vmove(&tau(i), &t2(1), ap::vlen(i,n-2));
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                v = ap::vdotproduct(tau.getvector(i, n-2), a.getcolumn(i, i+1, n-1));
                alpha = -0.5*taui*v;
                ap::vadd(tau.getvector(i, n-2), a.getcolumn(i, i+1, n-1), alpha);
                
                //
                // Apply the transformation as a rank-2 update:
                //     A := A - v * w' - w * v'
                //
                //
                ap::vmove(t.getvector(1, n-i-1), a.getcolumn(i, i+1, n-1));
                ap::vmove(&t2(1), &tau(i), ap::vlen(1,n-i-1));
                symmetricrank2update(a, isupper, i+1, n-1, t, t2, t3, ap::real_value_type(-1));
                a(i+1,i) = e(i);
            }
            d(i) = a(i,i);
            tau(i) = taui;
        }
        d(n-1) = a(n-1,n-1);
    }
}


/*************************************************************************
Unpacking matrix Q which reduces symmetric matrix to a tridiagonal
form.

Input parameters:
    A       -   the result of a SMatrixTD subroutine
    N       -   size of matrix A.
    IsUpper -   storage format (a parameter of SMatrixTD subroutine)
    Tau     -   the result of a SMatrixTD subroutine

Output parameters:
    Q       -   transformation matrix.
                array with elements [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************/
void smatrixtdunpackq(const ap::real_2d_array& a,
     const int& n,
     const bool& isupper,
     const ap::real_1d_array& tau,
     ap::real_2d_array& q)
{
    int i;
    int j;
    ap::real_1d_array v;
    ap::real_1d_array work;

    if( n==0 )
    {
        return;
    }
    
    //
    // init
    //
    q.setbounds(0, n-1, 0, n-1);
    v.setbounds(1, n);
    work.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
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
    if( isupper )
    {
        for(i = 0; i <= n-2; i++)
        {
            
            //
            // Apply H(i)
            //
            ap::vmove(v.getvector(1, i+1), a.getcolumn(i+1, 0, i));
            v(i+1) = 1;
            applyreflectionfromtheleft(q, tau(i), v, 0, i, 0, n-1, work);
        }
    }
    else
    {
        for(i = n-2; i >= 0; i--)
        {
            
            //
            // Apply H(i)
            //
            ap::vmove(v.getvector(1, n-i-1), a.getcolumn(i, i+1, n-1));
            v(1) = 1;
            applyreflectionfromtheleft(q, tau(i), v, i+1, n-1, 0, n-1, work);
        }
    }
}
