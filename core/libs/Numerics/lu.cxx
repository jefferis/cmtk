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


#include "lu.h"

static const int lunb = 8;

void rmatrixlu2(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots);

/*************************************************************************
LU decomposition of a general matrix of size MxN

The subroutine calculates the LU decomposition of a rectangular general
matrix with partial pivoting (with row permutations).

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices L and U in compact form (see below).
            Array whose indexes range within [0..M-1, 0..N-1].
    Pivots - permutation matrix in compact form (see below).
            Array whose index ranges within [0..Min(M-1,N-1)].

Matrix A is represented as A = P * L * U, where P is a permutation matrix,
matrix L - lower triangular (or lower trapezoid, if M>N) matrix,
U - upper triangular (or upper trapezoid, if M<N) matrix.

Let M be equal to 4 and N be equal to 3:

                   (  1          )    ( U11 U12 U13  )
A = P1 * P2 * P3 * ( L21  1      )  * (     U22 U23  )
                   ( L31 L32  1  )    (         U33  )
                   ( L41 L42 L43 )

Matrix L has size MxMin(M,N), matrix U has size Min(M,N)xN, matrix P(i) is
a permutation of the identity matrix of size MxM with numbers I and Pivots[I].

The algorithm returns array Pivots and the following matrix which replaces
matrix A and contains matrices L and U in compact form (the example applies
to M=4, N=3).

 ( U11 U12 U13 )
 ( L21 U22 U23 )
 ( L31 L32 U33 )
 ( L41 L42 L43 )

As we can see, the unit diagonal isn't stored.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
void rmatrixlu(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    ap::real_2d_array b;
    ap::real_1d_array t;
    ap::integer_1d_array bp;
    int minmn;
    int i;
    int ip;
    int j;
    int j1;
    int j2;
    int cb;
    int nb;
    ap::real_value_type v;

#ifndef NO_AP_ASSERT
    ap::ap_error::make_assertion(lunb>=1, "RMatrixLU internal error");
#endif

    nb = lunb;
    
    //
    // Decide what to use - blocked or unblocked code
    //
    if( n<=1||ap::minint(m, n)<=nb||nb==1 )
    {
        
        //
        // Unblocked code
        //
        rmatrixlu2(a, m, n, pivots);
    }
    else
    {
        
        //
        // Blocked code.
        // First, prepare temporary matrix and indices
        //
        b.setbounds(0, m-1, 0, nb-1);
        t.setbounds(0, n-1);
        pivots.setbounds(0, ap::minint(m, n)-1);
        minmn = ap::minint(m, n);
        j1 = 0;
        j2 = ap::minint(minmn, nb)-1;
        
        //
        // Main cycle
        //
        while(j1<minmn)
        {
            cb = j2-j1+1;
            
            //
            // LU factorization of diagonal and subdiagonal blocks:
            // 1. Copy columns J1..J2 of A to B
            // 2. LU(B)
            // 3. Copy result back to A
            // 4. Copy pivots, apply pivots
            //
            for(i = j1; i <= m-1; i++)
            {
                ap::vmove(&b(i-j1, 0), &a(i, j1), ap::vlen(0,cb-1));
            }
            rmatrixlu2(b, m-j1, cb, bp);
            for(i = j1; i <= m-1; i++)
            {
                ap::vmove(&a(i, j1), &b(i-j1, 0), ap::vlen(j1,j2));
            }
            for(i = 0; i <= cb-1; i++)
            {
                ip = bp(i);
                pivots(j1+i) = j1+ip;
                if( bp(i)!=i )
                {
                    if( j1!=0 )
                    {
                        
                        //
                        // Interchange columns 0:J1-1
                        //
                        ap::vmove(&t(0), &a(j1+i, 0), ap::vlen(0,j1-1));
                        ap::vmove(&a(j1+i, 0), &a(j1+ip, 0), ap::vlen(0,j1-1));
                        ap::vmove(&a(j1+ip, 0), &t(0), ap::vlen(0,j1-1));
                    }
                    if( j2<n-1 )
                    {
                        
                        //
                        // Interchange the rest of the matrix, if needed
                        //
                        ap::vmove(&t(j2+1), &a(j1+i, j2+1), ap::vlen(j2+1,n-1));
                        ap::vmove(&a(j1+i, j2+1), &a(j1+ip, j2+1), ap::vlen(j2+1,n-1));
                        ap::vmove(&a(j1+ip, j2+1), &t(j2+1), ap::vlen(j2+1,n-1));
                    }
                }
            }
            
            //
            // Compute block row of U
            //
            if( j2<n-1 )
            {
                for(i = j1+1; i <= j2; i++)
                {
                    for(j = j1; j <= i-1; j++)
                    {
                        v = a(i,j);
                        ap::vsub(&a(i, j2+1), &a(j, j2+1), ap::vlen(j2+1,n-1), v);
                    }
                }
            }
            
            //
            // Update trailing submatrix
            //
            if( j2<n-1 )
            {
                for(i = j2+1; i <= m-1; i++)
                {
                    for(j = j1; j <= j2; j++)
                    {
                        v = a(i,j);
                        ap::vsub(&a(i, j2+1), &a(j, j2+1), ap::vlen(j2+1,n-1), v);
                    }
                }
            }
            
            //
            // Next step
            //
            j1 = j2+1;
            j2 = ap::minint(minmn, j1+nb)-1;
        }
    }
}

/*************************************************************************
Level 2 BLAS version of RMatrixLU

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
void rmatrixlu2(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    int i;
    int j;
    int jp;
    ap::real_1d_array t1;
    ap::real_value_type s;

    pivots.setbounds(0, ap::minint(m-1, n-1));
    t1.setbounds(0, ap::maxint(m-1, n-1));

#ifndef NO_AP_ASSERT
    ap::ap_error::make_assertion(m>=0&&n>=0, "Error in LUDecomposition: incorrect function arguments");
#endif
    
    //
    // Quick return if possible
    //
    if( m==0||n==0 )
    {
        return;
    }
    for(j = 0; j <= ap::minint(m-1, n-1); j++)
    {
        
        //
        // Find pivot and test for singularity.
        //
        jp = j;
        for(i = j+1; i <= m-1; i++)
        {
            if( fabs(a(i,j))>fabs(a(jp,j)) )
            {
                jp = i;
            }
        }
        pivots(j) = jp;
        if( a(jp,j)!=0 )
        {
            
            //
            //Apply the interchange to rows
            //
            if( jp!=j )
            {
                ap::vmove(&t1(0), &a(j, 0), ap::vlen(0,n-1));
                ap::vmove(&a(j, 0), &a(jp, 0), ap::vlen(0,n-1));
                ap::vmove(&a(jp, 0), &t1(0), ap::vlen(0,n-1));
            }
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if( j<m )
            {
                jp = j+1;
                s = 1/a(j,j);
                ap::vmul(a.getcolumn(j, jp, m-1), s);
            }
        }
        if( j<ap::minint(m, n)-1 )
        {
            
            //
            //Update trailing submatrix.
            //
            jp = j+1;
            for(i = j+1; i <= m-1; i++)
            {
                s = a(i,j);
                ap::vsub(&a(i, jp), &a(j, jp), ap::vlen(jp,n-1), s);
            }
        }
    }
}



