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
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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

#include "sevd.h"

/*************************************************************************
Finding the eigenvalues and eigenvectors of a symmetric matrix

The algorithm finds eigen pairs of a symmetric matrix by reducing it to
tridiagonal form and using the QL/QR algorithm.

Input parameters:
    A       -   symmetric matrix which is given by its upper or lower
                triangular part.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.
    ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                If ZNeeded is equal to:
                 * 0, the eigenvectors are not returned;
                 * 1, the eigenvectors are returned.

Output parameters:
    D       -   eigenvalues in ascending order.
                Array whose index ranges within [0..N-1].
    Z       -   if ZNeeded is equal to:
                 * 0, Z hasn’t changed;
                 * 1, Z contains the eigenvectors.
                Array whose indexes range within [0..N-1, 0..N-1].
                The eigenvectors are stored in the matrix columns.

Result:
    True, if the algorithm has converged.
    False, if the algorithm hasn't converged (rare case).

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************/
bool smatrixevd(ap::real_2d_array a,
     int n,
     int zneeded,
     bool isupper,
     ap::real_1d_array& d,
     ap::real_2d_array& z)
{
    bool result;
    ap::real_1d_array tau;
    ap::real_1d_array e;

    ap::ap_error::make_assertion(zneeded==0||zneeded==1, "SMatrixEVD: incorrect ZNeeded");
    smatrixtd(a, n, isupper, tau, d, e);
    if( zneeded==1 )
    {
        smatrixtdunpackq(a, n, isupper, tau, z);
    }
    result = smatrixtdevd(d, e, n, zneeded, z);
    return result;
}
