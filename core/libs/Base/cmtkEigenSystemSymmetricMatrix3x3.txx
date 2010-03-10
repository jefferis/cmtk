/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class TFloat>
EigenSystemSymmetricMatrix3x3<TFloat>
::EigenSystemSymmetricMatrix3x3( const Matrix3x3<TFloat>& matrix, const bool sortAbsolute )
{
  TFloat e[3];
  for (int i = 0; i < 3; i++) 
    {
    for (int j = 0; j < 3; j++) 
      {
      this->m_Eigenvectors[i][j] = matrix[i][j];
      }
    }
  tred2( this->m_Eigenvectors, this->m_Eigenvalues, e);
  tql2( this->m_Eigenvectors, this->m_Eigenvalues, e, sortAbsolute );
}

template<class TFloat>
TFloat
EigenSystemSymmetricMatrix3x3<TFloat>
::hypot2( const TFloat& x, const TFloat& y) 
{
  return sqrt(x*x+y*y);
}

template<class TFloat>
void 
EigenSystemSymmetricMatrix3x3<TFloat>
::tred2( TFloat V[3][3], TFloat d[3], TFloat e[3] ) 
{
//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.
  
  for (int j = 0; j < 3; j++) 
    {
    d[j] = V[2][j];
    }
  
  // Householder reduction to tridiagonal form.
  
  for (int i = 2; i > 0; i--) 
    {    
    // Scale to avoid under/overflow.
    
    TFloat scale = 0.0;
    TFloat h = 0.0;
    for (int k = 0; k < i; k++) 
      {
      scale = scale + fabs(d[k]);
      }
    if (scale == 0.0) 
      {
      e[i] = d[i-1];
      for (int j = 0; j < i; j++) 
	{
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
	}
      } 
    else 
      {
      // Generate Householder vector.
      
      for (int k = 0; k < i; k++) 
	{
        d[k] /= scale;
        h += d[k] * d[k];
	}
      TFloat f = d[i-1];
      TFloat g = sqrt(h);
      if (f > 0) 
	{
        g = -g;
	}
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (int j = 0; j < i; j++) 
	{
	e[j] = 0.0;
	}
      
      // Apply similarity transformation to remaining columns.

      for (int j = 0; j < i; j++) 
	{
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (int k = j+1; k <= i-1; k++) 
	  {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
	  }
        e[j] = g;
	}
      f = 0.0;
      for (int j = 0; j < i; j++) 
	{
        e[j] /= h;
        f += e[j] * d[j];
	}
      TFloat hh = f / (h + h);
      for (int j = 0; j < i; j++) 
	{
        e[j] -= hh * d[j];
	}
      for (int j = 0; j < i; j++) 
	{
	f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) 
	  {
	  V[k][j] -= (f * e[k] + g * d[k]);
	  }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
	}
      }
    d[i] = h;
    }
  
  // Accumulate transformations.
  
  for (int i = 0; i < 2; i++) 
    {
    V[2][i] = V[i][i];
    V[i][i] = 1.0;
    TFloat h = d[i+1];
    if (h != 0.0)
      {
      for (int k = 0; k <= i; k++) 
	{
        d[k] = V[k][i+1] / h;
	}
      for (int j = 0; j <= i; j++) 
	{
        TFloat g = 0.0;
        for (int k = 0; k <= i; k++) 
	  {
	  g += V[k][i+1] * V[k][j];
	  }
        for (int k = 0; k <= i; k++) 
	  {
          V[k][j] -= g * d[k];
	  }
	}
      }
    for (int k = 0; k <= i; k++) 
      {
      V[k][i+1] = 0.0;
      }
    }
  for (int j = 0; j < 3; j++) 
    {
    d[j] = V[2][j];
    V[2][j] = 0.0;
    }
  V[2][2] = 1.0;
  e[0] = 0.0;
} 

template<class TFloat>
void 
EigenSystemSymmetricMatrix3x3<TFloat>
::tql2(TFloat V[3][3], TFloat d[3], TFloat e[3], const bool sortAbsolute) 
{
//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.
  
  for (int i = 1; i < 3; i++) 
    {
    e[i-1] = e[i];
    }
  e[2] = 0.0;
  
  TFloat f = 0.0;
  TFloat tst1 = 0.0;
  TFloat eps = pow(2.0,-52.0);
  for (int l = 0; l < 3; l++) 
    {
    // Find small subdiagonal element
    
    tst1 = ( tst1 > fabs( d[l]) + fabs(e[l]) ) ? tst1 : fabs( d[l]) + fabs(e[l]) ;
    int m = l;
    while (m < 3) 
      {
      if (fabs(e[m]) <= eps*tst1) 
	{
        break;
	}
      m++;
      }
    
    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.
    
    if (m > l) 
      {
      int iter = 0;
      do 
	{
	iter = iter + 1;  // (Could check iteration count here.)
	
        // Compute implicit shift
	
        TFloat g = d[l];
        TFloat p = (d[l+1] - g) / (2.0 * e[l]);
        TFloat r = hypot2(p,1.0);
        if (p < 0) 
	  {
	  r = -r;
	  }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        TFloat dl1 = d[l+1];
        TFloat h = g - d[l];
        for (int i = l+2; i < 3; i++) 
	  {
          d[i] -= h;
	  }
        f = f + h;
	
        // Implicit QL transformation.
	
        p = d[m];
        TFloat c = 1.0;
        TFloat c2 = c;
        TFloat c3 = c;
        TFloat el1 = e[l+1];
        TFloat s = 0.0;
        TFloat s2 = 0.0;
        for (int i = m-1; i >= l; i--) 
	  {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);
	  
          // Accumulate transformation.
	  
          for (int k = 0; k < 3; k++) 
	    {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
	    }
	  }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
	
        // Check for convergence.
	
	} 
      while (fabs(e[l]) > eps*tst1);
      }
    d[l] = d[l] + f;
    e[l] = 0.0;
    }
  
  // Sort eigenvalues and corresponding vectors.
  
  for (int i = 0; i < 2; i++) 
    {
    int k = i;
    TFloat p = d[i];
    for (int j = i+1; j < 3; j++) 
      {
      const bool swap = sortAbsolute ? (fabs(d[j]) < fabs(p)) : (d[j] < p);
      if ( swap ) 
	{
	k = j;
        p = d[j];
	}
      }
    if (k != i) 
      {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < 3; j++) 
	{
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
	}
      }
    }
}

} // namespace cmtk
