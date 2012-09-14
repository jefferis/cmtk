/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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
//  $Revision: 1317 $
//
//  $LastChangedDate: 2010-03-10 14:08:16 -0800 (Wed, 10 Mar 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class TFloat>
EigenValuesSymmetricMatrix<TFloat>
::EigenValuesSymmetricMatrix( const SymmetricMatrix<TFloat>& matrix )
  : m_Eigenvalues( matrix.Dim() )
{
  const int n = static_cast<int>( matrix.Dim() );

  /*  Convert SymmetricMatrix to ap::real_2d_array
   */
  ap::real_2d_array apMatrix;
  apMatrix.setbounds(0, n-1, 0, n-1 );
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      apMatrix(i,j) = static_cast<double>( matrix(i,j) );
  
  ap::real_1d_array apEigenvalues;
  apEigenvalues.setbounds( 0, n-1 );

  ap::real_2d_array apEigenvectorsDummy;

  /*  Run AlgLib eigensystem computation
   */
  if ( ! smatrixevd(apMatrix, (int)n, 1 /*upper storage*/, false /*eigenvectors not needed*/, apEigenvalues, apEigenvectorsDummy) )
    {
    StdErr << "WARNING: smatrixevd did not converge\n";
    }
  
  /*  Convert ap::real_1d_array and ap::real_2d_array
   *  back to our data types.
   */
  for (int i = 0; i < n; i++)
    this->m_Eigenvalues[i] = static_cast<TFloat>( apEigenvalues(i) );
}

} // namespace cmtk
