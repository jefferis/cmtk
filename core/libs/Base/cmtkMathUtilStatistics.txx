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

#include <sevd.h>
#include <studentttests.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

namespace
MathUtil
{

template<class T>                                                                     
T                                                                                     
Mean                                                                                  
( const unsigned int nValues, const T* values )                                       
{                                                                                     
  T mean = 0.0;                                                                       
                                                                                      
  for ( unsigned int j = 0; j < nValues; j++ )                                        
    mean += values[j];                                                                
  mean /= nValues;                                                                    
                                                                                      
  return mean;                                                                        
}                                                                                     
                                                                                      
template<class T>                                                                     
T                                                                                     
Variance                                                                              
( const unsigned int nValues, const T* values, const T mean, const bool unbiased )           
{                                                                                            
  T sumOfSquares = 0.0;                                                                          
                                                                                             
  T sum = 0.0;                                                                               
  for ( unsigned int j = 0; j < nValues; j++ ) 
    { 
    const T s = values[j] - mean;                                                                  
    sum += s;                                                                                
    sumOfSquares += s*s;                                                                         
    }                                                                                          
  
  if ( unbiased )                                                                            
    return (sumOfSquares - sum*sum/nValues) / (nValues-1);                                   
  else                                                                                       
    return (sumOfSquares - sum*sum/nValues) / (nValues);                                     
}        

template<class T>
T
Correlation
( const size_t n, const T* x, const T* y )
{
  // compute means
  T meanx = 0, meany = 0;
  for ( size_t i = 0; i < n; ++i )
    {
      meanx += x[i];
      meany += y[i];
    }
  meanx /= n;
  meany /= n;

  // compute parameter correlations
  T c = 0, xSq = 0, ySq = 0;
  T dx, dy;
  for ( size_t i = 0; i < n; ++i )
    {
    c += (dx=x[i]-meanx) * (dy=y[i]-meany);
    xSq += dx * dx;
    ySq += dy * dy;
    }
  
  return static_cast<T>( c / (sqrt(xSq*ySq)+1e-20) );
}
  
template<class T>
T
TTest
( const unsigned int nValuesX, const T* valuesX, const unsigned int nValuesY, const T* valuesY, T& t )
{
  T averageX, averageY;
  return TTest( nValuesX, valuesX, nValuesY, valuesY, t, averageX, averageY );
}

template<class T>
T
PairedTTest
( const unsigned int nValues, const T* valuesX, const T* valuesY, T& t, T& avgX, T& avgY )
{
  avgX = Mean( nValues, valuesX );
  avgY = Mean( nValues, valuesY );

  T SSD = 0;
  for ( size_t i = 0; i < nValues; ++i )
    SSD += Square( (valuesX[i]-avgX) - (valuesY[i]-avgY) );

  t = (avgX - avgY) * sqrt( (nValues * (nValues-1)) / SSD );
  
  double s = studenttdistribution(nValues-1, t);
  double p1 = 2 * ap::minreal(s, 1-s);

  return (T) p1; // probability
}

/** Two-sample t-test
 */
template<class T>
T
TTest
( const unsigned int nValuesX, const T* valuesX, const unsigned int nValuesY, const T* valuesY, T& t, T& avgX, T& avgY )
{
  
  ap::real_1d_array apValuesX;
  apValuesX.setbounds( 0, nValuesX-1 );
  for (int i = 0; i < nValuesX; i++)
    apValuesX(i) = (double)(1.0 * valuesX[i]);

  ap::real_1d_array apValuesY;
  apValuesY.setbounds( 0, nValuesY-1 );
  for (int i = 0; i < nValuesY; i++)
    apValuesY(i) = (double)(1.0 * valuesY[i]);
  
  double t_temp, p1, p2, p3;

  avgX = MathUtil::Mean<T>( nValuesX, valuesX );
  avgY = MathUtil::Mean<T>( nValuesY, valuesY );
  
  studentttest2( apValuesX, nValuesX, apValuesY, nValuesY, t_temp, p1, p2, p3 );

  t = (T) t_temp;

  return (T) p1; // probability

}

/** One-sample t-test
 */
template<class T>
T
TTest
( const unsigned int nValuesX, const T* valuesX, T& t, T& avgX )
{

  avgX = Mean( nValuesX, valuesX );
  T varianceX = Variance( nValuesX, valuesX, avgX, true /*unbiased*/ );

  t = avgX * nValuesX / sqrt( varianceX );

  double s = studenttdistribution(nValuesX-1, t);
  double p1 = 2 * ap::minreal(s, 1-s);

  return (T) p1; // probability
}

} // namespace MathUtil
} // namespace cmtk
