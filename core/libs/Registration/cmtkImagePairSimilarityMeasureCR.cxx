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

#include <cmtkImagePairSimilarityMeasureCR.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairSimilarityMeasureCR::ImagePairSimilarityMeasureCR
( const UniformVolume::SmartPtr& refVolume, const UniformVolume::SmartPtr& fltVolume, const Interpolators::InterpolationEnum interpolation )
    : ImagePairSimilarityMeasure( refVolume, fltVolume, interpolation )
{ 
  NumBinsX = std::max<unsigned>( std::min<unsigned>( refVolume->GetNumberOfPixels(), 128 ), 8 );
  HistogramI.Resize( NumBinsX );
  
  NumBinsY = std::max<unsigned>( std::min<unsigned>( fltVolume->GetNumberOfPixels(), 128 ), 8 );
  HistogramJ.Resize( NumBinsY );
  
  HistogramI.SetRange( refVolume->GetData()->GetRange() );
  
  SumJ.resize( NumBinsX );
  SumJ2.resize( NumBinsX );
  
  fltVolume->GetData()->GetStatistics( MuJ, SigmaSqJ );
  
  HistogramJ.SetRange( fltVolume->GetData()->GetRange() );
  
  SumI.resize( NumBinsY );
  SumI2.resize( NumBinsY );
  
  refVolume->GetData()->GetStatistics( MuI, SigmaSqI );
}

ImagePairSimilarityMeasureCR::ReturnType
ImagePairSimilarityMeasureCR::Get () const
{
  const double invSampleCount = 1.0 / HistogramI.SampleCount();
  // initialize variable for the weighted sum of the sigma^2 values over all
  // reference intensity classes.
  double sumSigmaSquare = 0;
  // run over all bins, i.e., reference classes
  for ( unsigned int j = 0; j < NumBinsX; ++j ) 
    {
    // are there any values in the current class?
    if ( HistogramI[j] ) 
      {
      // compute mean floating value for this reference class
      const double mu = SumJ[j] / HistogramI[j];
      // compute variance of floating values for this reference class
      const double sigmaSq = ( mu*mu*HistogramI[j] - 2.0*mu*SumJ[j] + SumJ2[j] ) / HistogramI[j]; 
      // update sum over all classes with weighted sigma^2 for this class.
      sumSigmaSquare += (invSampleCount * HistogramI[j]) * sigmaSq;
      }
    }
  
  // compute (supposedly) correlation ratio
  Self::ReturnType cr = static_cast<Self::ReturnType>( 1.0 - (1.0 /  SigmaSqJ ) * sumSigmaSquare );
  
  sumSigmaSquare = 0;
  for ( unsigned int i = 0; i < NumBinsY; ++i ) 
    {
    if ( HistogramJ[i] ) 
      {
      const double mu = SumI[i] / HistogramJ[i];
      const double sigmaSq = ( mu*mu*HistogramJ[i] - 2.0*mu*SumI[i] + SumI2[i] ) / HistogramJ[i]; 
      // update sum over all classes with weighted sigma^2 for this class.
      sumSigmaSquare += (invSampleCount * HistogramJ[i]) * sigmaSq;
      }
    }
  
  // add reverse correlation ratio
  cr += static_cast<Self::ReturnType>(1.0 - (1.0 /  SigmaSqI ) * sumSigmaSquare);
  
  return cr;
}

} // namespace cmtk
