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

#include <cmtkLabelCombinationMultiClassSTAPLE.h>
#include <cmtkLabelCombinationVoting.h>

#include <cmtkProgress.h>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

LabelCombinationMultiClassSTAPLE
::LabelCombinationMultiClassSTAPLE
( const std::vector<TypedArray::SmartPtr>& data, const int maxIterations )
{
  const size_t numberOfInputs = data.size();
  const size_t numberOfPixels = data[ 0 ]->GetDataSize();

  int numberOfClasses = 1;
  for ( size_t k = 0; k < numberOfInputs; ++k )
    {
    const Types::DataItemRange range = data[k]->GetRange();
    numberOfClasses = std::max( numberOfClasses, 1+static_cast<int>( range.m_UpperBound ) );
    }

  // allocate priors vector
  this->m_Priors.resize( numberOfClasses );

  // init priors
  size_t totalMass = 0;
  std::fill( this->m_Priors.begin(), this->m_Priors.end(), static_cast<RealValueType>( 0.0 ) );
  for ( size_t k = 0; k < numberOfInputs; ++k )
    {
    Types::DataItem lVal;
    for ( size_t n = 0; n < numberOfPixels; ++n )
      {
      if ( data[k]->Get( lVal, n ) )
	{
	this->m_Priors[static_cast<int>(lVal)]++;
	++totalMass;
	}
      }
    }
  for ( int l = 0; l < numberOfClasses; ++l )
    this->m_Priors[l] /= totalMass;

  // initialize result using simple voting.
  { LabelCombinationVoting voting( data ); this->m_Result = voting.GetResult(); } // use local scope to free voting object storage right away

  // allocate current and updated confusion matrix arrays
  this->m_Confusion.resize( numberOfInputs );
  this->m_ConfusionNew.resize( numberOfInputs );
  for ( size_t k = 0; k < numberOfInputs; ++k )
    {
    this->m_Confusion[k].Resize( 1+numberOfClasses, numberOfClasses );
    this->m_ConfusionNew[k].Resize( 1+numberOfClasses, numberOfClasses );
    }

  // initialize confusion matrices from voting result
  for ( size_t k = 0; k < numberOfInputs; ++k )
    {
    this->m_Confusion[k].SetAll( 0.0 );

    for ( size_t n = 0; n < numberOfPixels; ++n )
      {
      Types::DataItem lValue, vValue;
      if ( data[k]->Get( lValue, n ) )
	{
	if ( this->m_Result->Get( vValue, n ) && (vValue >= 0) )
	  ++(this->m_Confusion[k][static_cast<int>(lValue)][static_cast<int>(vValue)]);
	}
      }
    }
  
  // normalize matrix rows to unit probability sum
  for ( size_t k = 0; k < numberOfInputs; ++k )
    {
    for ( int inLabel = 0; inLabel <= numberOfClasses; ++inLabel )
      {
      // compute sum over all output labels for given input label
      float sum = 0;
      for ( int outLabel = 0; outLabel < numberOfClasses; ++outLabel )
	{
	sum += this->m_Confusion[k][inLabel][outLabel];
	}
      
      // make sure that this input label did in fact show up in the input!!
      if ( sum > 0 )
	{
	// normalize
	for ( int outLabel = 0; outLabel < numberOfClasses; ++outLabel )
	  {
	  this->m_Confusion[k][inLabel][outLabel] /= sum;
	  }
	}
      }
    }
  
  // allocate array for pixel class weights
  std::vector<float> W( numberOfClasses );

  Progress::Begin( 0, maxIterations, 1, "Multi-label STAPLE" );
  
  // main EM loop
  for ( int it = 0; it < maxIterations; ++it )
    {
    Progress::SetProgress( it );

    // reset updated confusion matrices.
    for ( size_t k = 0; k < numberOfInputs; ++k )
      {
      this->m_ConfusionNew[k].SetAll( 0.0 );
      }

    for ( size_t n = 0; n < numberOfPixels; ++n )
      {
      // the following is the E step
      for ( int ci = 0; ci < numberOfClasses; ++ci )
	W[ci] = this->m_Priors[ci];
      
      for ( size_t k = 0; k < numberOfInputs; ++k )
	{
	Types::DataItem lValue;
	if ( data[k]->Get( lValue, n ) )
	  {
	  for ( int ci = 0; ci < numberOfClasses; ++ci )
	    {
	    W[ci] *= this->m_Confusion[k][static_cast<int>(lValue)][ci];
	    }
	  }
	}
      
      // the following is the M step
      float sumW = W[0];
      for ( int ci = 1; ci < numberOfClasses; ++ci )
	sumW += W[ci];
      
      if ( sumW )
	{
	for ( int ci = 0; ci < numberOfClasses; ++ci )
	  W[ci] /= sumW;
	}
      
      for ( size_t k = 0; k < numberOfInputs; ++k )
	{
	Types::DataItem lValue;
	if ( data[k]->Get( lValue, n ) )
	  {
	  for ( int ci = 0; ci < numberOfClasses; ++ci )
	    {
	    this->m_ConfusionNew[k][static_cast<int>(lValue)][ci] += W[ci];	    
	    }
	  }
	}
      }

    // Normalize matrix elements of each of the updated confusion matrices
    // with sum over all expert decisions.
    for ( size_t k = 0; k < numberOfInputs; ++k )
      {
      // compute sum over all output classifications
      for ( int ci = 0; ci < numberOfClasses; ++ci ) 
	{
	float sumW = this->m_ConfusionNew[k][0][ci]; 
	for ( int j = 1; j <= numberOfClasses; ++j )
	  sumW += this->m_ConfusionNew[k][j][ci];
	
	// normalize with for each class ci
	if ( sumW )
	  {
	  for ( int j = 0; j <= numberOfClasses; ++j )
	    this->m_ConfusionNew[k][j][ci] /= sumW;
	  }
	}
      }
  
    // now we're applying the update to the confusion matrices and compute the
    // maximum parameter change in the process.
    for ( size_t k = 0; k < numberOfInputs; ++k )
      for ( int j = 0; j <= numberOfClasses; ++j )
	for ( int ci = 0; ci < numberOfClasses; ++ci )
	  {
	  this->m_Confusion[k][j][ci] = this->m_ConfusionNew[k][j][ci];
	  }    
    } // main EM loop

  // assemble output
  for ( size_t n = 0; n < numberOfPixels; ++n )
    {
    // basically, we'll repeat the E step from above
    for ( int ci = 0; ci < numberOfClasses; ++ci )
      W[ci] = this->m_Priors[ci];
    
    for ( size_t k = 0; k < numberOfInputs; ++k )
      {
      Types::DataItem lValue;
      if ( data[k]->Get( lValue, n ) )
	{
	for ( int ci = 0; ci < numberOfClasses; ++ci )
	  {
	  W[ci] *= this->m_Confusion[k][static_cast<int>(lValue)][ci];
	  }
	}
      }
    
    // now determine the label with the maximum W
    int winningLabel = -1;
    float winningLabelW = 0;
    for ( int ci = 0; ci < numberOfClasses; ++ci )
      {
      if ( W[ci] > winningLabelW )
	{
	winningLabelW = W[ci];
	winningLabel = ci;
	}
      else
	if ( ! (W[ci] < winningLabelW ) )
	  {
	  winningLabel = -1;
	  }
      }
    
    this->m_Result->Set( winningLabel, n );
    }

  Progress::Done();
}

} // namespace cmtk
