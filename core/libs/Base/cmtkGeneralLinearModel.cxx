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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include "cmtkGeneralLinearModel.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <System/cmtkProgress.h>

#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

#define TOL 1.0e-5
GeneralLinearModel::GeneralLinearModel
( const size_t nParameters, const size_t nData, const double* designMatrix ) :
  NParameters( nParameters ),
  NData( nData ),
  DesignMatrix( nData, nParameters, designMatrix ),
  Up( nParameters ),
  Vp( nParameters ),
  Wp( nParameters ),
  VariableMean( nParameters ),
  VariableSD( nParameters )
{
  U = new Matrix2D<double>( NData, NParameters );
  V = new Matrix2D<double>( NParameters, NParameters );
  W = new std::vector<double>( NParameters );

  double wmax, thresh;

  std::vector<double> data( this->NData );
  for ( size_t j=0; j < NParameters; ++j ) 
    {
    // set up data vector for parameter 'j'
    for ( size_t i=0; i < this->NData; ++i ) 
      {
      data[i] = DesignMatrix[i][j];
      (*U)[i][j] = DesignMatrix[i][j];
      }

    // compute variance
    this->VariableMean[j] = MathUtil::Mean<double>( data );
    this->VariableSD[j] = MathUtil::Variance<double>( data, this->VariableMean[j] );
    
    // convert variance to standard deviation
    this->VariableSD[j] = sqrt( this->VariableSD[j] );
    }
  
  // perform SVD of design matrix
  MathUtil::SVD( this->U, this->NData, this->NParameters, this->W, this->V );
 
  // prepare partial regressions, each with one of the parameters omitted
  for ( size_t p=0; p < this->NParameters; ++p) 
    {
    Up[p] = new Matrix2D<double>( NData, NParameters-1 );
    Vp[p] = new Matrix2D<double>( NParameters-1, NParameters-1 );
    Wp[p] = new std::vector<double>( NParameters-1 );
    
    // create partial design matrix, omitting parameter 'p'
    for ( size_t i=0; i < this->NData; ++i ) 
      {
      size_t jj = 0;
      for ( size_t j=0; j < this->NParameters; ++j ) 
	{
	if ( j != p ) 
	  {
	  (*(this->Up[p]))[i][jj] = DesignMatrix[i][j];
	  ++jj;
	  }
	}
      }
    
    MathUtil::SVD( this->Up[p], this->NData, this->NParameters-1, this->Wp[p], this->Vp[p] );
    }
  
  wmax=0.0;
  for ( size_t j=0;j<NParameters;j++)
    if ((*W)[j] > wmax) wmax=(*W)[j];
  thresh=TOL*wmax;
  for ( size_t j=0;j<NParameters;j++)
    if ((*W)[j] < thresh) (*W)[j]=0.0;
}


Matrix2D<double>*
GeneralLinearModel::GetCorrelationMatrix() const
{
  Matrix2D<double>* CC = new Matrix2D<double>( this->NParameters, this->NParameters );

  std::vector<double> pi( this->NData );
  std::vector<double> pj( this->NData );

  for ( size_t i = 0; i < this->NParameters; ++i )
    {
    for ( size_t n = 0; n < this->NData; ++n ) 
      {
      pi[n] = this->DesignMatrix[n][i];
      }
    
    for ( size_t j = 0; j < this->NParameters; ++j ) 
      {
      if ( i <= j )
	{
	for ( size_t n = 0; n < this->NData; ++n ) 
	  {
	  pj[n] = this->DesignMatrix[n][j];
	  }
	(*CC)[i][j] = MathUtil::Correlation( pi, pj );
	}
      else
	{
	(*CC)[i][j] = (*CC)[j][i];
	}
      }
    }
  
  return CC;
}

GeneralLinearModel::~GeneralLinearModel()
{
  for ( size_t p=0; p < this->NParameters; ++p) 
    {
    delete this->Wp[p];
    delete this->Vp[p];
    delete this->Up[p];
    }
  delete this->W; 
  delete this->V; 
  delete this->U; 
}

void
GeneralLinearModel::InitResults( const size_t nPixels )
{
  Model.clear();
  TStat.clear();
  for (size_t p = 0; (p<this->NParameters); ++p ) 
    {
    TypedArray::SmartPtr nextModel( TypedArray::Create( TYPE_FLOAT, nPixels ) );
    Model.push_back( nextModel );
    
    TypedArray::SmartPtr nextTStat( TypedArray::Create( TYPE_FLOAT, nPixels ) );
    TStat.push_back( nextTStat );
    }
  
  FStat = TypedArray::SmartPtr( TypedArray::Create( TYPE_FLOAT, nPixels ) );
}

void
GeneralLinearModel::FitModel
( std::vector<TypedArray::SmartPtr>& y, const bool normalizeParameters )
{
  assert( y.size() == this->NData );
  const size_t nPixels = y[0]->GetDataSize();
  this->InitResults( nPixels );
  
  std::vector<double> lm_params( this->NParameters );
  std::vector<double> b( this->NData );
  std::vector<double> valueYhat( this->NData );
  
  // number of degrees of freedom for t-statistics
  // note: we omit "-1" because the constant in our model is either suppressed or an explicit parameter
  const int df = this->NData - this->NParameters;

  const size_t pixelUpdateIncrement = 10000;
  Progress::Begin( 0, nPixels, pixelUpdateIncrement, "Linear model fitting" );

  for ( size_t n = 0; n < nPixels; ++n ) 
    {
    if ( ! (n % pixelUpdateIncrement) )
      if ( Progress::SetProgress( n ) != Progress::OK ) break;
    
    bool missing = false;
    Types::DataItem value;
    for (size_t i = 0; (i<this->NData) && !missing; i++) 
      if ( y[i]->Get( value, n ) && finite( value ) )
	b[i] = value;
      else
	missing = true;

    if ( missing )
      {
      for (size_t p = 0; (p<this->NParameters); ++p ) 
	{
	this->Model[p]->SetPaddingAt( n );
	this->TStat[p]->SetPaddingAt( n );
	}
      } 
    else 
      {
      // use SVD of design matrix to compute model parameters lm_params[] from data b[]
      MathUtil::SVDLinearRegression( this->U, this->NData, this->NParameters, this->W, this->V, &b[0], lm_params );

      // compute variance of data
      double varY, avgY;
      avgY = MathUtil::Mean<double>( this->NData, &b[0] );
      varY = MathUtil::Variance<double>( this->NData, &b[0], avgY );
      
      // copy model parameters into output
      for (size_t p = 0; (p<this->NParameters); ++p ) 
	{
	value = lm_params[p];
	if ( normalizeParameters )
	  // Cohen & Cohen, Eq. (3.5.2)
//	  Model[p]->Set( lm_params[p] * this->GetNormFactor( p ) / sqrt( varY ), n );
	  value *= this->GetNormFactor( p );

	if ( finite( value ) )
	  this->Model[p]->Set( value, n );
	else
	  this->Model[p]->SetPaddingAt( n );
	}
      
      // compute variance of approximated data using entire model
      double varYhat, avgYhat;
      for (size_t i = 0; i<this->NData; i++) 
	{ 
	valueYhat[i] = 0.0;
	for (size_t pi = 0; (pi<this->NParameters); ++pi )
	  valueYhat[i] += lm_params[pi] * this->DesignMatrix[i][pi];
	}
      avgYhat = MathUtil::Mean<double>( this->NData, &valueYhat[0] );
      varYhat = MathUtil::Variance<double>( this->NData, &valueYhat[0], avgYhat );
      
      // compute multiple R square
      const double R2 = varYhat / varY;
      this->FStat->Set( (R2*df) / ((1-R2)*this->NParameters), n );
      
      std::vector<double> lm_params_P( this->NParameters-1 );
      std::vector<double> valueYhatp( this->NData );
      
      // for each parameter, evaluate R^2_i for model without parameter Xi
      for (size_t p = 0; p < this->NParameters; ++p ) 
	{
	// exclude constant parameter
//	if ( this->VariableSD[p] > 0 )
	  {
//	  // use SVD of partial design matrix to compute partial regression
          MathUtil::SVDLinearRegression( this->Up[p], this->NData, this->NParameters-1, this->Wp[p], this->Vp[p], &b[0], lm_params_P );

	  // compute variance of data
	  for (size_t i = 0; i < this->NData; i++) 
	    { 
	    valueYhatp[i] = 0.0;
	    size_t pip = 0;
	    for (size_t pi = 0; pi < this->NParameters; ++pi ) 
	      {
	      if ( p != pi ) 
		{
		valueYhatp[i] += lm_params_P[pip] * this->DesignMatrix[i][pi];
		++pip;
		}
	      }
	    }
	  
	  double varYhatp, avgYhatp;
          avgYhatp = MathUtil::Mean<double>( valueYhatp );
          varYhatp = MathUtil::Variance<double>( valueYhatp, avgYhatp );
	  
	  // copmpute R^2_p
	  const double R2p = varYhatp / varY;
//	  assert( (R2p >= 0) && (R2p < 1) );
	  // compute sr_p^2 from R^2 and R^2_p
	  const double srp = sqrt( R2 - R2p );
	  //	      assert( (sr2p >= 0) && (sr2p <= 1 ) );
	  // compute T-statistics
	  double tStat = static_cast<double>( srp * sqrt( df / (1.0-R2) ) );
	  // export T-statistics (set to zero if NAN)
	  if ( ! MathUtil::IsFinite( tStat ) ) 
	    tStat = 0;
	  this->TStat[p]->Set( tStat, n ); 
	  }
	}
      }
    }
  
  Progress::Done();
}

} // namespace cmtk
