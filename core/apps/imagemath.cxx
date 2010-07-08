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

#include <cmtkconfig.h>

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>

#include <cmtkUniformVolume.h>
#include <cmtkFileFormat.h>
#include <cmtkVolumeIO.h>

#include <cmtkLabelCombinationVoting.h>
#include <cmtkLabelCombinationSTAPLE.h>
#include <cmtkLabelCombinationMultiClassSTAPLE.h>

#include <math.h>
#include <cmtkMathFunctionWrappers.h>
#include <cmtkTypedArrayFunctionHistogramMatching.h>

#include <queue>
#include <vector>

namespace cmtk
{

/// Square function.
double Square( const double x )
{
  return x*x; 
}

/// Logit function.
double Logit( const double x )
{
  return log(x / (1.0-x)); 
}

/// Logistic function.
double Logistic( const double x )
{
  return 1.0/(1.0+exp(-x));
}

} // namespace cmtk

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace imagemath
{
#endif

/// Flag for verbose operation.
bool Verbose = false;

/// The output data type of image operations.
cmtk::ScalarDataType ResultType = cmtk::TYPE_FLOAT;

/// Value used for padding input data, if PaddingFlag is true.
cmtk::Types::DataItem PaddingValue;

/// If this flag is set, PaddingValue defines the padding value for images read from files.
bool PaddingFlag = false;

void
CallbackSetPaddingValue( const double paddingValue )
{
  PaddingValue = paddingValue;
  PaddingFlag = true;
}

void
CallbackUnsetPaddingValue()
{
  PaddingFlag = false;
}

/// Type for "stack" of images.
typedef std::deque<cmtk::UniformVolume::SmartPtr> ImageStackType;

/// Operating stack of images.
ImageStackType ImageStack;

/// If this flag is set, the next single-image operation will be applied to all images on the stack.
bool ApplyNextToAll = false;

bool
CheckStackOneImage( const char* function )
{
  if ( ImageStack.empty() )
    {
    cmtk::StdErr << "ERROR: stack is empty in function '" << function << "'\n";
    return false;
    }
  return true;
}

bool
CheckStackTwoMatchingImages( const char* function )
{
  if ( ImageStack.size() < 2 )
    {
    cmtk::StdErr << "ERROR: at least two images are required to perform operation '" << function << "'\n";
    return false;
    }

  cmtk::UniformVolume::SmartPtr top = ImageStack.front();
  ImageStack.pop_front();
  const bool pixelCountsMatch = (top->GetNumberOfPixels() == ImageStack.front()->GetNumberOfPixels() );
  ImageStack.push_front( top );

  if ( !pixelCountsMatch )
    {
    cmtk::StdErr << "ERROR: two top images do not have equal pixel counts as required for operation '" << function << "'\n";
    return false;
    }
  
  return true;
}

void
CallbackIn( const char** argv, int& argsUsed )
{
  argsUsed = 0;
  while ( argv[argsUsed][0] != '-' )
    {
    cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( argv[argsUsed], Verbose ) );
    if ( !volume || !volume->GetData() )
      {
      cmtk::StdErr << "ERROR: could not read input image " << argv[argsUsed] << "\n";
      exit( 1 );
      }

    if ( PaddingFlag )
      volume->GetData()->SetPaddingValue( PaddingValue );

    ImageStack.push_front( volume );
    ++argsUsed;
    }
}

void
CallbackOut( const char* argv )
{
  if ( CheckStackOneImage( "Out" ) )
    {
    cmtk::VolumeIO::Write( *(ImageStack.front()), argv, Verbose );
    }
}

void
CallbackPop()
{
  if ( CheckStackOneImage( "Pop" ) )
    {
    ImageStack.pop_front();
    }
}

void
CallbackDup()
{
  if ( CheckStackOneImage( "Dup" ) )
    ImageStack.push_front( cmtk::UniformVolume::SmartPtr( ImageStack.front()->Clone() ) );
}
    
void
CallbackAll()
{
  ApplyNextToAll = true;
}
    
void
CallbackFill( const double value)
{
  if ( CheckStackOneImage( "Fill" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->Fill( value );

      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackAbs()
{
  if ( CheckStackOneImage( "Abs" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Abs );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackTrunc()
{
  if ( CheckStackOneImage( "Trunc" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Trunc );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackLog()
{
  if ( CheckStackOneImage( "Log" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Log );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackLogit()
{
  if ( CheckStackOneImage( "Logit" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Logit );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackLogistic()
{
  if ( CheckStackOneImage( "Logistic" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Logistic );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackExp()
{
  if ( CheckStackOneImage( "Exp" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Exp );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}

void
CallbackSqr()
{
  if ( CheckStackOneImage( "Sqr" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Square );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}

void
CallbackSqrt()
{
  if ( CheckStackOneImage( "Sqrt" ) )
    {
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
      (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Sqrt );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
        
void
CallbackThreshBelow( const char* argv )
{
  if ( CheckStackOneImage( "ThreshBelow" ) )
    {
    const float threshold = atof( argv );
    
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      cmtk::TypedArray::SmartPtr data = (*it)->GetData();
      
      const cmtk::Types::DataItemRange range = data->GetRange();
      data->Threshold( cmtk::Types::DataItemRange( threshold, range.m_UpperBound ) );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackThreshAbove( const char* argv )
{
  if ( CheckStackOneImage( "ThreshAbove" ) )
    {
    const float threshold = atof( argv );
    
    ImageStackType::iterator it = ImageStack.begin();
    while ( it != ImageStack.end() )
      {
      cmtk::TypedArray::SmartPtr data = (*it)->GetData();
      
      const cmtk::Types::DataItemRange range = data->GetRange();
      data->Threshold( cmtk::Types::DataItemRange( range.m_LowerBound, threshold ) );
      
      if ( ApplyNextToAll )
	++it;
      else
	it = ImageStack.end();
      }
    ApplyNextToAll = false;
    }
}
    
void
CallbackScalarMul( const double c )
{
  if ( ! CheckStackOneImage( "ScalarMul" ) )
    return;
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    cmtk::UniformVolume& p = **it;
    const size_t numberOfPixels = p.GetNumberOfPixels();
    
    cmtk::TypedArray::SmartPtr mul( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
    
#pragma omp parallel for
    for ( size_t i = 0; i < numberOfPixels; ++i )
      {
      cmtk::Types::DataItem pv;
      if ( p.GetDataAt( pv, i ) )
	{
	mul->Set( c * pv, i );
	}
      else
	{
	mul->SetPaddingAt( i );
	}
      }
    
    p.SetData( mul );
    
    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}

void
CallbackScalarAdd( const double c )
{
  if ( ! CheckStackOneImage( "ScalarAdd" ) )
    return;
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    cmtk::UniformVolume& p = **it;
    const size_t numberOfPixels = p.GetNumberOfPixels();
    
    cmtk::TypedArray::SmartPtr add( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
#pragma omp parallel for
    for ( size_t i = 0; i < numberOfPixels; ++i )
      {
      cmtk::Types::DataItem pv;
      if ( p.GetDataAt( pv, i ) )
	{
	add->Set( c + pv, i );
	}
      else
	{
	add->SetPaddingAt( i );
	}
      }
    
    p.SetData( add );
    
    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}

void
CallbackScalarXor( const long int c )
{
  if ( ! CheckStackOneImage( "ScalarXor" ) )
    return;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();

  const size_t numberOfPixels = p->GetNumberOfPixels();

  cmtk::TypedArray::SmartPtr out( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv;
    if ( p->GetDataAt( pv, i ) )
      {
      const long int iv = static_cast<long int>( pv );
      out->Set( iv ^ c, i );
      }
    else
      {
      out->SetPaddingAt( i );
      }
    }
  
  p->SetData( out );
  ImageStack.push_front( p );
}

void
CallbackScalarAnd( const long int c )
{
  if ( ! CheckStackOneImage( "ScalarAnd" ) )
    return;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();

  const size_t numberOfPixels = p->GetNumberOfPixels();

  cmtk::TypedArray::SmartPtr out( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv;
    if ( p->GetDataAt( pv, i ) )
      {
      const long int iv = static_cast<long int>( pv );
      out->Set( iv & c, i );
      }
    else
      {
      out->SetPaddingAt( i );
      }
    }
  
  p->SetData( out );
  ImageStack.push_front( p );
}

void
CallbackOneOver()
{
  if ( ! CheckStackOneImage( "OneOver" ) )
    return;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();

  const size_t numberOfPixels = p->GetNumberOfPixels();

  cmtk::TypedArray::SmartPtr inv( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv;
    if ( p->GetDataAt( pv, i ) )
      {
      inv->Set( 1.0 / pv, i );
      }
    else
      {
      inv->SetPaddingAt( i );
      }
    }
  
  p->SetData( inv );
  ImageStack.push_front( p );
}

void
CallbackAdd()
{
  if ( ! CheckStackTwoMatchingImages( "Add" ) )
    return;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();
  cmtk::UniformVolume::SmartPtr q = ImageStack.front();
  ImageStack.pop_front();
  
  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr add( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv, qv;
    if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) )
      {
      add->Set( pv+qv, i );
      }
    else
      {
      add->SetPaddingAt( i );
      }
    }
  
  p->SetData( add );
  ImageStack.push_front( p );
}

void
CallbackMul()
{
  if ( ! CheckStackTwoMatchingImages( "Mul" ) )
    return;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();
  cmtk::UniformVolume::SmartPtr q = ImageStack.front();
  ImageStack.pop_front();
  
  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr mul( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv, qv;
    if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) )
      {
      mul->Set( pv*qv, i );
      }
    else
      {
      mul->SetPaddingAt( i );
      }
    }
  
  p->SetData( mul );
  ImageStack.push_front( p );
}

void
CallbackDiv()
{
  if ( ! CheckStackTwoMatchingImages( "Div" ) )
    return;

  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();
  cmtk::UniformVolume::SmartPtr q = ImageStack.front();
  ImageStack.pop_front();
  
  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr div( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv, qv;
    if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) && (qv != 0) )
      {
      div->Set( pv/qv, i );
      }
    else
      {
      div->SetPaddingAt( i );
      }
    }
  
  p->SetData( div );
  ImageStack.push_front( p );
}

void
CallbackAtan2()
{
  if ( ! CheckStackTwoMatchingImages( "Atan2" ) )
    return;

  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();
  cmtk::UniformVolume::SmartPtr q = ImageStack.front();
  ImageStack.pop_front();
  
  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv, qv;
    if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) && (qv != 0) )
      {
      result->Set( atan2( qv, pv ), i );
      }
    else
      {
      result->SetPaddingAt( i );
      }
    }
  
  p->SetData( result );
  ImageStack.push_front( p );
}

void
CallbackMatchHistograms()
{
  if ( ImageStack.size() < 2 )
    {
    cmtk::StdErr << "ERROR: need at least two images on stack for histogram intensity matching\n";
    return;
    }
  
  cmtk::UniformVolume::SmartPtr ref = ImageStack.front();
  ImageStack.pop_front();
  ImageStack.front()->GetData()->ApplyFunctionObject( cmtk::TypedArrayFunctionHistogramMatching( *(ImageStack.front()->GetData()), *(ref->GetData()) ) );
}

void
CallbackSum()
{
  while ( ImageStack.size() > 1 )
    {
    CallbackAdd();
    }
}

void
CallbackAverage()
{
  const size_t nimages = ImageStack.size();
  CallbackSum();
  ImageStack.front()->GetData()->Rescale( 1.0 / nimages );
}

void
CallbackVoteCombination()
{
  if ( ! CheckStackTwoMatchingImages( "Vote" ) )
    return;

  cmtk::UniformVolume::SmartPtr grid = ImageStack.front();
  
  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "Vote" ) )
	return;
    
    dataPtrs.push_back( ImageStack.front()->GetData() );
    ImageStack.pop_front();
    }

  cmtk::LabelCombinationVoting voting( dataPtrs );
  grid->SetData( voting.GetResult() );
  ImageStack.push_front( grid );
}

void
CallbackSTAPLE( const long int maxIterations )
{
  if ( ! CheckStackTwoMatchingImages( "STAPLE" ) )
    return;
  
  cmtk::UniformVolume::SmartPtr imageGrid = ImageStack.front();

  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "STAPLE" ) )
	return;
    
    dataPtrs.push_back( ImageStack.front()->GetData() );
    ImageStack.pop_front();
    }

  cmtk::LabelCombinationSTAPLE staple( dataPtrs, maxIterations, ResultType );
  imageGrid->SetData( staple.GetResult() );
  ImageStack.push_front( imageGrid );

  if ( Verbose )
    {
    cmtk::StdErr.printf( "p  " );
    for ( size_t i = 0; i < dataPtrs.size(); ++i )
      {
      cmtk::StdErr.printf( "%.3f ", staple.GetPValue( i ) );
      }
    cmtk::StdErr.printf( "\nq  " );
    for ( size_t i = 0; i < dataPtrs.size(); ++i )
      {
      cmtk::StdErr.printf( "%.3f ", staple.GetQValue( i ) );
      }
    cmtk::StdErr.printf( "\n" );
    }
}

void
CallbackMultiClassSTAPLE( const long int maxIterations )
{
  if ( ! CheckStackTwoMatchingImages( "MultiClassSTAPLE" ) )
    return;
  
  cmtk::UniformVolume::SmartPtr imageGrid = ImageStack.front();

  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "MultiClassSTAPLE" ) )
	return;
    
    dataPtrs.push_back( ImageStack.front()->GetData() );
    ImageStack.pop_front();
    }

  cmtk::LabelCombinationMultiClassSTAPLE mstaple( dataPtrs, maxIterations );
  imageGrid->SetData( mstaple.GetResult() );
  ImageStack.push_front( imageGrid );
}

void
CallbackStackEntropyLabels()
{
  if ( ! CheckStackTwoMatchingImages( "StackEntropyLabels" ) )
    return;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "StackEntropyLabels" ) )
	return;
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr entropyArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );

#pragma omp parallel for  
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    std::map<int,unsigned int> labelCount;

    size_t totalCount = 0;
    for ( size_t curVol = 0; curVol < volPtrs.size(); ++curVol )
      {
      cmtk::Types::DataItem v;
      if ( volPtrs[ curVol ]->GetDataAt( v, i ) ) 
        {
        ++labelCount[ static_cast<int>( v ) ];
	++totalCount;
        }
      }

    if ( totalCount )
      {
      const double factor = 1.0 / totalCount;
      
      double entropy = 0;
      for ( std::map<int,unsigned int>::const_iterator it = labelCount.begin(); it != labelCount.end(); ++it )
	{
	const double p = factor * it->second;
	entropy += p * log( p );
	}
      entropyArray->Set( -entropy / log( 2.0 ), i );
      }
    else
      entropyArray->SetPaddingAt( i );       
    }

  volPtrs[0]->SetData( entropyArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackMaxIndex()
{
  if ( ! CheckStackTwoMatchingImages( "MaxIndex" ) )
    return;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs( ImageStack.size() );

  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "MaxIndex" ) )
	return;
    
    volPtrs[ImageStack.size()-1] = ImageStack.front();
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr maxArray( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, numberOfPixels ) );

#pragma omp parallel for  
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    float maxValue = 0;
    short maxIndex = -2;
    cmtk::Types::DataItem v;
    
    for ( size_t curVol = 0; curVol < volPtrs.size(); ++curVol )
      {
      if ( volPtrs[ curVol ]->GetDataAt( v, i ) ) 
        {
	if ( maxIndex < -1 )
	  {
	  maxValue = v;
	  maxIndex = 0;
	  }
	else
	  {
	  if ( v > maxValue )
	    {
	    maxValue = v;
	    maxIndex = curVol;
	    }
	  else
	    if ( v == maxValue )
	      {
	      maxIndex = -1;
	      }
	  }
        }
      }
    
    maxArray->Set( maxIndex, i ); 
    }
  
  volPtrs[0]->SetData( maxArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackMaxValue()
{
  if ( ! CheckStackTwoMatchingImages( "MaxValue" ) )
    return;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "MaxValue" ) )
	return;
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr maxArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );

#pragma omp parallel for  
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem maxValue = 0;
    bool maxValueValid = false;
    cmtk::Types::DataItem v;

    for ( size_t curVol = 0; curVol < volPtrs.size(); ++curVol )
      {
      if ( volPtrs[ curVol ]->GetDataAt( v, i ) ) 
        {
	if ( maxValueValid )
	  {
	  maxValue = std::max( maxValue, v );
	  }
	else
	  {
	  maxValueValid = true;
	  maxValue = v;
	  }
        }
      }
    
    
    maxArray->Set( maxValue, i ); 
    }
  
  volPtrs[0]->SetData( maxArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackMinValue()
{
  if ( ! CheckStackTwoMatchingImages( "MinValue" ) )
    return;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "MinValue" ) )
	return;
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr minArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );

#pragma omp parallel for  
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem minValue = 0;
    bool minValueValid = false;
    cmtk::Types::DataItem v;

    for ( size_t curVol = 0; curVol < volPtrs.size(); ++curVol )
      {
      if ( volPtrs[ curVol ]->GetDataAt( v, i ) ) 
        {
	if ( minValueValid )
	  {
	  minValue = std::min( minValue, v );
	  }
	else
	  {
	  minValueValid = true;
	  minValue = v;
	  }
        }
      }
    
    
    minArray->Set( minValue, i ); 
    }
  
  volPtrs[0]->SetData( minArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackContractLabels()
{
  if ( ! CheckStackTwoMatchingImages( "ContractLabels" ) )
    return;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "ContractLabels" ) )
	return;
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr outArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );

#pragma omp parallel for  
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem v = 0;

    for ( size_t curVol = 0; curVol < volPtrs.size(); ++curVol )
      {
      if ( volPtrs[ curVol ]->GetDataAt( v, i ) && v ) 
	break;
      }
        
    outArray->Set( v, i );
    }
  
  volPtrs[0]->SetData( outArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackCombinePCA()
{
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "CombinePCA" ) )
	return;
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[0]->GetNumberOfPixels();
  const size_t numberOfImages = volPtrs.size();
  
  cmtk::Vector<cmtk::Types::DataItem> meanVector( numberOfImages );
  for ( size_t image = 0; image < numberOfImages; ++image )
    {
    cmtk::Types::DataItem mean = 0;

    for ( size_t pixel = 0; pixel < numberOfPixels; ++pixel ) 
      {      
      cmtk::Types::DataItem value;
      if ( volPtrs[image]->GetDataAt( value, pixel ) )
	{
	mean += value;
	}
      }
    
    meanVector[image] = mean / numberOfPixels;
    }
  
  cmtk::Matrix2D<cmtk::Types::DataItem> cc( numberOfImages, numberOfImages );

  for ( size_t imageY = 0; imageY < numberOfImages; ++imageY )
    {
    for ( size_t imageX = 0; imageX < numberOfImages; ++imageX )
      {
      // compute upper half of matrix, use symmetry to fill lower half
      if ( imageY <= imageX )
	{
	cmtk::Types::DataItem ccXY = 0;
	
	for ( size_t pixel = 0; pixel < numberOfPixels; ++pixel )
	  {
	  cmtk::Types::DataItem valueX, valueY;
	  if ( volPtrs[imageX]->GetDataAt( valueX, pixel ) && volPtrs[imageY]->GetDataAt( valueY, pixel ) )
	    {
	    ccXY += ( valueX - meanVector[imageX] ) * ( valueY - meanVector[imageY] );
	    }
	  }
	
	cc[imageX][imageY] = ccXY / numberOfImages;
	}
      else
	{
	cc[imageX][imageY] = cc[imageY][imageX];
	}      
      }
    }
  
  cmtk::Matrix2D<cmtk::Types::DataItem> eigensystem( numberOfImages, numberOfImages );
  std::vector<cmtk::Types::DataItem> eigenvalues( numberOfImages );
  
  cmtk::MathUtil::ComputeEigensystem( cc, eigensystem, eigenvalues );
  
  // find maximum eigenvalue
  int maxL = 0;
  cmtk::Types::DataItem maxLambda = fabs( eigenvalues[0] );
  for ( size_t l = 1; l < numberOfImages; ++l )
    {
    if ( fabs( eigenvalues[l] ) > maxLambda )
      {
      maxLambda = fabs( eigenvalues[l] );
      maxL = l;
      }
    }
  
  // get normalized pricipal component
  cmtk::Vector<cmtk::Types::DataItem> ev( numberOfImages );
  for ( size_t l = 0; l < numberOfImages; ++l )
    {
    ev[l] = eigensystem[l][maxL];
    }
  ev *= 1.0 / ev.EuclidNorm();
  
  if ( Verbose )
    {
    cmtk::StdErr << "Principal eigenvector (normalized):\n";
    for ( size_t l = 0; l < numberOfImages; ++l )
      {
      cmtk::StdErr.printf( "%f\t", ev[l] );
      }
    cmtk::StdErr << "\n";
    }
  
  // project all pixel vectors onto dominant component  
  cmtk::TypedArray::SmartPtr output( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
  cmtk::Vector<cmtk::Types::DataItem> v( numberOfImages );
  for ( size_t pixel = 0; pixel < numberOfPixels; ++pixel )
    {
    for ( size_t image = 0; image < numberOfImages; ++image )
      {
      if ( ! volPtrs[image]->GetDataAt( v[image], pixel ) )
	v[image] = 0;
      v[image] -= meanVector[image];	    
      }
    
    output->Set( v * ev, pixel );
    }
  volPtrs[0]->SetData( output );
  ImageStack.push_front( volPtrs[0] );
}

int
main( int argc, char *argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image operations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Perform operations on images using stack-based postfix notation.\n\n"
		       "Images can be read from files and pushed onto the stack. Images on the stack can be processed and combined via different operators. "
		       "Results of all operations are put back onto the stack, where they can be further processed or written back to image files." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Be verbose" );

    cl.BeginGroup( "Input/output", "Input/output operations" );
    cl.AddCallback( Key( "in" ), CallbackIn, "Read input image(s) to top of stack" );
    cl.AddCallback( Key( "out" ), CallbackOut, "Write output image from top of stack (but leave it on the stack)" );
    cl.AddCallback( Key( "set-padding-value" ), CallbackSetPaddingValue, "Set the value that is interpreted as padding value in subsequently read images." );
    cl.AddCallback( Key( "unset-padding" ), CallbackUnsetPaddingValue, "Disable padding. All values in subsequently read images will be interpreted as actual data." );
    cl.EndGroup();

    cl.BeginGroup( "Internal", "Internal settings" );
    cl.AddSwitch( Key( "float" ), &ResultType, cmtk::TYPE_FLOAT, "Use single precision for computations and results" );
    cl.AddSwitch( Key( "double" ), &ResultType, cmtk::TYPE_DOUBLE, "Use double precision for computations and results" );
    cl.EndGroup();

    cl.BeginGroup( "Stack", "Stack operations" );
    cl.AddCallback( Key( "pop" ), CallbackPop, "Pop (discard) top image from stack." );
    cl.AddCallback( Key( "dup" ), CallbackDup, "Duplicate image on top of the stack." );
    cl.AddCallback( Key( "all" ), CallbackAll, "Apply next single-image operation to all images on the stack." );
    cl.EndGroup();

    cl.BeginGroup( "Single image", "Single-image operators" );
    cl.AddCallback( Key( "fill" ), CallbackFill, "Fill top image with constant value (i.e., assign value to all pixels)" );
    cl.AddCallback( Key( "abs" ), CallbackAbs, "Apply abs() function to top image" );
    cl.AddCallback( Key( "log" ), CallbackLog, "Apply log() function to top image" );
    cl.AddCallback( Key( "logit" ), CallbackLogit, "Apply log(x/(1-x)) function to top image" );
    cl.AddCallback( Key( "logistic" ), CallbackLogistic, "Apply 1/(1+exp(-x)) function to top image" );
    cl.AddCallback( Key( "exp" ), CallbackExp, "Apply exp() function to top image" );
    cl.AddCallback( Key( "sqr" ), CallbackSqr, "Apply square operator to top image" );
    cl.AddCallback( Key( "sqrt" ), CallbackSqrt, "Apply square root operator to top image" );
    cl.AddCallback( Key( "trunc" ), CallbackTrunc, "Truncate all values in top image to integer" );
    cl.AddCallback( Key( "one-over" ), CallbackOneOver, "For each pixel, replace its value x with 1.0/x" );
    cl.AddCallback( Key( "scalar-mul" ), CallbackScalarMul, "Multiply top image with a scalar value" );
    cl.AddCallback( Key( "scalar-add" ), CallbackScalarAdd, "Add a scalar to each pixel of the top image" );
    cl.AddCallback( Key( "scalar-xor" ), CallbackScalarXor, "Bitwise exclusive-or between top level and given scalar value" );
    cl.AddCallback( Key( "scalar-and" ), CallbackScalarAnd, "Bitwise and operation between top level and given scalar value" );

    cl.AddCallback( Key( "threshold-below" ), CallbackThreshBelow, "Set values below given threshold to threshold." );
    cl.AddCallback( Key( "threshold-above" ), CallbackThreshAbove, "Set values above given threshold to threshold." );
    cl.EndGroup();

    cl.BeginGroup( "Two images", "Image pair operators" );
    cl.AddCallback( Key( "add" ), CallbackAdd, "Add top and second image, place result on stack" );
    cl.AddCallback( Key( "mul" ), CallbackMul, "Multiply top and second image, place result on stack" );
    cl.AddCallback( Key( "div" ), CallbackDiv, "Divide top image by second image, place result on stack" );
    cl.AddCallback( Key( "atan2" ), CallbackAtan2, "Compute atan2() function from tup two image pixel pairs, place result on stack" );    
    cl.AddCallback( Key( "match-histograms" ), CallbackMatchHistograms, "Scale intensities one image to match intensities another. The image pushed onto the stack last provides the reference intensity distribution, the preceding image will be modified. Both input images will be removed from the stack and the modified image will be pushed onto the stack." );
    cl.EndGroup();

    cl.BeginGroup( "Contract multiple images", "Operators that contract the entire stack into a single image" );
    cl.AddCallback( Key( "sum" ), CallbackSum, "Sum all images on stack, place result on stack" );
    cl.AddCallback( Key( "average" ), CallbackAverage, "Average all images on stack, place result on stack" );
    cl.AddCallback( Key( "combine-pca" ), CallbackCombinePCA, "Combine images using PCA by projecting onto direction of largest correlation" );
    cl.AddCallback( Key( "max-value" ), CallbackMaxValue, "For each pixel, compute maximum VALUE over all images, place result on stack" );
    cl.AddCallback( Key( "min-value" ), CallbackMinValue, "For each pixel, compute minimum VALUE over all images, place result on stack" );
    cl.AddCallback( Key( "max-index" ), CallbackMaxIndex, "For each pixel, compute INDEX of image with maximum value, place result on stack" );
    cl.AddCallback( Key( "contract-labels" ), CallbackContractLabels, "Contract multiple label maps into one by selecting the first (over all images on the stack) non-zero label at each pixel" );
    cl.EndGroup();

    cl.BeginGroup( "Contract multiple label images", "Operators that contract a stack of label images into a single label image" );
    cl.AddCallback( Key( "vote" ), CallbackVoteCombination, "Merge all images on stack with voting, place result on stack" );
    cl.AddCallback( Key( "staple" ), CallbackSTAPLE, "Combine binary maps on the stack using [arg] iterations of the STAPLE algorithm. "
		    "The result of this operation is the spatial map of 'weights' W, which are the probabilities of image foreground at each pixel. In 'verbose' "
		    "mode, estimated expert parameters p (sensitivity) and q (specificity) are also written to standard output." );
    cl.AddCallback( Key( "mstaple" ), CallbackMultiClassSTAPLE, "Combine multi-label maps on the stack using [arg] iterations of the multi-class STAPLE algorithm."
		    "The result of this operation is the combined maximum-likeliood multi-label map." );
    cl.AddCallback( Key( "stack-entropy-labels" ), CallbackStackEntropyLabels, "Compute stack entropy at each pixel from integer (label) input images" );
    cl.EndGroup();

    cl.Parse();

    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }
  
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace imagemath
} // namespace apps
} // namespace cmtk
#endif
