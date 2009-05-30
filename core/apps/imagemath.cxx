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
#include <stack>

#ifdef _MSC_VER
double trunc( const double x )
{
	return static_cast<double>( static_cast<long int>( x ) );
}
#endif

namespace cmtk
{
double Square( const double x )
{
  return x*x; 
}
}

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace imagemath
{
#endif
bool Verbose = false;

cmtk::ScalarDataType ResultType = cmtk::TYPE_FLOAT;

std::stack<cmtk::UniformVolume::SmartPtr> ImageStack;

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
    cmtk::StdErr << "ERROR: two images are required to perform operation '" << function << "'\n";
    return false;
    }

  cmtk::UniformVolume::SmartPtr top = ImageStack.top();
  ImageStack.pop();
  const bool pixelCountsMatch = (top->GetNumberOfPixels() == ImageStack.top()->GetNumberOfPixels() );
  ImageStack.push( top );

  if ( !pixelCountsMatch )
    {
    cmtk::StdErr << "ERROR: two top images do not have equal pixel counts as required for operation '" << function << "'\n";
    return false;
    }
  
  return true;
}

const char*
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

    ImageStack.push( volume );
    ++argsUsed;
    }
  return NULL;
}

const char*
CallbackOut( const char* argv )
{
  if ( CheckStackOneImage( "Out" ) )
    {
    cmtk::VolumeIO::Write( ImageStack.top(), argv, Verbose );
    }
  
  return NULL;
}

const char*
CallbackPop()
{
  if ( CheckStackOneImage( "Pop" ) )
    {
    ImageStack.pop();
    }

  return NULL;
}

const char*
CallbackDup()
{
  if ( CheckStackOneImage( "Dup" ) )
    ImageStack.push( cmtk::UniformVolume::SmartPtr( ImageStack.top()->Clone() ) );
  
  return NULL;
}
    
const char*
CallbackAbs()
{
  if ( CheckStackOneImage( "Abs" ) )
    {
    ImageStack.top()->SetData( cmtk::TypedArray::SmartPtr( ImageStack.top()->GetData()->Convert( ResultType ) ) );
    ImageStack.top()->GetData()->ApplyFunctionDouble( fabs );
    }

  return NULL;
}
    
const char*
CallbackTrunc()
{
  if ( CheckStackOneImage( "Trunc" ) )
    {
    ImageStack.top()->SetData( cmtk::TypedArray::SmartPtr( ImageStack.top()->GetData()->Convert( ResultType ) ) );
	ImageStack.top()->GetData()->ApplyFunctionDouble( trunc );
    }

  return NULL;
}
    
const char*
CallbackLog()
{
  if ( CheckStackOneImage( "Log" ) )
    {
    ImageStack.top()->SetData( cmtk::TypedArray::SmartPtr( ImageStack.top()->GetData()->Convert( ResultType ) ) );
    ImageStack.top()->GetData()->ApplyFunctionDouble( log );
    }

  return NULL;
}
    
const char*
CallbackExp()
{
  if ( CheckStackOneImage( "Exp" ) )
    {
    ImageStack.top()->SetData( cmtk::TypedArray::SmartPtr( ImageStack.top()->GetData()->Convert( ResultType ) ) );
    ImageStack.top()->GetData()->ApplyFunctionDouble( exp );
    }

  return NULL;
}

const char*
CallbackSqr()
{
  if ( CheckStackOneImage( "Sqr" ) )
    {
    ImageStack.top()->SetData( cmtk::TypedArray::SmartPtr( ImageStack.top()->GetData()->Convert( ResultType ) ) );
    ImageStack.top()->GetData()->ApplyFunctionDouble( cmtk::Square );
    }

  return NULL;
}

const char*
CallbackSqrt()
{
  if ( CheckStackOneImage( "Sqrt" ) )
    {
    ImageStack.top()->SetData( cmtk::TypedArray::SmartPtr( ImageStack.top()->GetData()->Convert( ResultType ) ) );
    ImageStack.top()->GetData()->ApplyFunctionDouble( sqrt );
    }

  return NULL;
}
        
const char*
CallbackThreshMin( const char* argv )
{
  if ( CheckStackOneImage( "ThreshMin" ) )
    {
    const float threshold = atof( argv );
    
    cmtk::TypedArray::SmartPtr data = ImageStack.top()->GetData();
    
    cmtk::Types::DataItem min, max;
    data->GetRange( min, max );
    data->Threshold( threshold, max );
    }

  return NULL;
}
    
const char*
CallbackScalarMul( const char* argv )
{
  if ( ! CheckStackOneImage( "ScalarMul" ) )
    return NULL;
  
  const float c = atof( argv );
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.top();
  ImageStack.pop();

  const size_t numberOfPixels = p->GetNumberOfPixels();

  cmtk::TypedArray::SmartPtr mul( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv;
    if ( p->GetDataAt( pv, i ) )
      {
      mul->Set( c * pv, i );
      }
    else
      {
      mul->SetPaddingAt( i );
      }
    }
  
  p->SetData( mul );
  ImageStack.push( p );

  return NULL;
}

const char*
CallbackScalarAdd( const char* argv )
{
  if ( ! CheckStackOneImage( "ScalarAdd" ) )
    return NULL;
  
  const float c = atof( argv );
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.top();
  ImageStack.pop();

  const size_t numberOfPixels = p->GetNumberOfPixels();

  cmtk::TypedArray::SmartPtr add( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  
#pragma omp parallel for
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    cmtk::Types::DataItem pv;
    if ( p->GetDataAt( pv, i ) )
      {
      add->Set( c + pv, i );
      }
    else
      {
      add->SetPaddingAt( i );
      }
    }
  
  p->SetData( add );
  ImageStack.push( p );

  return NULL;
}

const char*
CallbackOneOver()
{
  if ( ! CheckStackOneImage( "OneOver" ) )
    return NULL;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.top();
  ImageStack.pop();

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
  ImageStack.push( p );

  return NULL;
}

const char*
CallbackAdd()
{
  if ( ! CheckStackTwoMatchingImages( "Add" ) )
    return NULL;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.top();
  ImageStack.pop();
  cmtk::UniformVolume::SmartPtr q = ImageStack.top();
  ImageStack.pop();
  
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
  ImageStack.push( p );

  return NULL;
}

const char*
CallbackMul()
{
  if ( ! CheckStackTwoMatchingImages( "Mul" ) )
    return NULL;
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.top();
  ImageStack.pop();
  cmtk::UniformVolume::SmartPtr q = ImageStack.top();
  ImageStack.pop();
  
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
  ImageStack.push( p );

  return NULL;
}

const char*
CallbackDiv()
{
  if ( ! CheckStackTwoMatchingImages( "Div" ) )
    return NULL;

  cmtk::UniformVolume::SmartPtr p = ImageStack.top();
  ImageStack.pop();
  cmtk::UniformVolume::SmartPtr q = ImageStack.top();
  ImageStack.pop();
  
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
  ImageStack.push( p );

  return NULL;
}

const char*
CallbackAtan2()
{
  if ( ! CheckStackTwoMatchingImages( "Atan2" ) )
    return NULL;

  cmtk::UniformVolume::SmartPtr p = ImageStack.top();
  ImageStack.pop();
  cmtk::UniformVolume::SmartPtr q = ImageStack.top();
  ImageStack.pop();
  
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
  ImageStack.push( p );

  return NULL;
}

const char*
CallbackSum()
{
  while ( ImageStack.size() > 1 )
    {
    CallbackAdd();
    }
  return NULL;
}

const char*
CallbackAverage()
{
  const size_t nimages = ImageStack.size();
  CallbackSum();
  ImageStack.top()->GetData()->Rescale( 1.0 / nimages );
  return NULL;
}

const char*
CallbackVoteCombination()
{
  if ( ! CheckStackTwoMatchingImages( "Vote" ) )
    return NULL;

  cmtk::UniformVolume::SmartPtr grid = ImageStack.top();
  
  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "Vote" ) )
	return NULL;
    
    dataPtrs.push_back( ImageStack.top()->GetData() );
    ImageStack.pop();
    }

  cmtk::LabelCombinationVoting voting( dataPtrs );
  grid->SetData( voting.GetResult() );
  ImageStack.push( grid );

  return NULL;
}

const char*
CallbackSTAPLE( const long int maxIterations )
{
  if ( ! CheckStackTwoMatchingImages( "STAPLE" ) )
    return NULL;
  
  cmtk::UniformVolume::SmartPtr imageGrid = ImageStack.top();

  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "STAPLE" ) )
	return NULL;
    
    dataPtrs.push_back( ImageStack.top()->GetData() );
    ImageStack.pop();
    }

  cmtk::LabelCombinationSTAPLE staple( dataPtrs, maxIterations, ResultType );
  imageGrid->SetData( staple.GetResult() );
  ImageStack.push( imageGrid );

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

  return NULL;
}

const char*
CallbackMultiClassSTAPLE( const long int maxIterations )
{
  if ( ! CheckStackTwoMatchingImages( "MultiClassSTAPLE" ) )
    return NULL;
  
  cmtk::UniformVolume::SmartPtr imageGrid = ImageStack.top();

  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "MultiClassSTAPLE" ) )
	return NULL;
    
    dataPtrs.push_back( ImageStack.top()->GetData() );
    ImageStack.pop();
    }

  cmtk::LabelCombinationMultiClassSTAPLE mstaple( dataPtrs, maxIterations );
  imageGrid->SetData( mstaple.GetResult() );
  ImageStack.push( imageGrid );

  return NULL;
}

const char*
CallbackStackEntropyLabels()
{
  if ( ! CheckStackTwoMatchingImages( "StackEntropyLabels" ) )
    return NULL;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "StackEntropyLabels" ) )
	return NULL;
    
    volPtrs.push_back( ImageStack.top() );
    ImageStack.pop();
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
  ImageStack.push( volPtrs[0] );
  
  return NULL;
}

const char*
CallbackMaxIndex()
{
  if ( ! CheckStackTwoMatchingImages( "MaxIndex" ) )
    return NULL;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs( ImageStack.size() );

  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "MaxIndex" ) )
	return NULL;
    
    volPtrs[ImageStack.size()-1] = ImageStack.top();
    ImageStack.pop();
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
  ImageStack.push( volPtrs[0] );
  
  return NULL;
}

const char*
CallbackMaxValue()
{
  if ( ! CheckStackTwoMatchingImages( "MaxValue" ) )
    return NULL;
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "MaxValue" ) )
	return NULL;
    
    volPtrs.push_back( ImageStack.top() );
    ImageStack.pop();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr maxArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );

#pragma omp parallel for  
  for ( size_t i = 0; i < numberOfPixels; ++i )
    {
    float maxValue = 0;
    bool maxValueValid = false;
    cmtk::Types::DataItem v;

    for ( size_t curVol = 0; curVol < volPtrs.size(); ++curVol )
      {
      if ( volPtrs[ curVol ]->GetDataAt( v, i ) ) 
        {
	if ( maxValueValid )
	  {
	  maxValue = std::max<float>( maxValue, v );
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
  ImageStack.push( volPtrs[0] );
  
  return NULL;
}

const char*
CallbackCombinePCA()
{
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      if ( ! CheckStackTwoMatchingImages( "CombinePCA" ) )
	return NULL;
    
    volPtrs.push_back( ImageStack.top() );
    ImageStack.pop();
    }
  
  const size_t numberOfPixels = volPtrs[0]->GetNumberOfPixels();
  const size_t numberOfImages = volPtrs.size();
  
  cmtk::Vector<cmtk::Types::DataItem> meanVector( numberOfImages );
  for ( size_t image = 0; image < numberOfImages; ++image )
    {
    cmtk::Types::DataItem mean = 0;
    size_t numberOfPixelsNotPadding = 0;

    for ( size_t pixel = 0; pixel < numberOfPixels; ++pixel ) 
      {      
      cmtk::Types::DataItem value;
      if ( volPtrs[image]->GetDataAt( value, pixel ) )
	{
	mean += value;
	++numberOfPixelsNotPadding;
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
  cmtk::Array<cmtk::Types::DataItem> eigenvalues( numberOfImages );
  
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
  ImageStack.push( volPtrs[0] );

  return NULL;
}

int
main( int argc, char *argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image operations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Perform operations on images using stack-based postfix notation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] [operations]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Be verbose" );

    cl.BeginGroup( "Input/output", "Input/output operations" );
    cl.AddCallback( Key( "in" ), CallbackIn, "Read input image(s) to top of stack" );
    cl.AddCallback( Key( "out" ), CallbackOut, "Write output image from top of stack (but leave it on the stack)" );
    cl.EndGroup();

    cl.BeginGroup( "Internal", "Internal settings" );
    cl.AddSwitch( Key( "float" ), &ResultType, cmtk::TYPE_FLOAT, "Use single precision for computations and results [default]" );
    cl.AddSwitch( Key( "double" ), &ResultType, cmtk::TYPE_DOUBLE, "Use double precision for computations and results" );
    cl.EndGroup();

    cl.BeginGroup( "Stack", "Stack operations" );
    cl.AddCallback( Key( "pop" ), CallbackPop, "Pop (discard) top image from stack" );
    cl.AddCallback( Key( "dup" ), CallbackDup, "Duplicate image on top of the stack" );
    cl.EndGroup();

    cl.BeginGroup( "Single image", "Single-image operators" );
    cl.AddCallback( Key( "abs" ), CallbackAbs, "Apply abs() function to top image" );
    cl.AddCallback( Key( "log" ), CallbackLog, "Apply log() function to top image" );
    cl.AddCallback( Key( "exp" ), CallbackExp, "Apply exp() function to top image" );
    cl.AddCallback( Key( "sqr" ), CallbackSqr, "Apply square operator to top image" );
    cl.AddCallback( Key( "sqrt" ), CallbackSqrt, "Apply square root operator to top image" );
    cl.AddCallback( Key( "trunc" ), CallbackTrunc, "Truncate all values in top image to integer" );
    cl.AddCallback( Key( "one-over" ), CallbackOneOver, "For each pixel, replace its value x with 1.0/x" );
    cl.AddCallback( Key( "scalar-mul" ), CallbackScalarMul, "Multiply top image with a scalar value" );
    cl.AddCallback( Key( "scalar-add" ), CallbackScalarAdd, "Add a scalar to each pixel of the top image" );

    cl.AddCallback( Key( "threshold-min" ), CallbackThreshMin, "Apply lower threshold to image at top of stack" );
    cl.EndGroup();

    cl.BeginGroup( "Two images", "Image pair operators" );
    cl.AddCallback( Key( "add" ), CallbackAdd, "Add top and second image, place result on stack" );
    cl.AddCallback( Key( "mul" ), CallbackMul, "Multiply top and second image, place result on stack" );
    cl.AddCallback( Key( "div" ), CallbackDiv, "Divide top image by second image, place result on stack" );
    cl.AddCallback( Key( "atan2" ), CallbackAtan2, "Compute atan2() function from tup two image pixel pairs, place result on stack" );
    cl.EndGroup();

    cl.BeginGroup( "Multiple images", "Multi-image operators (operate on the entire stack)" );
    cl.AddCallback( Key( "sum" ), CallbackSum, "Sum all images on stack, place result on stack" );
    cl.AddCallback( Key( "average" ), CallbackAverage, "Average all images on stack, place result on stack" );
    cl.AddCallback( Key( "combine-pca" ), CallbackCombinePCA, "Combine images using PCA by projecting onto direction of largest correlation" );
    cl.AddCallback( Key( "max-value" ), CallbackMaxValue, "For each pixel, compute maximum VALUE over all images, place result on stack" );
    cl.AddCallback( Key( "max-index" ), CallbackMaxIndex, "For each pixel, compute INDEX of image with maximum value, place result on stack" );
    cl.EndGroup();

    cl.BeginGroup( "Multiple label images", "Multi-image operators for label images" );
    cl.AddCallback( Key( "vote" ), CallbackVoteCombination, "Merge all images on stack with voting, place result on stack" );
    cl.AddCallback( Key( "staple" ), CallbackSTAPLE, "Combine binary masks on the stack using [arg] iterations of the STAPLE algorithm" );
    cl.AddCallback( Key( "mstaple" ), CallbackMultiClassSTAPLE, "Combine multi-label masks on the stack using [arg] iterations of the multi-class STAPLE algorithm" );
    cl.AddCallback( Key( "stack-entropy-labels" ), CallbackStackEntropyLabels, "Compute stack entropy at each pixel from integer (label) input images" );
    cl.EndGroup();

    cl.Parse();

    }
  catch ( cmtk::CommandLine::Exception e )
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
