/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkThreads.h>
#include <System/cmtkDebugOutput.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkMathFunctionWrappers.h>
#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkVolumeIO.h>

#include <Segmentation/cmtkLabelCombinationVoting.h>
#include <Segmentation/cmtkLabelCombinationSTAPLE.h>
#include <Segmentation/cmtkLabelCombinationMultiClassSTAPLE.h>

#include <math.h>
#include <queue>
#include <vector>

#ifdef CMTK_USE_GCD
#  include <dispatch/dispatch.h>
#endif

/// The output data type of image operations.
cmtk::ScalarDataType ResultType = cmtk::TYPE_FLOAT;

/// Value used for padding input data, if PaddingFlag is true.
cmtk::Types::DataItem PaddingValue;

/// If this flag is set, PaddingValue defines the padding value for images read from files.
bool PaddingFlag = false;

void
CallbackSetPaddingValue( const double paddingValue )
{
  cmtk::DebugOutput( 2 ) << "CallbackSetPaddingValue\n";
  PaddingValue = paddingValue;
  PaddingFlag = true;
}

void
CallbackUnsetPaddingValue()
{
  cmtk::DebugOutput( 2 ) << "CallbackUnsetPaddingValue\n";
  PaddingFlag = false;
}

/// Type for "stack" of images.
typedef std::deque<cmtk::UniformVolume::SmartPtr> ImageStackType;

/// Operating stack of images.
ImageStackType ImageStack;

/// If this flag is set, the next single-image operation will be applied to all images on the stack.
bool ApplyNextToAll = false;

void
CheckStackOneImage( const char* function )
{
  if ( ImageStack.empty() )
    {
    cmtk::StdErr << "ERROR: stack is empty in function '" << function << "'\n";
    throw cmtk::ExitException( 1 );
    }
}

void
CheckStackTwoMatchingImages( const char* function )
{
  if ( ImageStack.size() < 2 )
    {
    cmtk::StdErr << "ERROR: at least two images are required to perform operation '" << function << "'\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr top = ImageStack.front();
  ImageStack.pop_front();
  const bool pixelCountsMatch = (top->GetNumberOfPixels() == ImageStack.front()->GetNumberOfPixels() );
  ImageStack.push_front( top );

  if ( !pixelCountsMatch )
    {
    cmtk::StdErr << "ERROR: two top images do not have equal pixel counts as required for operation '" << function << "'\n";
    throw cmtk::ExitException( 1 );
    }
}

void
CallbackIn( const char** argv, int& argsUsed )
{
  cmtk::DebugOutput( 2 ) << "CallbackIn\n";
  argsUsed = 0;
  while ( argv[argsUsed][0] != '-' )
    {
    cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( argv[argsUsed] ) );
    if ( !volume || !volume->GetData() )
      {
      cmtk::StdErr << "ERROR: could not read input image " << argv[argsUsed] << "\n";
      throw cmtk::ExitException( 1 );
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
  cmtk::DebugOutput( 2 ) << "CallbackOut\n";
  CheckStackOneImage( "Out" );
  cmtk::VolumeIO::Write( *(ImageStack.front()), argv );
}

void
CallbackPop()
{
  cmtk::DebugOutput( 2 ) << "CallbackPop\n";
  CheckStackOneImage( "Pop" );
  ImageStack.pop_front();
}

void
CallbackDup()
{
  cmtk::DebugOutput( 2 ) << "CallbackDup\n";
  CheckStackOneImage( "Dup" );
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
  cmtk::DebugOutput( 2 ) << "CallbackFill\n";
  CheckStackOneImage( "Fill" );

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
    
void
CallbackAbs()
{
  cmtk::DebugOutput( 2 ) << "CallbackAbs\n";
  CheckStackOneImage( "Abs" );

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
    
void
CallbackTrunc()
{
  CheckStackOneImage( "Trunc" );

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
    
void
CallbackLog()
{
  CheckStackOneImage( "Log" );

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
    
void
CallbackLogit()
{
  CheckStackOneImage( "Logit" );

  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
    (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Logit );
    
    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}
    
void
CallbackLogistic()
{
  CheckStackOneImage( "Logistic" );
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
    (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Logistic );
    
    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}
    
void
CallbackExp()
{
  CheckStackOneImage( "Exp" );

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

void
CallbackSqr()
{
  CheckStackOneImage( "Sqr" );

  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    (*it)->SetData( cmtk::TypedArray::SmartPtr( (*it)->GetData()->Convert( ResultType ) ) );
    (*it)->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Square );
    
    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}

void
CallbackSqrt()
{
  CheckStackOneImage( "Sqrt" );

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
        
void
CallbackThreshBelow( const char* argv )
{
  CheckStackOneImage( "ThreshBelow" );

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
    
void
CallbackThreshAbove( const char* argv )
{
  CheckStackOneImage( "ThreshAbove" );

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
    
void
CallbackScalarMul( const double c )
{
  CheckStackOneImage( "ScalarMul" );
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    cmtk::UniformVolume& p = **it;
    const size_t numberOfPixels = p.GetNumberOfPixels();
    
    cmtk::TypedArray::SmartPtr mul( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
    cmtk::TypedArray* mulRef = mul.GetPtr(); // need this to work around GCD / clang limitations
    
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numberOfPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem pv;
      if ( p.GetDataAt( pv, i ) )
	{
	mulRef->Set( c * pv, i );
	}
      else
	{
	mulRef->SetPaddingAt( i );
	}
      }
#ifdef CMTK_USE_GCD
		    } );
#endif
    
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
  CheckStackOneImage( "ScalarAdd" );
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    cmtk::UniformVolume& p = **it;
    const size_t numberOfPixels = p.GetNumberOfPixels();
    
    cmtk::TypedArray::SmartPtr add( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
    cmtk::TypedArray* addRef = add.GetPtr(); // need this to work around GCD / clang limitations

#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numberOfPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem pv;
      if ( p.GetDataAt( pv, i ) )
	{
	addRef->Set( c + pv, i );
	}
      else
	{
	addRef->SetPaddingAt( i );
	}
      }
#ifdef CMTK_USE_GCD
		    });
#endif    

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
  CheckStackOneImage( "ScalarXor" );
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    cmtk::UniformVolume& p = **it;
    const size_t numberOfPixels = p.GetNumberOfPixels();
    
    cmtk::TypedArray::SmartPtr out( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
    cmtk::TypedArray* outRef = out.GetPtr();  // need this to work around GCD / clang limitations
    
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numberOfPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem pv;
      if ( p.GetDataAt( pv, i ) )
	{
	const long int iv = static_cast<long int>( pv );
	outRef->Set( iv ^ c, i );
	}
      else
	{
	outRef->SetPaddingAt( i );
	}
      }
#ifdef CMTK_USE_GCD
		    });
#endif

    p.SetData( out );

    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}

void
CallbackScalarAnd( const long int c )
{
  CheckStackOneImage( "ScalarAnd" );
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    cmtk::UniformVolume& p = **it;
    const size_t numberOfPixels = p.GetNumberOfPixels();
    
    cmtk::TypedArray::SmartPtr out( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
    cmtk::TypedArray* outRef = out.GetPtr(); // need this to work around GCD / clang limitations
    
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numberOfPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem pv;
      if ( p.GetDataAt( pv, i ) )
	{
	const long int iv = static_cast<long int>( pv );
	outRef->Set( iv & c, i );
	}
      else
	{
	outRef->SetPaddingAt( i );
	}
      }
#ifdef CMTK_USE_GCD
		    });
#endif    

    p.SetData( out );

    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}

void
CallbackOneOver()
{
  CheckStackOneImage( "OneOver" );
  
  ImageStackType::iterator it = ImageStack.begin();
  while ( it != ImageStack.end() )
    {
    cmtk::UniformVolume& p = **it;
    const size_t numberOfPixels = p.GetNumberOfPixels();
    
    cmtk::TypedArray::SmartPtr inv( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
    cmtk::TypedArray* invRef = inv.GetPtr(); // need this to work around GCD / clang limitations
    
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numberOfPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem pv;
      if ( p.GetDataAt( pv, i ) )
	{
	invRef->Set( 1.0 / pv, i );
	}
      else
	{
	invRef->SetPaddingAt( i );
	}
      }
#ifdef CMTK_USE_GCD
		    });
#endif
    
    p.SetData( inv );

    if ( ApplyNextToAll )
      ++it;
    else
      it = ImageStack.end();
    }
  ApplyNextToAll = false;
}

void
CallbackAdd()
{
  CheckStackTwoMatchingImages( "Add" );
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();
  cmtk::UniformVolume::SmartPtr q = ImageStack.front();
  ImageStack.pop_front();
  
  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr add( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* addRef = add.GetPtr(); // need this to work around GCD / clang limitations
  
#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
		      for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    cmtk::Types::DataItem pv, qv;
    if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) )
      {
      addRef->Set( pv+qv, i );
      }
    else
      {
      addRef->SetPaddingAt( i );
      }
    }
#ifdef CMTK_USE_GCD
		    });
#endif
  
  p->SetData( add );
  ImageStack.push_front( p );
}

void
CallbackMul()
{
  CheckStackTwoMatchingImages( "Mul" );
  
  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();

  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr mul( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* mulRef = mul.GetPtr(); // need this to work around GCD / clang limitations
  
  while ( ! ImageStack.empty() )
    {
    cmtk::UniformVolume::SmartPtr q = ImageStack.front();
    ImageStack.pop_front();
    
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numberOfPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem pv, qv;
      if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) )
	{
	mulRef->Set( pv*qv, i );
	}
      else
	{
	mulRef->SetPaddingAt( i );
	}
      }
#ifdef CMTK_USE_GCD
		    });
#endif    

    if ( !ApplyNextToAll )
      break;
    }

  ApplyNextToAll = false;
    
  p->SetData( mul );
  ImageStack.push_front( p );
}

void
CallbackDiv()
{
  CheckStackTwoMatchingImages( "Div" );

  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();
  cmtk::UniformVolume::SmartPtr q = ImageStack.front();
  ImageStack.pop_front();
  
  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr div( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* divRef = div.GetPtr(); // need this to work around GCD / clang limitations
  
#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
		      for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    cmtk::Types::DataItem pv, qv;
    if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) && (qv != 0) )
      {
      divRef->Set( pv/qv, i );
      }
    else
      {
      divRef->SetPaddingAt( i );
      }
    }
#ifdef CMTK_USE_GCD
		  });
#endif  

  p->SetData( div );
  ImageStack.push_front( p );
}

void
CallbackAtan2()
{
  CheckStackTwoMatchingImages( "Atan2" );

  cmtk::UniformVolume::SmartPtr p = ImageStack.front();
  ImageStack.pop_front();
  cmtk::UniformVolume::SmartPtr q = ImageStack.front();
  ImageStack.pop_front();
  
  const size_t numberOfPixels = p->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* resultRef = result.GetPtr(); // need this to work around GCD / clang limitations
  
#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
		      for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    cmtk::Types::DataItem pv, qv;
    if ( p->GetDataAt( pv, i ) && q->GetDataAt( qv, i ) && (qv != 0) )
      {
      resultRef->Set( atan2( qv, pv ), i );
      }
    else
      {
      resultRef->SetPaddingAt( i );
      }
    }
#ifdef CMTK_USE_GCD
		  });
#endif  

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

  cmtk::TypedArray::SmartPtr result( ImageStack.front()->GetData()->Convert( ResultType ) );
  result->ApplyFunctionObject( cmtk::TypedArrayFunctionHistogramMatching( *(ImageStack.front()->GetData()), *(ref->GetData()) ) );
  ImageStack.front()->SetData( result );
}

void
CallbackMatchMeanSDev()
{
  if ( ImageStack.size() < 2 )
    {
    cmtk::StdErr << "ERROR: need at least two images on stack for histogram intensity matching\n";
    return;
    }
  
  cmtk::UniformVolume::SmartPtr ref = ImageStack.front();
  ImageStack.pop_front();

  cmtk::Types::DataItem rMean, rVar;
  ref->GetData()->GetStatistics( rMean, rVar );

  cmtk::Types::DataItem mMean, mVar;
  ImageStack.front()->GetData()->GetStatistics( mMean, mVar );

  const cmtk::Types::DataItem scale = sqrt( rVar / mVar );
  const cmtk::Types::DataItem offset = rMean - scale * mMean;

  cmtk::TypedArray::SmartPtr result( ImageStack.front()->GetData()->Convert( ResultType ) );
  result->Rescale( scale, offset );
  ImageStack.front()->SetData( result );
}

void
CallbackMatchMeanSDevThree()
{
  if ( ImageStack.size() < 3 )
    {
    cmtk::StdErr << "ERROR: need at least three images on stack for histogram intensity matching using separate reference images\n";
    return;
    }
  
  cmtk::UniformVolume::SmartPtr ref = ImageStack.front();
  ImageStack.pop_front();

  cmtk::UniformVolume::SmartPtr mod = ImageStack.front();
  ImageStack.pop_front();

  cmtk::Types::DataItem rMean, rVar;
  ref->GetData()->GetStatistics( rMean, rVar );

  cmtk::Types::DataItem mMean, mVar;
  mod->GetData()->GetStatistics( mMean, mVar );
  
  const cmtk::Types::DataItem scale = sqrt( rVar / mVar );
  const cmtk::Types::DataItem offset = rMean - scale * mMean;
  
  cmtk::TypedArray::SmartPtr result( ImageStack.front()->GetData()->Convert( ResultType ) );
  result->Rescale( scale, offset );
  ImageStack.front()->SetData( result );
}

void
CallbackMaskAverage()
{
  if ( ImageStack.size() < 2 )
    {
    cmtk::StdErr << "ERROR: need at last two images on stack for mask averaging\n";
    return;
    }
  
  cmtk::UniformVolume::SmartPtr msk = ImageStack.front();
  ImageStack.pop_front();
  
  const cmtk::TypedArray& mskData = *(msk->GetData());
  const cmtk::Types::DataItemRange labelRange = mskData.GetRange();
  
  const size_t nLabels = static_cast<size_t>( labelRange.Width()+1 );
  std::vector<double> means( nLabels, 0.0 );
  std::vector<size_t> count( nLabels, 0 );

  cmtk::TypedArray::SmartPtr result( ImageStack.front()->GetData()->Convert( ResultType ) );
  ImageStack.front()->SetData( result );

  cmtk::TypedArray& avgData = *(ImageStack.front()->GetData());

  // first pass: compute mean for each labelled ROI
  const size_t n = mskData.GetDataSize();
  for ( size_t idx = 0; idx < n; ++idx )
    {
    cmtk::Types::DataItem l;
    mskData.Get( l, idx );

    cmtk::Types::DataItem v;
    if ( avgData.Get( v, idx ) )
      {
      const size_t ll = static_cast<size_t>( l-labelRange.m_LowerBound );
      means[ll] += v;
      ++count[ll];
      }
    }

  // second pass: replace values with computed ROI means
#pragma omp parallel for
  for ( int idx = 0; idx < static_cast<int>( n ); ++idx )
    {
    cmtk::Types::DataItem l;
    mskData.Get( l, idx );

    const size_t ll = static_cast<size_t>( l-labelRange.m_LowerBound );
    avgData.Set( means[ll] / count[ll], idx );
    }
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
CallbackProduct()
{
  while ( ImageStack.size() > 1 )
    {
    CallbackMul();
    }
}

void
CallbackAverage()
{
  const size_t nimages = ImageStack.size();
  if ( ! nimages )
    {
    cmtk::StdErr << "ERROR: need at least two images on stack for averaging\n";
    return;
    }

  CallbackSum();
  ImageStack.front()->GetData()->Rescale( 1.0 / nimages );
}

void
CallbackVariance()
{
  CheckStackTwoMatchingImages( "Variance" );
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "Variance" );
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const std::vector<cmtk::UniformVolume::SmartPtr>& vPtrs = volPtrs; // using const reference, we prevent GCD from segfaulting

  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr varianceArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* varianceArrayRef = varianceArray.GetPtr(); // need this to work around GCD / clang limitations

#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for  
  for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    double sum = 0;
    double sum2 = 0;
    size_t nvalues = 0;
    for ( size_t curVol = 0; curVol < vPtrs.size(); ++curVol )
      {
      cmtk::Types::DataItem v;
      if ( vPtrs[ curVol ]->GetDataAt( v, i ) ) 
        {
	++nvalues;
	sum += v;
	sum2 += v*v;
        }
      }
    
    if ( nvalues )
      {
      const double mean = sum / nvalues;
      varianceArrayRef->Set( (nvalues * mean * mean - 2 * mean * sum + sum2)/nvalues, i );      
      }
    else
      varianceArrayRef->SetPaddingAt( i );       
    }
#ifdef CMTK_USE_GCD
		  });
#endif

  volPtrs[0]->SetData( varianceArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackVoteCombination()
{
  CheckStackTwoMatchingImages( "Vote" );

  cmtk::UniformVolume::SmartPtr grid = ImageStack.front();
  
  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "Vote" );
    
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
  CheckStackTwoMatchingImages( "STAPLE" );
  
  cmtk::UniformVolume::SmartPtr imageGrid = ImageStack.front();

  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "STAPLE" );
    
    dataPtrs.push_back( ImageStack.front()->GetData() );
    ImageStack.pop_front();
    }

  cmtk::LabelCombinationSTAPLE staple( dataPtrs, maxIterations, ResultType );
  imageGrid->SetData( staple.GetResult() );
  ImageStack.push_front( imageGrid );

  if ( cmtk::DebugOutput::GetGlobalLevel() > 0 )
    {
    cmtk::DebugOutput( 1 ) << "p  ";
    for ( size_t i = 0; i < dataPtrs.size(); ++i )
      {
      cmtk::DebugOutput( 1 ).GetStream().printf( "%.3f ", staple.GetPValue( i ) );
      }
    cmtk::DebugOutput( 1 ).GetStream().printf( "\nq  " );
    for ( size_t i = 0; i < dataPtrs.size(); ++i )
      {
      cmtk::DebugOutput( 1 ).GetStream().printf( "%.3f ", staple.GetQValue( i ) );
      }
    cmtk::DebugOutput( 1 ).GetStream().printf( "\n" );
    }
}

void
CallbackMultiClassSTAPLE( const long int maxIterations )
{
  CheckStackTwoMatchingImages( "MultiClassSTAPLE" );
  
  cmtk::UniformVolume::SmartPtr imageGrid = ImageStack.front();

  std::vector<cmtk::TypedArray::SmartPtr> dataPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "MultiClassSTAPLE" );
    
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
  CheckStackTwoMatchingImages( "StackEntropyLabels" );
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "StackEntropyLabels" );
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const std::vector<cmtk::UniformVolume::SmartPtr>& vPtrs = volPtrs; // using const reference, we prevent GCD from segfaulting

  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr entropyArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* entropyArrayRef = entropyArray.GetPtr(); // need this to work around GCD / clang limitations

#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for  
  for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    std::map<int,unsigned int> labelCount;

    size_t totalCount = 0;
    for ( size_t curVol = 0; curVol < vPtrs.size(); ++curVol )
      {
      cmtk::Types::DataItem v;
      if ( vPtrs[ curVol ]->GetDataAt( v, i ) ) 
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
      entropyArrayRef->Set( -entropy / log( 2.0 ), i );
      }
    else
      entropyArrayRef->SetPaddingAt( i );       
    }
#ifdef CMTK_USE_GCD
		  });
#endif

  volPtrs[0]->SetData( entropyArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackMaxIndex()
{
  CheckStackTwoMatchingImages( "MaxIndex" );
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs( ImageStack.size() );

  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "MaxIndex" );
    
    volPtrs[ImageStack.size()-1] = ImageStack.front();
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr maxArray( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, numberOfPixels ) );
  cmtk::TypedArray* maxArrayRef = maxArray.GetPtr(); // need this to work around GCD / clang limitations

  const std::vector<cmtk::UniformVolume::SmartPtr>& vPtrs = volPtrs; // using const reference, we prevent GCD from segfaulting

#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for  
  for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    float maxValue = 0;
    short maxIndex = -2;
    cmtk::Types::DataItem v;
    
    for ( size_t curVol = 0; curVol < vPtrs.size(); ++curVol )
      {
      if ( vPtrs[ curVol ]->GetDataAt( v, i ) ) 
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
    
    maxArrayRef->Set( maxIndex, i ); 
    }
#ifdef CMTK_USE_GCD
		  });
#endif
  
  volPtrs[0]->SetData( maxArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackMaxValue()
{
  CheckStackTwoMatchingImages( "MaxValue" );
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "MaxValue" );
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr maxArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* maxArrayRef = maxArray.GetPtr(); // need this to work around GCD / clang limitations

  const std::vector<cmtk::UniformVolume::SmartPtr>& vPtrs = volPtrs; // using const reference, we prevent GCD from segfaulting

#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for  
  for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    cmtk::Types::DataItem maxValue = 0;
    bool maxValueValid = false;
    cmtk::Types::DataItem v;

    for ( size_t curVol = 0; curVol < vPtrs.size(); ++curVol )
      {
      if ( vPtrs[ curVol ]->GetDataAt( v, i ) ) 
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
    
    
    maxArrayRef->Set( maxValue, i ); 
    }
#ifdef CMTK_USE_GCD
		  });
#endif
  
  volPtrs[0]->SetData( maxArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackMinValue()
{
  CheckStackTwoMatchingImages( "MinValue" );
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "MinValue" );
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr minArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* minArrayRef = minArray.GetPtr(); // need this to work around GCD / clang limitations

  const std::vector<cmtk::UniformVolume::SmartPtr>& vPtrs = volPtrs; // using const reference, we prevent GCD from segfaulting

#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for  
  for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    cmtk::Types::DataItem minValue = 0;
    bool minValueValid = false;
    cmtk::Types::DataItem v;

    for ( size_t curVol = 0; curVol < vPtrs.size(); ++curVol )
      {
      if ( vPtrs[ curVol ]->GetDataAt( v, i ) ) 
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
       
    minArrayRef->Set( minValue, i );
    }
#ifdef CMTK_USE_GCD
		  });
#endif
  
  volPtrs[0]->SetData( minArray );
  ImageStack.push_front( volPtrs[0] );
}

void
CallbackContractLabels()
{
  CheckStackTwoMatchingImages( "ContractLabels" );
  
  std::vector<cmtk::UniformVolume::SmartPtr> volPtrs;
  while ( ImageStack.size() > 0 ) 
    {
    if ( ImageStack.size() > 1 )
      CheckStackTwoMatchingImages( "ContractLabels" );
    
    volPtrs.push_back( ImageStack.front() );
    ImageStack.pop_front();
    }
  
  const size_t numberOfPixels = volPtrs[ 0 ]->GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr outArray( cmtk::TypedArray::Create( ResultType, numberOfPixels ) );
  cmtk::TypedArray* outArrayRef = outArray.GetPtr();

  const std::vector<cmtk::UniformVolume::SmartPtr>& vPtrs = volPtrs; // using const reference, we prevent GCD from segfaulting

#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numberOfPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for  
  for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
#endif
    {
    cmtk::Types::DataItem v = 0;

    for ( size_t curVol = 0; curVol < vPtrs.size(); ++curVol )
      {
      if ( vPtrs[ curVol ]->GetDataAt( v, i ) && v ) 
	break;
      }
        
    outArrayRef->Set( v, i );
    }
#ifdef CMTK_USE_GCD
		  });
#endif
  
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
      CheckStackTwoMatchingImages( "CombinePCA" );
    
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
  
  if ( cmtk::DebugOutput::GetGlobalLevel() > 0 )
    {
    cmtk::DebugOutput( 1 ) << "Principal eigenvector (normalized):\n";
    for ( size_t l = 0; l < numberOfImages; ++l )
      {
      cmtk::DebugOutput( 1 ).GetStream().printf( "%f\t", ev[l] );
      }
    cmtk::DebugOutput( 1 ) << "\n";
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
doMain( const int argc, const char *argv[] )
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image operations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Perform operations on images using stack-based postfix notation.\n\n"
		       "Images can be read from files and pushed onto the stack. Images on the stack can be processed and combined via different operators. "
		       "Results of all operations are put back onto the stack, where they can be further processed or written back to image files." );

    typedef cmtk::CommandLine::Key Key;
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

    cl.AddCallback( Key( "thresh-below" ), CallbackThreshBelow, "Set values below given threshold to threshold." );
    cl.AddCallback( Key( "thresh-above" ), CallbackThreshAbove, "Set values above given threshold to threshold." );
    cl.EndGroup();

    cl.BeginGroup( "Two images", "Image pair operators" );
    cl.AddCallback( Key( "add" ), CallbackAdd, "Add top and second image, place result on stack" );
    cl.AddCallback( Key( "mul" ), CallbackMul, "Multiply top and second image, place result on stack" );
    cl.AddCallback( Key( "div" ), CallbackDiv, "Divide top image by second image, place result on stack" );
    cl.AddCallback( Key( "atan2" ), CallbackAtan2, "Compute atan2() function from tup two image pixel pairs, place result on stack" );    
    cl.AddCallback( Key( "match-histograms" ), CallbackMatchHistograms, "Scale intensities in one image to match intensities of another. The last image pushed onto the stack provides the reference intensity distribution, the preceding image will be modified. Both input images are removed from the stack and the modified image is pushed onto the stack." );
    cl.AddCallback( Key( "match-mean-sdev" ), CallbackMatchMeanSDev, "Scale intensities of one image to match mean and standard deviation of another. The last image pushed onto the stack provides the reference intensity distribution, the preceding image will be modified. Both input images are removed from the stack and the modified image is pushed onto the stack." );
    cl.AddCallback( Key( "match-mean-sdev3" ), CallbackMatchMeanSDevThree, "Scale intensities of an image by a factor and offset computed from two other images to match their mean and standard deviations. The last image pushed onto the stack provides the reference intensity distribution, the preceding image provides the intensity distribution to match to the reference image's, and the third image on the stack will be modified. All three input images are removed from the stack and the modified image is pushed onto the stack." );
    cl.AddCallback( Key( "mask-average" ), CallbackMaskAverage, "Mask averaging: the top image is taken as a multi-label mask. The pixels in the second image are averaged by mask labels, and then replaced with the average value for each mask label." );
    cl.EndGroup();

    cl.BeginGroup( "Contract multiple images", "Operators that contract the entire stack into a single image" );
    cl.AddCallback( Key( "sum" ), CallbackSum, "Sum all images on stack, place result on stack" );
    cl.AddCallback( Key( "product" ), CallbackProduct, "Compute product of all images on stack, place result on stack" );
    cl.AddCallback( Key( "average" ), CallbackAverage, "Average all images on stack, place result on stack" );
    cl.AddCallback( Key( "variance" ), CallbackVariance, "For each pixel, compute variance over all images on stack, place result on stack" );
    cl.AddCallback( Key( "combine-pca" ), CallbackCombinePCA, "Combine images using PCA by projecting onto direction of largest correlation" );
    cl.AddCallback( Key( "max-value" ), CallbackMaxValue, "For each pixel, compute maximum VALUE over all images, place result on stack" );
    cl.AddCallback( Key( "min-value" ), CallbackMinValue, "For each pixel, compute minimum VALUE over all images, place result on stack" );
    cl.AddCallback( Key( "max-index" ), CallbackMaxIndex, "For each pixel, compute INDEX of image with maximum value, place result on stack" );
    cl.EndGroup();

    cl.BeginGroup( "Contract multiple label images", "Operators that contract a stack of label images into a single label image" );
    cl.AddCallback( Key( "vote" ), CallbackVoteCombination, "Merge all images on stack with voting, place result on stack" );
    cl.AddCallback( Key( "staple" ), CallbackSTAPLE, "Combine binary maps on the stack using [arg] iterations of the STAPLE algorithm. "
		    "The result of this operation is the spatial map of 'weights' W, which are the probabilities of image foreground at each pixel. In 'verbose' "
		    "mode, estimated expert parameters p (sensitivity) and q (specificity) are also written to standard output." );
    cl.AddCallback( Key( "contract-labels" ), CallbackContractLabels, "Contract multiple label maps into one by selecting the first (over all images on the stack) non-zero label at each pixel" );
    cl.AddCallback( Key( "mstaple" ), CallbackMultiClassSTAPLE, "Combine multi-label maps on the stack using [arg] iterations of the multi-class STAPLE algorithm."
		    "The result of this operation is the combined maximum-likeliood multi-label map." );
    cl.AddCallback( Key( "stack-entropy-labels" ), CallbackStackEntropyLabels, "Compute stack entropy at each pixel from integer (label) input images" );
    cl.EndGroup();

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  if ( ImageStack.size() > 1 )
    {
    cmtk::StdErr << "WARNING: more than one image left on stack. Are you sure you aren't missing something?\n";
    }
  
  return 0;
}

#include "cmtkSafeMain"
