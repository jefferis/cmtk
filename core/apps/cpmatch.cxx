/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>
#include <cmtkVolumeIO.h>

#include <cassert>

cmtk::TypedArray*
ComplexToPhase( const cmtk::TypedArray* real, const cmtk::TypedArray* imag )
{
  assert( real->GetDataSize() == imag->GetDataSize() );

  const size_t nSamples = real->GetDataSize();
  cmtk::TypedArray* phase = cmtk::TypedArray::Create( cmtk::TYPE_DOUBLE, nSamples );
#pragma omp parallel for
  for ( size_t n = 0; n < nSamples; ++n )
    {
    cmtk::Types::DataItem r, i;
    if ( real->Get( r, n ) && imag->Get( i, n ) )
      {
      phase->Set( static_cast<cmtk::Types::DataItem>( atan2( i, r ) ), n );
      }
    else
      phase->SetPaddingAt( n );
    }

  return phase;
}

double
EstimatePhaseDifference
( const cmtk::TypedArray* phaseRef, const cmtk::TypedArray* phaseTst )
{
  assert( phaseRef->GetDataSize() == phaseTst->GetDataSize() );

  const double twoPi = 2.0 * M_PI;

  double sum = 0;
  
  const size_t nSamples = phaseRef->GetDataSize();
  for ( size_t n = 0; n < nSamples; ++n )
    {
    cmtk::Types::DataItem pR, pT;
    if ( phaseRef->Get( pR, n ) && phaseTst->Get( pT, n ) )
      {
      const double diff = (pT - pR);
      sum += diff;
//      sum += (diff < -M_PI) ? diff + twoPi : diff;
      }
    }

  return sum / nSamples;
}

int
main( int argc, char* argv[] )
{
  bool verbose = false;

  try
    {
    cmtk::CommandLine cl( argc, argv );

    cl.BeginGroup( "General", "General Parameters" );
    cl.AddSwitch( cmtk::CommandLine::Key( 'v', "verbose" ), &verbose, true, "Verbose operation" );
    cl.EndGroup();

    cl.Parse();
    
    const char* pathRefReal = cl.GetNext();
    const char* pathRefImag = cl.GetNext();
    
    const char* pathChgReal = cl.GetNext();
    const char* pathChgImag = cl.GetNext();
    
    const char* pathOutReal = cl.GetNext();
    const char* pathOutImag = cl.GetNext();

    cmtk::UniformVolume::SmartPtr imageRefReal( cmtk::VolumeIO::Read( pathRefReal, verbose ) );
    cmtk::UniformVolume::SmartPtr imageRefImag( cmtk::VolumeIO::Read( pathRefImag, verbose ) );
    
    cmtk::UniformVolume::SmartPtr imageChgReal( cmtk::VolumeIO::Read( pathChgReal, verbose ) );
    cmtk::UniformVolume::SmartPtr imageChgImag( cmtk::VolumeIO::Read( pathChgImag, verbose ) );

    cmtk::TypedArray::SmartPtr dataRefPhase( ComplexToPhase( imageRefReal->GetData(), imageRefImag->GetData() ) );
    cmtk::TypedArray::SmartPtr dataChgPhase( ComplexToPhase( imageChgReal->GetData(), imageChgImag->GetData() ) );
    
    const double phaseDifference = EstimatePhaseDifference( dataRefPhase, dataChgPhase );
    
    cmtk::StdErr << phaseDifference * 180 / M_PI << "\n";    
    }
  catch ( cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    exit( 1 );
    }
    
  return 0;
}
