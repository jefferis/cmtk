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

#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkGeneralLinearModel.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkVolumeIO.h>

#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <fstream>

void
ParseArgsString( const std::string& in, std::vector<double>& args )
{
  std::string arg = in;
  while ( arg.length() )
    {
    const size_t colon = arg.find( ':' );
  
    if ( colon == std::string::npos )
      {
      args.push_back( atof( arg.c_str() ) );
      arg.clear();
      }
    else
      {
      const std::string next = arg.substr( 0, colon );
      args.push_back( atof( next.c_str() ) );
      arg.erase( 0, colon+1 );
      }
    }
}

void
ParseArgsWarpString( const std::string& in, std::string& fname, std::vector<double>& args )
{
  size_t colon = in.find( ':' );
  if ( colon == std::string::npos )
    {
    cmtk::StdErr << "ERROR: could not parse argument string '" << in << "'; must be ID:x0:[x1...]\n";
    throw cmtk::ExitException( 1 );
    }

  fname = in.substr( 0, colon );
  ParseArgsString( in.substr( colon+1 ), args );
}

int
doMain( const int argc, const char* argv[] )
{
  std::list<std::string> pathList;

  const char* outFileName = NULL;
  const char* controlFileName = NULL;

  int order = 1;
  size_t nParameters = 0;
  size_t nActualParameters = 0;
  std::vector<double> design;
  std::vector<double> meanParam;

  const char* sParameters = NULL;
  std::vector<double> param;

  const char* pathPrintf = "%s";

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Regression" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Linear (and higher-order polynomial) regression of deformation fields and images." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] controlFile" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( 'o', "output" ), &outFileName, "File name for output." );
    cl.AddOption( Key( 'p', "parameters" ), &sParameters, "Model parameters for regression instantiation" );
    cl.AddOption( Key( 's', "substitution" ), &pathPrintf, "Printf format string for ID-to-path substitition" );
    cl.AddOption( Key( 'O', "order" ), &order, "Polynonial order of the regression [default: 1=linear]" );

    cl.Parse( argc, argv );

    controlFileName = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    }

  std::ifstream controlStream( controlFileName );
  if ( controlStream.fail() )
    {
    cmtk::StdErr << "ERROR: could not open control file " << controlFileName << "\n";
    throw cmtk::ExitException( 1 );
    }

  while ( controlStream.good() )
    {
    std::string line;
    controlStream >> line;
    
    cmtk::DebugOutput( 7 ) << line << "\n";

    if ( line.length() )
      {
      std::string fname;
      std::vector<double> args;
      ParseArgsWarpString( line, fname, args );
      
      if ( nActualParameters )
	{
	if ( nActualParameters != args.size() )
	  {
	  cmtk::StdErr << "ERROR: number of model parameters changed from " << nActualParameters << " to " << args.size() << " in argument " << line << "\n";
	  throw cmtk::ExitException( 1 );
	  }
	}
      else
	{
	nActualParameters = args.size();
	nParameters = nActualParameters;
	if ( order > 1 )
	  nParameters += nActualParameters * (nActualParameters+1) / 2;
	meanParam.resize( nParameters );
	}

      size_t ofs = 0;
      for ( size_t i = 0; i < args.size(); ++i )
	{
	design.push_back( args[i] );
	meanParam[ofs++] += args[i];
	}
      if ( order > 1 )
	{
	for ( size_t j = 0; j < nActualParameters; ++j )
	  {
	  for ( size_t i = j; i < nActualParameters; ++i )
	    {
	    design.push_back( args[i] * args[j] );
	    meanParam[ofs++] += args[i] * args[j];
	    }
	  }
	}
      design.push_back( 1.0 );
      
      char path[PATH_MAX];
      snprintf( path, PATH_MAX, pathPrintf, fname.c_str() );
      pathList.push_back( std::string( path ) );
      }
    }

  const size_t nData = pathList.size();
  cmtk::GeneralLinearModel glm( 1+nParameters, nData, &(design[0]) );

  for ( size_t i = 0; i < meanParam.size(); ++i )
    {
    meanParam[i] /= nData;
    }
  
  if ( sParameters )
    {
    ParseArgsString( std::string( sParameters ), param );

    if ( order > 1 )
      {
      size_t ofs = nActualParameters;
      param.resize( nParameters );
      for ( size_t j = 0; j < nActualParameters; ++j )
	{
	for ( size_t i = j; i < nActualParameters; ++i )
	  {
	  param[ofs++] += param[i] * param[j];
	  }
	}
      }
    }
  else
    {
    param = meanParam;
    }

  if ( cmtk::FileFormat::Identify( pathList.begin()->c_str() ) == cmtk::FILEFORMAT_TYPEDSTREAM )
    // xform mode
    {
    cmtk::DebugOutput( 1 ) << "INFO: operating in XFORM mode\n";

    std::vector<cmtk::SplineWarpXform::SmartPtr> vWarpXform;
    for ( std::list<std::string>::const_iterator it = pathList.begin(); it != pathList.end(); ++it )
      {
      cmtk::SplineWarpXform::SmartPtr xform = cmtk::SplineWarpXform::SmartPtr::DynamicCastFrom( cmtk::Xform::SmartPtr( cmtk::XformIO::Read( it->c_str() ) ) );
      if ( ! xform )
	{
	cmtk::StdErr << "ERROR: transformation '" << *it << "' is either invalid or not a spline warp xform\n";
	throw cmtk::ExitException( 1 );
	}
      vWarpXform.push_back( xform );
      }
    
    cmtk::SplineWarpXform::SmartPtr referenceWarp = vWarpXform[0];
    for ( size_t i = 0; i < vWarpXform.size(); ++i )
      {
      vWarpXform[i]->ReplaceInitialAffine( NULL );
      }
    
    const size_t nWarpParameters = vWarpXform[0]->m_NumberOfParameters;
    std::vector<double> mean( nWarpParameters );
    std::fill( mean.begin(), mean.end(), 0 );
    
    std::vector<cmtk::TypedArray::SmartPtr> vData;
    for ( size_t i = 0; i <vWarpXform.size(); ++i )
      {
      cmtk::Types::Coordinate* deformation = vWarpXform[i]->GetPureDeformation();
      cmtk::TypedArray::SmartPtr data( cmtk::TypedArray::Create( cmtk::TYPE_COORDINATE, deformation, nWarpParameters ) );
      for ( size_t n = 0; n < nWarpParameters; ++n )
	{
	mean[n] += deformation[n];
	}
      vData.push_back( data );
      
      vWarpXform[i] = cmtk::SplineWarpXform::SmartPtr::Null(); // no longer needed
      }
    
    for ( size_t n = 0; n < mean.size(); ++n )
      mean[n] /= nData;
    
    glm.FitModel( vData, false /*normalizeParameters*/ );
    
    for ( size_t p = 0; p < nParameters; ++p )
      {
      cmtk::TypedArray::SmartPtr mode = glm.GetModel( p );
      
      for ( size_t n = 0; n < nWarpParameters; ++n )
	{
	cmtk::Types::DataItem value;
	if ( mode->Get( value, n ) )
	  {
	  mean[n] += value * (param[p] - meanParam[p]);
	  }
	}
      }
    
    for ( size_t n = 0; n < nWarpParameters; ++n )
      {
      referenceWarp->SetParameter( n, mean[n] );
      }
    
    if ( outFileName )
      cmtk::XformIO::Write( referenceWarp, outFileName );
    }
  else
    // image mode
    {
    cmtk::DebugOutput( 1 ) << "INFO: operating in IMAGE mode\n";

    std::vector<cmtk::UniformVolume::SmartPtr> vVolume;
    for ( std::list<std::string>::const_iterator it = pathList.begin(); it != pathList.end(); ++it )
      {
      cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( it->c_str() ) );
      if ( ! volume )
	{
	cmtk::StdErr << "ERROR: image '" << *it << "' could not be read\n";
	throw cmtk::ExitException( 1 );
	}
      vVolume.push_back( volume );
      }
    
    cmtk::UniformVolume::SmartPtr referenceVolume = vVolume[0];
    
    const size_t nVolumeParameters = referenceVolume->GetNumberOfPixels();
    std::vector<double> mean( nVolumeParameters );
    std::fill( mean.begin(), mean.end(), 0.0 );
    
    std::vector<cmtk::TypedArray::SmartPtr> vData;
    for ( size_t i = 0; i <vVolume.size(); ++i )
      {
      vData.push_back( vVolume[i]->GetData() );
      }
    
    glm.FitModel( vData, false /*normalizeParameters*/ );
    
    for ( size_t p = 0; p < nParameters; ++p )
      {
      cmtk::TypedArray::SmartPtr mode = glm.GetModel( p );
      
      for ( size_t n = 0; n < nVolumeParameters; ++n )
	{
	cmtk::Types::DataItem value;
	if ( mode->Get( value, n ) )
	  {
	  mean[n] += value * param[p];
	  }
	}
      }
    
    cmtk::TypedArray::SmartPtr constant = glm.GetModel( nParameters );
    
    for ( size_t n = 0; n < nVolumeParameters; ++n )
      {
      cmtk::Types::DataItem value;
      if ( constant->Get( value, n ) )
	{
	mean[n] += value;
	}
      }
    
    referenceVolume->CreateDataArray( cmtk::TYPE_FLOAT );
    for ( size_t n = 0; n < nVolumeParameters; ++n )
      {
      referenceVolume->SetDataAt( mean[n], n );
      }
    
    if ( outFileName )
      cmtk::VolumeIO::Write( *referenceVolume, outFileName );
    }

  return 0;
}

#include "cmtkSafeMain"
