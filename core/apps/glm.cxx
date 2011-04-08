/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
#include <System/cmtkProgressConsole.h>
#include <System/cmtkMemory.h>

#include <Base/cmtkTypedArray.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkGeneralLinearModel.h>

#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkMathFunctionWrappers.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <list>
#include <set>
#include <string>
#include <stdio.h>
#include <limits.h>

bool DumpStatistics = true;

std::set<int> IgnoreSet;
std::set<std::string> SelectSet;
std::vector<std::string> FieldNames;

void
CallbackIgnore( const char* argv )
{
  IgnoreSet.insert( atoi( argv ) );
}

void
CallbackSelect( const char* argv )
{
  SelectSet.insert( argv );
}

bool ExcludeConstant = false;
bool NormalizeParameters = false;
bool ExponentialModel = false;

std::list<const char*> CtlFileName;
std::list<const char*> ImgFilePatt;

const char* OutputFilePatt = "model_%s_%02d_%s.nii";

bool CropImages = false;
cmtk::DataGrid::RegionType CropImagesRegion;

void
CallbackCropImages( const char* arg )
{
  int cropFrom[3], cropTo[3];
  CropImages = (6 == sscanf( arg, "%d,%d,%d,%d,%d,%d", cropFrom, cropFrom+1, cropFrom+2, cropTo,cropTo+1,cropTo+2 ) );

  if ( CropImages )
    {
    CropImagesRegion = cmtk::DataGrid::RegionType( cmtk::DataGrid::IndexType( cropFrom ), cmtk::DataGrid::IndexType( cropTo ) );
    }
}

void
Import
( const char* ctlFileName, const char* imgFilePatt, 
  cmtk::UniformVolume::SmartPtr& refVolume, 
  std::vector<cmtk::TypedArray::SmartPtr>& imagesData, 
  size_t& nParameters, cmtk::Types::DataItem*& parameters, size_t& nParametersTotal )
{
  nParameters = 0;
  std::vector<cmtk::Types::DataItem> vParam;

  std::ifstream ctlFile( ctlFileName );
  if ( ! ctlFile.is_open() )
    {
    std::cerr << "Error opening control/parameter file " << ctlFileName << std::endl;
    throw cmtk::ExitException( 1 );
    }

  FieldNames.clear();
  std::string line;
  getline( ctlFile, line );
  std::istringstream lineStream( line );
  std::string nextFieldName;
  lineStream >> nextFieldName; // skip ID field

  while ( ! lineStream.eof() )
    {
    lineStream >> nextFieldName;

    if ( !SelectSet.empty() )
      // positive parameter selection by name
      {
      if ( SelectSet.find( nextFieldName ) == SelectSet.end() ) 
	{
	IgnoreSet.insert( FieldNames.size() );
	}
      }

    FieldNames.push_back( nextFieldName );
    }
  
  FieldNames.push_back( "CONST" );

  cmtk::DebugOutput( 2 ) << "\n\nImporting image files.\n";
  
  nParametersTotal = 0;
  while ( ! ctlFile.eof() ) 
    {
    getline( ctlFile, line );
    if ( line.length() && line[0] != '#' ) 
      {
      std::istringstream lineStream( line );

      std::string imgName;
      lineStream >> imgName;

      char imgPath[PATH_MAX];
      if ( snprintf( imgPath, sizeof( imgPath ), imgFilePatt, imgName.c_str() ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	throw cmtk::ExitException( 1 );
	}

      cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( imgPath ) );
      if ( !volume ) 
	{
	cmtk::StdErr << "ERROR: Could not read image file " << imgPath << "\n";
	throw cmtk::ExitException( 1 );
	}
      
      if ( CropImages )
	{
	volume->CropRegion() = CropImagesRegion;
	volume = cmtk::UniformVolume::SmartPtr( volume->GetCroppedVolume() );
	}
      
      cmtk::TypedArray::SmartPtr imageData( volume->GetData() );
      if ( !imageData ) 
	{
	cmtk::StdErr << "ERROR: no image data in file " << imgPath << "\n";
	throw cmtk::ExitException( 1 );
	}
      if ( ExponentialModel ) 
	imageData->ApplyFunctionDouble( cmtk::Wrappers::Log );
      
      if ( ! refVolume ) refVolume = volume;
      
      imagesData.push_back( imageData );

      cmtk::Types::DataItem param;
      if ( nParametersTotal ) 
	{
	for ( size_t p = 0; p < nParametersTotal; ++p ) 
	  {
	  if ( lineStream.eof() )
	    {
	    cmtk::StdErr << "ERROR: insufficient number of model parameters in line '" << line << "'\n";
	    throw cmtk::ExitException( 1 );
	    }
	  lineStream >> param;
	  // parameter exclusion by index
	  if ( IgnoreSet.find( p ) == IgnoreSet.end() ) 
	    {
	    vParam.push_back( param );
	    }
	  }
	if ( ! ExcludeConstant )
	  vParam.push_back( 1.0 );
	} 
      else 
	{
	while ( ! lineStream.eof() ) 
	  {
	  lineStream >> param;
	  if ( IgnoreSet.find( nParametersTotal ) == IgnoreSet.end() )
	    {
	    vParam.push_back( param );
	    }
	  else
	    {
	    cmtk::DebugOutput( 1 ) << "INFORMATION: Ignoring parameter #" << nParametersTotal << "\n";
	    }
	  ++nParametersTotal;
	  }
	if ( ! ExcludeConstant )
	  vParam.push_back( 1.0 );
	    
	// set nParameters from parameter array to account for ignored 
	// parameters.
	nParameters = vParam.size();
	cmtk::DebugOutput( 1 ) << "NParameters = " << nParameters << "\n";
	}
      }
    }
  
  if ( nParameters * imagesData.size() != vParam.size() ) 
    {
    cmtk::StdErr << "ERROR: number of parameters does not equal expected number"
		 << "(" << nParameters * imagesData.size() << " != " << vParam.size() << ")\n";
    throw cmtk::ExitException( 1 );
    }
  
  parameters = cmtk::Memory::AllocateArray<cmtk::Types::DataItem>( vParam.size() );
  for ( size_t i = 0; i < vParam.size(); ++i )
    parameters[i] = vParam[i];
  
  if ( cmtk::DebugOutput::GetGlobalLevel() > 0 )
    {
    cmtk::DebugOutput( 1 ) << "\n\nDesign Matrix: \n";
    for ( size_t j = 0; j < FieldNames.size(); ++j )
      {
      if ( IgnoreSet.find( j ) != IgnoreSet.end() ) continue;
      cmtk::DebugOutput( 1 ).GetStream().printf( "%8s\t", FieldNames[j].c_str() );
      }
    cmtk::DebugOutput( 1 ) << "\n";
    
    size_t ofs = 0;
    for ( size_t i = 0; i < imagesData.size(); ++i )
      {
      for ( size_t j = 0; j < nParameters; ++j, ++ofs )
	{
	cmtk::DebugOutput( 1 ).GetStream().printf( "%8.2f\t", parameters[ofs] );
	}
      cmtk::DebugOutput( 1 ) << "\n";
      }
    }
  
  if ( ! ExcludeConstant )
    ++nParametersTotal;
}

int
doMain( const int argc, const char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "General Linear Model" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Statistical modeling of pixel intensities in multiple images using a General Linear Model." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] ctlfile imgfile_pattern [ctlfile imgfile_pattern ...]" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Statistics and Modeling" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'x', "exclude-constant" ), &ExcludeConstant, true, "Exclude automatic constant parameter from model." );
    cl.AddSwitch( Key( 'n', "normalize" ), &NormalizeParameters, true, "Normalize model parameters w.r.t. data variances." );
    cl.AddSwitch( Key( 'e', "exp" ), &ExponentialModel, true, "Use exponential model rather than linear model." );
      
    cl.AddCallback( Key( 'i', "ignore-parameter" ), CallbackIgnore, "Ignore parameter with given NUMBER (0..n-1). Can be repeated." );
    cl.AddCallback( Key( 's', "select-parameter" ), CallbackSelect, "Select parameter with given NAME for model. Can be repeated." );
    cl.AddCallback( Key( 'c', "crop" ), CallbackCropImages, "To save space/time, crop images: x0,y0,z0,x1,y1,z2" );
    
    cl.AddOption( Key( 'O', "output-pattern" ), &OutputFilePatt, "Filename pattern for output (default: 'model_%s_%02d_%s.nii') with %d for parameter number" );

    cl.Parse( argc, argv );

    CtlFileName.push_back( cl.GetNext() );
    ImgFilePatt.push_back( cl.GetNext() );

    const char* next = cl.GetNextOptional();
    while ( next )
      {
      CtlFileName.push_back( next );
      ImgFilePatt.push_back( cl.GetNext() );      
      next = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e;
    }

  if ( CtlFileName.empty() || ImgFilePatt.empty() ) 
    {
    cmtk::StdErr << "ERROR: need both a control file name and an image file pattern\n";
    throw cmtk::ExitException( 1 );
    }
  
  cmtk::ProgressConsole progressIndicator;
  
  if ( cmtk::DebugOutput::GetGlobalLevel() > 0 )
    {
    std::set<std::string>::const_iterator it = SelectSet.begin();
    if ( it != SelectSet.end() )
      {
      cmtk::DebugOutput( 1 ) << "Selected model parameters:" << "\n";
      while ( it != SelectSet.end() )
	{
        cmtk::DebugOutput( 1 ) << "\t" << *it;
	++it;
	}
      cmtk::DebugOutput( 1 ) << "\n";
      }
    }

  cmtk::UniformVolume::SmartPtr refVolume( NULL );
  std::vector<size_t> nParameters( CtlFileName.size(), 0 );
  std::vector<size_t> nParametersTotal( CtlFileName.size(), 0 );
  std::vector<cmtk::Types::DataItem*> parameters( CtlFileName.size(), NULL );
  std::vector< std::vector<cmtk::TypedArray::SmartPtr> > ImagesData;

  std::vector<cmtk::TypedArray::SmartPtr> allImages;

  std::list<const char*>::const_iterator itCtlFile = CtlFileName.begin(), itImgFile = ImgFilePatt.begin();
  for ( size_t idx = 0; (itCtlFile != CtlFileName.end()) && (itImgFile != ImgFilePatt.end()); ++itCtlFile, ++itImgFile, ++idx )
    {
    std::vector<cmtk::TypedArray::SmartPtr> imagesDataNext;
    Import( *itCtlFile, *itImgFile, refVolume, imagesDataNext, nParameters[idx], parameters[idx], nParametersTotal[idx] );
    ImagesData.push_back( imagesDataNext );    
    
    if ( nParametersTotal[idx] != nParametersTotal[0] )
      {
      cmtk::StdErr << "Total number of parameters for control file #" << idx << "does not match that for file #0\n";
      throw cmtk::ExitException( 1 );
      }
    if ( nParameters[idx] != nParameters[0] )
      {
      cmtk::StdErr << "Number of active parameters for control file #" << idx << "does not match that for file #0\n";
      throw cmtk::ExitException( 1 );
      }

    const size_t totalNumberOfImages = allImages.size();
    allImages.resize( totalNumberOfImages + ImagesData[idx].size() );
    std::copy( ImagesData[idx].begin(), ImagesData[idx].end(), &allImages[totalNumberOfImages] );
    }
  
  double* allParameters = cmtk::Memory::AllocateArray<double>( nParameters[0] * allImages.size() );
  size_t idx = 0;
  for ( size_t ctl = 0; ctl < parameters.size(); ++ctl )
    {
    for ( size_t param = 0; param < ImagesData[ctl].size() * nParameters[0]; ++param, ++idx )
      {
      allParameters[idx] = parameters[ctl][param];
      }
    }
  
  cmtk::GeneralLinearModel glm( nParameters[0], allImages.size(), allParameters );
  
  if ( cmtk::DebugOutput::GetGlobalLevel() > 0 )
    {
    cmtk::DebugOutput( 1 ) << "\n\nSingular values:\n";
    size_t p = 0;
    for ( size_t pp = 0; pp < nParametersTotal[0]; ++pp ) 
      {
      // if this parameter is ignored, continue with next one.
      if ( IgnoreSet.find( pp ) != IgnoreSet.end() ) 
	continue;
      cmtk::DebugOutput( 1 ) << "\t" << glm.GetSingularValue( p++ );
      }
    cmtk::DebugOutput( 1 ) << "\n";
  
    cmtk::DebugOutput( 1 ) << "\n\nParameter correlation matrix:\n";
    cmtk::Matrix2D<double>* cc = glm.GetCorrelationMatrix();
    for ( size_t p = 0; p < nParameters[0]; ++p ) 
      {
      for ( size_t pp = 0; pp < nParameters[0]; ++pp ) 
	{
	cmtk::DebugOutput( 1 ).GetStream().printf( "%.2f\t", (*cc)[p][pp] );
	}
      cmtk::DebugOutput( 1 ) << "\n";
      }
    
    delete cc;
    }
  
  glm.FitModel( allImages, NormalizeParameters );

  if ( DumpStatistics ) 
    {
    cmtk::DebugOutput( 1 ) << "\n\nParameter normalization factors:\n";
    
    size_t p = 0;
    for ( size_t pp = 0; pp < nParametersTotal[0]; ++pp ) 
      {	  
      // if this parameter is ignored, continue with next one.
      if ( IgnoreSet.find( pp ) != IgnoreSet.end() ) continue;

      cmtk::StdOut.printf( "%d\t%f\n", (int)pp, glm.GetNormFactor( p ) );
      
      // increment actual parameter index.
      ++p;
      }
    }
  
  if ( OutputFilePatt ) 
    {
    cmtk::DebugOutput( 2 ) << "\n\nWriting output image files.\n";

    char outFileName[PATH_MAX];

    cmtk::TypedArray::SmartPtr fstatData = glm.GetFStat();
    if ( snprintf( outFileName, PATH_MAX, OutputFilePatt, "fstat", 0, "model" ) > PATH_MAX )
      {
      cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
      }
    else
      {
      refVolume->SetData( fstatData );
      cmtk::VolumeIO::Write( *refVolume, outFileName );
      }
      
    size_t p = 0;
    for ( size_t pp = 0; pp < nParametersTotal[0]; ++pp ) 
      {	  
      // if this parameter is ignored, continue with next one.
      if ( IgnoreSet.find( pp ) != IgnoreSet.end() ) continue;

      cmtk::TypedArray::SmartPtr modelData = glm.GetModel( p );
      if ( snprintf( outFileName, sizeof( outFileName ), OutputFilePatt, "param", pp, FieldNames[pp].c_str() ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	refVolume->SetData( modelData );
	cmtk::VolumeIO::Write( *refVolume, outFileName );
	}
	  
      cmtk::TypedArray::SmartPtr modelTStat = glm.GetTStat( p );
      if ( snprintf( outFileName, PATH_MAX, OutputFilePatt, "tstat", pp, FieldNames[pp].c_str() ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	refVolume->SetData( modelTStat );
	cmtk::VolumeIO::Write( *refVolume, outFileName );
	}
      
      // increment actual parameter index.
      ++p;
      }
    }
    
  return 0;
}

#include "cmtkSafeMain"
