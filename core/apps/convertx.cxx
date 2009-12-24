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
#include <cmtkTimers.h>

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>

#include <cmtkStudyList.h>
#include <cmtkClassStreamStudyList.h>

#include <math.h>
#include <cmtkMathFunctionWrappers.h>

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <set>
#include <map>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace convert
{
#endif

bool Verbose = false;

/// Image operation base class.
class ImageOperation
{
public:
  /// This class.
  typedef ImageOperation Self;

  /// Smart pointer.
  typedef cmtk::SmartPointer<Self> SmartPtr;
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr Apply( cmtk::UniformVolume::SmartPtr& volume ) 
  {
    return volume;
  }
};

/// List of image operations.
std::list<ImageOperation::SmartPtr> ImageOperationList;

/// Image operation: flip.
class ImageOperationConvertType
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationConvertType( const cmtk::ScalarDataType newType ) : m_NewType( newType ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {    
    switch ( this->m_NewType ) 
      {
      case cmtk::TYPE_CHAR:
      case cmtk::TYPE_BYTE:
      case cmtk::TYPE_SHORT:
      case cmtk::TYPE_USHORT:
      case cmtk::TYPE_INT:
      case cmtk::TYPE_FLOAT:
      case cmtk::TYPE_DOUBLE:
	if ( this->m_NewType != volume->GetData()->GetType() ) 
	  {
	  if ( Verbose )
	    cmtk::StdErr << "Converting to new data type.\n";
	  
	  volume->SetData( cmtk::TypedArray::SmartPtr( volume->GetData()->Convert( this->m_NewType ) ) );
	  }
	break;
      default:
	break;
      }
    return volume;
  }

  /// Create object to convert to "char" data.
  static void NewChar()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationConvertType( cmtk::TYPE_CHAR ) ) );
  }

  /// Create object to convert to "byte" data.
  static void NewByte()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationConvertType( cmtk::TYPE_BYTE ) ) );
  }

  /// Create object to convert to "short" data.
  static void NewShort()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationConvertType( cmtk::TYPE_SHORT ) ) );
  }

  /// Create object to convert to "unsigned short" data.
  static void NewUShort()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationConvertType( cmtk::TYPE_USHORT ) ) );
  }

  /// Create object to convert to "int" data.
  static void NewInt()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationConvertType( cmtk::TYPE_INT ) ) );
  }

  /// Create object to convert to "float" data.
  static void NewFloat()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationConvertType( cmtk::TYPE_FLOAT ) ) );
  }

  /// Create object to convert to "double" data.
  static void NewDouble()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationConvertType( cmtk::TYPE_DOUBLE ) ) );
  }
  
private:
  /// New data type.
  cmtk::ScalarDataType m_NewType;
};

/// Image operation: flip.
class ImageOperationFlip
  /// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationFlip( const int normalAxis ) : m_NormalAxis( normalAxis ) {}

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->ApplyMirrorPlane( this->m_NormalAxis );
    return volume;
  }

  /// Create x flip object.
  static void NewX()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationFlip( cmtk::AXIS_X ) ) );
  }

  /// Create y flip object.
  static void NewY()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationFlip( cmtk::AXIS_Y ) ) );
  }

  /// Create y flip object.
  static void NewZ()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationFlip( cmtk::AXIS_Z ) ) );
  }

private:
  /// The normal axis of the flip.
  int m_NormalAxis;
};

/// Apply mask image.
class ImageOperationApplyMask
/// Inherit generic image operation.
  : public ImageOperation
{
public:
  /// Constructor.
  ImageOperationApplyMask( const cmtk::UniformVolume::SmartPtr& maskVolume ) : m_MaskVolume( maskVolume ) {}

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    const std::string maskOrientation = this->m_MaskVolume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
    const std::string workingOrientation = volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
    if ( maskOrientation != workingOrientation )
      {
      if ( Verbose )
	{
	cmtk::StdErr << "INFO: reorienting mask from '" << maskOrientation << "' to '" << workingOrientation << "' to match working image.\n";
	}
      this->m_MaskVolume = cmtk::UniformVolume::SmartPtr( this->m_MaskVolume->GetReoriented( workingOrientation.c_str() ) );
      }
    
    for ( int dim = 0; dim < 3; ++dim )
      {
      if ( this->m_MaskVolume->m_Dims[dim] != volume->m_Dims[dim] )
	{
	cmtk::StdErr << "ERROR: mask volume dimensions do not match working volume dimensions.\n";
	exit( 1 );
	}
      }
    
    const cmtk::TypedArray* maskData = this->m_MaskVolume->GetData();
    cmtk::TypedArray::SmartPtr& volumeData = volume->GetData();

    const size_t nPixels = volume->GetNumberOfPixels();
    for ( size_t i = 0; i < nPixels; ++i )
      if ( maskData->IsPaddingOrZeroAt( i ) ) 
	volumeData->SetPaddingAt( i );
    return volume;
  }

  /// Create new mask operation.
  static void New( const char* maskFileName )
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationApplyMask( ReadMaskFile( maskFileName ) ) ) );
  }

  /// Create new inverse mask operation.
  static void NewInverse( const char* maskFileName )
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationApplyMask( ReadMaskFile( maskFileName, true /*inverse*/ ) ) ) );
  }

private:
  /// The mask volume.
  cmtk::UniformVolume::SmartPtr m_MaskVolume;

  /// Read the actual mask file.
  static cmtk::UniformVolume::SmartPtr ReadMaskFile( const char* maskFileName, const bool inverse = false )
  {
    cmtk::UniformVolume::SmartPtr maskVolume( cmtk::VolumeIO::ReadOriented( maskFileName, Verbose ) );
    if ( !maskVolume || !maskVolume->GetData() ) 
      {
      cmtk::StdErr << "ERROR: could not read mask from file " << maskFileName << "\nProgram will terminate now, just to be safe.\n";
      exit( 1 );
      }
    
    // binarize mask to 1/0, convert to char, and also consider "inverse" flag in the process.
    cmtk::TypedArray::SmartPtr& maskData = maskVolume->GetData();
    const size_t nPixels = maskData->GetDataSize();
    for ( size_t n = 0; n < nPixels; ++n )
      {
      if ( maskData->IsPaddingOrZeroAt( n ) != inverse ) 
	maskData->Set( n, 1 );
      else
	maskData->Set( n, 0 );
      }
    maskVolume->SetData( cmtk::TypedArray::SmartPtr( maskData->Convert( cmtk::TYPE_BYTE ) ) );
    
    return maskVolume;
  }
};

/// Image operation: erode or dilate.
class ImageOperationErodeDilate
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationErodeDilate( const int iterations ) : m_Iterations( iterations ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    if ( this->m_Iterations < 0 )
      volume->ApplyErode( -this->m_Iterations );
    else
      if ( this->m_Iterations > 0 )
	volume->ApplyDilate( this->m_Iterations );
    return volume;
  }

  /// Create new dilation operation.
  static void NewDilate( const long int iterations )
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationErodeDilate( iterations ) ) );
  }

  /// Create new erosion operation.
  static void NewErode( const long int iterations )
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationErodeDilate( -iterations ) ) );
  }
  
private:
  /// Number of iterations of erosion (if negative) or dilation (if positive).
  int m_Iterations;
};

/// Image operation: create binary or multi-valued boundary map.
class ImageOperationBoundaryMap
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationBoundaryMap( const bool multiValued ) : m_MultiValued( multiValued ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->SetData( cmtk::TypedArray::SmartPtr( volume->GetBoundaryMap( this->m_MultiValued ) ) );
    return volume;
  }
  
  /// Create new binary boundary map operation.
  static void New()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationBoundaryMap( false ) ) );
  }

  /// Create new multi-valued boundary map operation.
  static void NewMulti()
  {
    ImageOperationList.push_back( SmartPtr( new ImageOperationBoundaryMap( true ) ) );
  }

private:
  /// Multi-valued flag: if this is set, a multi-valued boundary map will be created, otherwise a binary map.
  bool m_MultiValued;
};

/// Image operation: grid downsampling.
class ImageOperationDownsample
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationDownsample( const int factorX, const int factorY, const int factorZ ) : m_FactorX( factorX ), m_FactorY( factorY ), m_FactorZ( factorZ ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    const int factors[3] = { this->m_FactorX, this->m_FactorY, this->m_FactorZ };
    return cmtk::UniformVolume::SmartPtr( volume->GetDownsampled( factors ) );
  }
  
  /// Create a new downsampler.
  static void New( const char* arg )
  {
    int factorsX = 1;
    int factorsY = 1;
    int factorsZ = 1;

    const size_t nFactors = sscanf( arg, "%d,%d,%d", &factorsX, &factorsY, &factorsZ );
    if ( nFactors == 1 )
      {
      factorsZ = factorsY = factorsX;
      }
    else
      {
      if ( nFactors != 3 )
	{
	cmtk::StdErr << "ERROR: downsampling factors must either be three integers, x,y,z, or a single integer\n";
	exit( 1 );
	}
      }
    ImageOperationList.push_back( SmartPtr( new ImageOperationDownsample( factorsX, factorsY, factorsZ ) ) );
  }
  
private:
  /// Downsampling factor in X direction.
  int m_FactorX;

  /// Downsampling factor in Y direction.
  int m_FactorY;

  /// Downsampling factor in Z direction.
  int m_FactorZ;
};

/// Image operation: grid downsampling.
class ImageOperationMedianFilter
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationMedianFilter( const int radiusX, const int radiusY, const int radiusZ ) : m_RadiusX( radiusX ), m_RadiusY( radiusY ), m_RadiusZ( radiusZ ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->ApplyMedianFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ );
    return volume;
  }
  
  /// Create a new downsampler.
  static void New( const char* arg )
  {
    int radiusX = 1;
    int radiusY = 1;
    int radiusZ = 1;

    const size_t nRadii = sscanf( arg, "%d,%d,%d", &radiusX, &radiusY, &radiusZ );
    if ( nRadii == 1 )
      {
      radiusZ = radiusY = radiusX;
      }
    else
      {
      if ( nRadii != 3 )
	{
	cmtk::StdErr << "ERROR: downsampling radii must either be three integers, x,y,z, or a single integer\n";
	exit( 1 );
	}
      }
    ImageOperationList.push_back( SmartPtr( new ImageOperationMedianFilter( radiusX, radiusY, radiusZ ) ) );
  }
  
private:
  /// Downsampling radius in X direction.
  int m_RadiusX;

  /// Downsampling radius in Y direction.
  int m_RadiusY;

  /// Downsampling radius in Z direction.
  int m_RadiusZ;
};

int
main( int argc, char* argv[] )
{
  cmtk::Types::DataItem paddingDataValue = 0.0;
  bool paddingDataFlag = false;

  const char* imagePathIn = NULL;
  const char* imagePathOut = NULL;

  try 
    {
    cmtk::CommandLine cl( argc, argv );  
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Convert between image file formats and data types. Also apply simple, general-purpose image operations in the process." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] infile outfile" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.BeginGroup( "Input", "Input Image Controls" );
    cl.AddOption( Key( 'N', "set-padding" ), &paddingDataValue, "Set Padding data for input image.", &paddingDataFlag );
    cl.EndGroup();

    cl.BeginGroup( "Conversion", "Data Type Conversion" );
    cl.AddCallback( Key( 'c', "char" ), &ImageOperationConvertType::NewChar, "8 bits, signed integer" );
    cl.AddCallback( Key( 'b', "byte" ), &ImageOperationConvertType::NewByte, "8 bits, unsigned integer" );
    cl.AddCallback( Key( 's', "short" ), &ImageOperationConvertType::NewShort, "16 bits, signed integer" );
    cl.AddCallback( Key( 'u', "ushort" ), &ImageOperationConvertType::NewUShort, "16 bits, unsigned integer" );
    cl.AddCallback( Key( 'i', "int" ), &ImageOperationConvertType::NewInt, "32 bits signed integer" );
    cl.AddCallback( Key( 'f', "float" ), &ImageOperationConvertType::NewFloat, "32 bits floating point" );
    cl.AddCallback( Key( 'd', "double" ), &ImageOperationConvertType::NewDouble, "64 bits floating point\n" );
    cl.EndGroup();
    
    cl.BeginGroup( "Flipping", "Image Flipping" );
    cl.AddCallback( Key( "flip-x" ), &ImageOperationFlip::NewX, "Flip (mirror) along x-direction" );
    cl.AddCallback( Key( "flip-y" ), &ImageOperationFlip::NewY, "Flip (mirror) along y-direction" );
    cl.AddCallback( Key( "flip-z" ), &ImageOperationFlip::NewZ, "Flip (mirror) along z-direction" );
    cl.EndGroup();

    cl.BeginGroup( "Masking", "Image Masking" );
    cl.AddCallback( Key( 'M', "mask" ), &ImageOperationApplyMask::New, "Binary mask file name: eliminate all image pixels where mask is 0." );
    cl.AddCallback( Key( "mask-inverse" ), &ImageOperationApplyMask::NewInverse, "Inverse binary mask file name eliminate all image pixels where mask is NOT 0." );
    cl.EndGroup();

    cl.BeginGroup( "Morphological", "Morphological Operations" );
    cl.AddCallback( Key( "erode" ), &ImageOperationErodeDilate::NewErode, "Morphological erosion operator" );
    cl.AddCallback( Key( "dilate" ), &ImageOperationErodeDilate::NewDilate, "Morphological dilation operator" );
    cl.AddCallback( Key( "boundary-map" ), &ImageOperationBoundaryMap::New, "Create boundary map" );
    cl.AddCallback( Key( "multi-boundary-map" ), &ImageOperationBoundaryMap::NewMulti, "Create multi-valued boundary map" );
    cl.EndGroup();

    cl.BeginGroup( "Filtering", "Filter Operations" );
    cl.AddCallback( Key( "median-filter" ), &ImageOperationMedianFilter::New, "Median filter. This operation takes the filter radius as the parameter. "
		    "A single integer defines the kernel radius in all three dimensions. Three comma-separated integers define separate radii for the three dimensions." );
    cl.EndGroup();

    cl.AddCallback( Key( "downsample" ), &ImageOperationDownsample::New, "Downsample image by factors 'x,y,z' or by single factor 'xyz'" );

    cl.AddParameter( &imagePathIn, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &imagePathOut, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );

    if ( ! cl.Parse() ) 
      return 1;
    }
  catch ( cmtk::CommandLine::Exception ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( imagePathIn, Verbose ) );
  if ( ! volume ) 
    {
    cmtk::StdErr << "ERROR: could not read image " << imagePathIn << "\n";
    exit( 1 );
    }
  else
    {
    cmtk::TypedArray::SmartPtr volumeData = volume->GetData();
    if ( ! volumeData ) 
      {
      cmtk::StdErr << "ERROR: image seems to contain no data.\n";
      exit( 1 );
      }
    }
  
  if ( paddingDataFlag ) 
    {
    volume->GetData()->SetPaddingValue( paddingDataValue );
    }
  
  for ( std::list<ImageOperation::SmartPtr>::iterator opIt = ImageOperationList.begin(); opIt != ImageOperationList.end(); ++opIt )
    {
    volume = (*opIt)->Apply( volume );
    }

  if ( Verbose )
    {
    cmtk::StdErr << "Writing to file " << imagePathOut << "\n";
    }
  
  cmtk::VolumeIO::Write( volume, imagePathOut, Verbose );
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace convert
} // namespace apps
} // namespace cmtk
#endif
