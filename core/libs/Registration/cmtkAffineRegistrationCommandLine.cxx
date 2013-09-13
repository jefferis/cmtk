/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkAffineRegistrationCommandLine.h"

#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkThreads.h>
#include <System/cmtkTimers.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkCompressedStream.h>
#include <System/cmtkMountPoints.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkAnatomicalOrientation.h>
#include <Base/cmtkTransformChangeToSpaceAffine.h>
#include <Base/cmtkTransformChangeFromSpaceAffine.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkClassStreamAffineXform.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkAffineXformITKIO.h>

#include <Registration/cmtkRegistrationCallback.h>
#include <Registration/cmtkProtocolCallback.h>

#ifdef CMTK_USE_SQLITE
#  include <Registration/cmtkImageXformDB.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef HAVE_SYS_UTSNAME_H
#  include <sys/utsname.h>
#endif

#include <iostream>
#include <stdio.h>
#include <string.h>

#ifdef _MSC_VER
#  include <direct.h>
#endif

namespace 
cmtk
{

/** \addtogroup Registration */
//@{

AffineRegistrationCommandLine
::AffineRegistrationCommandLine 
( const int argc, const char* argv[] ) 
{
  this->m_Metric = 0;

  this->m_AutoMultiLevels = 0;
  CoarsestResolution = -1;
  this->m_Exploration = 8;
  this->m_Accuracy = 0.1;
  this->m_Sampling = 1.0;
  InitXlate = 0;
  this->m_NoSwitch = 0;

  std::string InitialStudylist;
  std::string clArg1; // input studylist or reference image
  std::string clArg2; // empty or floating image

  try 
    {
    CommandLine cl( CommandLine::PROPS_XML );
    cl.SetProgramInfo( CommandLine::PRG_TITLE, "Rigid and affine registration" );
    cl.SetProgramInfo( CommandLine::PRG_DESCR, "This program performs rigid and affine image registration using multi-resolution optimization of voxel-based image similarity measures." );
    cl.SetProgramInfo( CommandLine::PRG_CATEG, "CMTK.Registration" );

    typedef CommandLine::Key Key;
    cl.BeginGroup( "Automation", "Automation Options" );
    cl.AddOption( Key( "auto-multi-levels" ), &this->m_AutoMultiLevels, "Automatic optimization and resolution parameter generation for <n> levels" );

    cl.BeginGroup( "Optimization", "Optimization settings" );
    cl.AddOption( Key( 'e', "exploration" ), &this->m_Exploration, "Exploration [initial optimizer step size]" );
    cl.AddOption( Key( 'a', "accuracy" ), &this->m_Accuracy, "Accuracy [final optimizer step size]" );
    cl.AddOption( Key( 'f', "stepfactor" ), &this->OptimizerStepFactor, "Factor for search step size reduction. Must be > 0.0 and < 1.0" );
    cl.AddOption( Key( "delta-f-threshold" ), &this->m_DeltaFThreshold, "Optional threshold to terminate optimization (level) if relative change of target function drops below this value." );
    cl.EndGroup();

    cl.BeginGroup( "Resolution", "Image resolution parameters" );
    cl.AddOption( Key( 's', "sampling" ), &this->m_Sampling, "Image sampling (finest resampled image resolution)" );
    cl.AddOption( Key( "coarsest" ), &this->CoarsestResolution, "Upper limit for image sampling in multiresolution hierarchy" );

    cl.AddSwitch( Key( "omit-original-data" ), &this->m_UseOriginalData, false, "Do not use original data in full resolution, omit final stage in multiresolution hierarchy, thus reducing computation time." );
    cl.EndGroup();

    cl.BeginGroup( "Transformation", "Transformation parameters" );
    cl.AddVector( Key( "dofs" ), this->NumberDOFs, "Add number of degrees of freedom [can be repeated]" );
    cl.AddVector( Key( "dofs-final" ), this->NumberDOFsFinal, "Add number of degrees of freedom for final level only [can be repeated]" );
    
    cl.AddSwitch( Key( 'n', "no-switch" ), &this->m_NoSwitch, true, "Do not auto-switch reference and floating image for improved computational performance" );
    cl.AddSwitch( Key( 'i', "initxlate" ), &InitXlate, true, "Initialized transformation by translating floating image FOV center onto reference image FOV center" );

    cl.AddOption( Key( "initial" ), &InitialStudylist, "Initialize transformation from given path" )->SetProperties( CommandLine::PROPS_XFORM );
    cl.AddSwitch( Key( "initial-is-inverse" ), &this->m_InitialXformIsInverse, true, "Invert initial transformation before initializing registration" );
    cl.EndGroup();

    cl.BeginGroup( "Image data", "Image data" );
    CommandLine::EnumGroup<int>::SmartPtr
      metricGroup = cl.AddEnum( "registration-metric", &this->m_Metric, "Registration metric for motion estimation by image-to-image registration." );
    metricGroup->AddSwitch( Key( "nmi" ), 0, "Normalized Mutual Information metric" );
    metricGroup->AddSwitch( Key( "mi" ), 1, "Standard Mutual Information metric" );
    metricGroup->AddSwitch( Key( "cr" ), 2, "Correlation Ratio metric" );
    metricGroup->AddSwitch( Key( "msd" ), 4, "Mean Squared Difference metric" );
    metricGroup->AddSwitch( Key( "ncc" ), 5, "Normalized Cross Correlation metric" );

    cl.AddSwitch( Key( "match-histograms" ), &this->m_MatchFltToRefHistogram, true, "Match floating image histogram to reference image histogram." );

    this->m_PreprocessorRef.AttachToCommandLine( cl );
    this->m_PreprocessorFlt.AttachToCommandLine( cl );

    cl.BeginGroup( "Output", "Output parameters" )->SetProperties( CommandLine::PROPS_NOXML );
    cl.AddOption( Key( 'o', "outlist" ), &this->Studylist, "Output path for final transformation" );
    cl.AddOption( Key( "out-matrix" ), &this->OutMatrixName, "Output path for final transformation in matrix format" );
    cl.AddOption( Key( "out-parameters" ), &this->OutParametersName, "Output path for final transformation in plain parameter list format" );
    cl.AddOption( Key( 'p', "protocol" ), &this->Protocol, "Optimization protocol output file name" );
    cl.AddOption( Key( 't', "time" ), &this->Time, "Computation time statistics output file name" );
    cl.EndGroup();

    cl.BeginGroup( "SlicerImport", "Import Results into Slicer" );
    cl.AddOption( Key( "out-itk" ), &this->m_OutputPathITK, "Output path for final transformation in ITK format" )
      ->SetProperties( CommandLine::PROPS_XFORM | CommandLine::PROPS_OUTPUT )
      ->SetAttribute( "reference", "FloatingImage" );
    cl.AddOption( Key( "write-reformatted" ), &this->m_ReformattedImagePath, "Write reformatted floating image." )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OUTPUT );
    cl.EndGroup();
    
#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &this->m_UpdateDB, "Path to image/transformation database that should be updated with the new registration and/or reformatted image." );
    cl.EndGroup();
#endif

    cl.AddParameter( &clArg1, "ReferenceImage", "Reference (fixed) image path" )->SetProperties( CommandLine::PROPS_IMAGE );
    cl.AddParameter( &clArg2, "FloatingImage", "Floating (moving) image path" )->SetProperties( CommandLine::PROPS_IMAGE | CommandLine::PROPS_OPTIONAL );

    cl.Parse( argc, argv );
    }
  catch ( const CommandLine::Exception& ex )
    {
    StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  if ( (OptimizerStepFactor <= 0) || (OptimizerStepFactor >= 1) ) 
    {
    StdErr << "ERROR: step factor value " << OptimizerStepFactor << " is invalid. Must be in range (0..1)\n";
    throw cmtk::ExitException( 1 );
    }

  this->SetInitialTransformation( AffineXform::SmartPtr( new AffineXform() ) );
    
  if ( ! clArg2.empty() ) 
    {
    Study1 = clArg1;
    Study2 = clArg2;
    } 
  else
    {
    this->m_InitialXformPath = clArg1;

    if ( ! InitialStudylist.empty() ) 
      {
      StdErr << "WARNING: transformation of input studylist will be overriden by transformation provided with '--initial'.\n";
      }
    
    DebugOutput( 1 ) << "Reading input studylist " << this->m_InitialXformPath << ".\n";
    
    ClassStreamInput typedStream( MountPoints::Translate(this->m_InitialXformPath), "registration" );
    if ( ! typedStream.IsValid() ) 
      {
      StdErr << "ERROR: could not open studylist archive " << this->m_InitialXformPath << ".\n";
      throw cmtk::ExitException( 1 );
      }

    typedStream.Seek ( "registration" );
    Study1 = typedStream.ReadStdString( "reference_study" );
    Study2 = typedStream.ReadStdString( "floating_study" );
    if ( ! Study2.empty() )
      {
      AffineXform::SmartPtr affineXform;
      typedStream >> affineXform;
      this->SetInitialTransformation( affineXform );
      }
    else
      {
      // legacy studylists have inverse transformation in them
      Study2 = typedStream.ReadStdString( "model_study" );
      AffineXform::SmartPtr affineXform;
      typedStream >> affineXform;
      this->SetInitialTransformation( affineXform->GetInverse() );
      }

    typedStream.Close();
    }

  if ( Study1.empty() )
    {
    StdErr << "ERROR: reference image path resolved to NULL.\n";
    throw cmtk::ExitException( 1 );
    }
  
  if ( Study2.empty() )
    {
    StdErr << "ERROR: floating image path resolved to NULL.\n";
    throw cmtk::ExitException( 1 );
    }
  
  UniformVolume::SmartPtr volume( VolumeIO::ReadOriented( this->Study1 ) );
  if ( !volume )
    {
    StdErr << "ERROR: volume " << this->Study1 << " could not be read\n";
    throw cmtk::ExitException( 1 );
    }
  this->SetVolume_1( UniformVolume::SmartPtr( this->m_PreprocessorRef.GetProcessedImage( volume ) ) );

  volume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( this->Study2 ) );
  if ( !volume )
    {
    StdErr << "ERROR: volume " << this->Study2 << " could not be read\n";
    throw cmtk::ExitException( 1 );
    }
  this->SetVolume_2(  UniformVolume::SmartPtr( this->m_PreprocessorFlt.GetProcessedImage( volume ) ) );

  if ( ! InitialStudylist.empty() ) 
    {
    this->m_InitialXformPath = InitialStudylist;
    Xform::SmartPtr xform( XformIO::Read( InitialStudylist ) );
    if ( ! xform ) 
      {
      StdErr << "ERROR: could not read transformation from " << InitialStudylist << "\n";
      throw cmtk::ExitException( 1 );
      }
    
    AffineXform::SmartPtr affine( AffineXform::SmartPtr::DynamicCastFrom( xform ) );
    if ( ! affine )
      {
      StdErr << "ERROR: transformation " << InitialStudylist << " is not affine.\n";
      throw cmtk::ExitException( 1 );
      }

    if ( affine->GetMetaInfo( META_SPACE ) != AnatomicalOrientation::ORIENTATION_STANDARD )
      {
      TransformChangeFromSpaceAffine toStandardSpace( *affine, *(this->m_Volume_1), *(this->m_Volume_2) );
      *affine = toStandardSpace.GetTransformation();
      affine->SetMetaInfo( META_SPACE, AnatomicalOrientation::ORIENTATION_STANDARD );
      }
    
    this->SetInitialTransformation( affine );
    }
  
  if ( InitXlate ) 
    {
    if ( ! this->m_InitialXformPath.empty() ) 
      {
      StdErr << "WARNING: Initial transformation was taken from studylist. Switch --initxlate / -i will be ignored.\n";
      } 
    else
      {
      this->SetInitialAlignCenters();
      }
    }

  if ( this->m_AutoMultiLevels > 0 )
    {
    const Types::Coordinate minDelta = std::min( this->m_Volume_1->GetMinDelta(), this->m_Volume_2->GetMinDelta() );
    const Types::Coordinate maxDelta = std::max( this->m_Volume_1->GetMaxDelta(), this->m_Volume_2->GetMaxDelta() );

    this->m_Accuracy = 0.1 * minDelta;
    this->m_Sampling = maxDelta;
    this->m_Exploration = maxDelta * (1<<(this->m_AutoMultiLevels-1));
    }
  
  if ( !Protocol.empty() ) 
    {
    RegistrationCallback::SmartPtr callback( new ProtocolCallback( Protocol ) );
    this->SetCallback( callback );
    }
}

CallbackResult
AffineRegistrationCommandLine::InitRegistration ()
{
  CallbackResult Result = AffineRegistration::InitRegistration();
  return Result;
}
	
void
AffineRegistrationCommandLine::OutputResultMatrix( const std::string& matrixName ) const
{
  const AffineXform::MatrixType& matrix = this->GetTransformation()->Matrix;
  
  FILE* mfile = fopen( matrixName.c_str(), "w" );
  if ( mfile )
    {
    for ( int i = 0; i < 4; ++i )
      {
      fprintf( mfile, "%e\t%e\t%e\t%e\n", static_cast<float>( matrix[0][i] ), static_cast<float>( matrix[1][i] ), static_cast<float>( matrix[2][i] ), static_cast<float>( matrix[3][i] ) );
      }
    fclose( mfile );
    }
}

void
AffineRegistrationCommandLine::OutputResultParameters
( const std::string& paramsName, const CoordinateVector& v ) const
{
  FILE* pfile = fopen( paramsName.c_str(), "w" );
  if ( pfile )
    {
    for ( unsigned int idx=0; idx < v.Dim; ++idx )
      fprintf( pfile, "#%d: %f\n", idx, v.Elements[idx] );
    fclose( pfile );
    }
}

void
AffineRegistrationCommandLine::OutputResultList( const std::string& studyList ) const
{
  ClassStreamOutput classStream( studyList, "studylist", ClassStreamOutput::MODE_WRITE );
  if ( !classStream.IsValid() ) return;
  
  classStream.Begin( "studylist" );
  classStream.WriteInt( "num_sources", 2 );
  classStream.End();
    
  classStream.Begin( "source" );
  classStream.WriteString( "studyname", CompressedStream::GetBaseName( Study1 ) );
  classStream.End();
    
  classStream.Begin( "source" );
  classStream.WriteString( "studyname", CompressedStream::GetBaseName( Study2 ) );
  classStream.End();
    
  classStream.Close();
    
  classStream.Open( studyList, "registration", ClassStreamOutput::MODE_WRITE );
    
  classStream.Begin( "registration" );
  classStream.WriteString( "reference_study", CompressedStream::GetBaseName( Study1 ) );
  classStream.WriteString( "floating_study", CompressedStream::GetBaseName( Study2 ) );
    
  classStream << *(this->GetTransformation());
    
  classStream.End();
  classStream.Close();
    
  classStream.Open( studyList, "settings", ClassStreamOutput::MODE_WRITE );
  classStream.WriteDouble( "exploration", this->m_Exploration );
  classStream.WriteDouble( "accuracy", this->m_Accuracy );
  classStream.WriteDouble( "min_sampling", this->m_Sampling );
  classStream.WriteDouble( "coarsest_resolution", CoarsestResolution );
  classStream.WriteInt( "metric", this->m_Metric );
  classStream.WriteDouble( "optimizer_step_factor", OptimizerStepFactor );
  classStream.WriteBool( "no_switch", this->m_NoSwitch );

  this->m_PreprocessorRef.WriteSettings( classStream );  
  this->m_PreprocessorFlt.WriteSettings( classStream );  

  classStream.Close();
    
  classStream.Open( studyList, "statistics", ClassStreamOutput::MODE_WRITE );
  classStream.WriteDouble( "time", this->GetTotalElapsedTime() );
  classStream.WriteDouble( "walltime", this->GetTotalElapsedWalltime() );
#ifdef CMTK_USE_PTHREADS
  classStream.WriteDouble( "thread_time", this->GetThreadTotalElapsedTime() );
#endif
    
#ifndef _MSC_VER
  struct utsname name;
  if ( uname( &name ) >= 0 ) 
    {
    classStream.WriteString( "host", name.nodename );
    classStream.WriteString( "system", name.sysname );
    }
#endif
  classStream.Close();
}

void
AffineRegistrationCommandLine::OutputResult ( const CoordinateVector* v, const CallbackResult irq )
{
  DebugOutput( 1 ) << "Resulting transformation parameters: \n";
  for ( unsigned int idx=0; idx<v->Dim; ++idx )
    DebugOutput( 1 ).GetStream().printf( "#%d: %f\n", idx, v->Elements[idx] );
  
  if ( !this->OutMatrixName.empty() )
    {
    if ( irq != CALLBACK_OK )
      this->OutputResultMatrix( this->OutMatrixName + "-partial" );
    else
      this->OutputResultMatrix( this->OutMatrixName );
    }

  if ( !this->OutParametersName.empty() )
    {
    if ( irq != CALLBACK_OK )
      this->OutputResultParameters( this->OutParametersName + "-partial", *v );
    else
      this->OutputResultParameters( this->OutParametersName, *v );
    }

  if ( !this->Studylist.empty() ) 
    {
    if ( irq != CALLBACK_OK )
      this->OutputResultList( this->Studylist + "-partial" );
    else
      this->OutputResultList( this->Studylist );
    }

  if ( !this->m_OutputPathITK.empty() ) 
    {
    TransformChangeToSpaceAffine toNative( *(this->GetTransformation()), *(this->m_Volume_1), *(this->m_Volume_2), AnatomicalOrientationBase::SPACE_ITK );
    if ( irq != CALLBACK_OK )
      AffineXformITKIO::Write( this->m_OutputPathITK + "-partial", toNative.GetTransformation() );
    else
      AffineXformITKIO::Write( this->m_OutputPathITK, toNative.GetTransformation() );
    }
  
  if ( !this->m_ReformattedImagePath.empty() )
    {
    if ( irq != CALLBACK_OK )
      VolumeIO::Write( *(this->GetReformattedFloatingImage()), this->m_ReformattedImagePath + "-partial" );
    else
      VolumeIO::Write( *(this->GetReformattedFloatingImage()), this->m_ReformattedImagePath );
    }

#ifdef CMTK_USE_SQLITE
  if ( (irq == CALLBACK_OK) && !this->m_UpdateDB.empty() )
    {
    try
      {
      ImageXformDB db( this->m_UpdateDB );
      
      if ( !this->m_ReformattedImagePath.empty() )
	{
	db.AddImage( this->m_ReformattedImagePath, this->m_ReferenceVolume->GetMetaInfo( META_FS_PATH ) );
	}
      
      if ( !this->Studylist.empty() )
	{
	if ( ! this->m_InitialXformPath.empty() )
	  {
	  db.AddRefinedXform( this->Studylist, true /*invertible*/, this->m_InitialXformPath, this->m_InitialXformIsInverse );
	  }
	else
	  {
	  db.AddImagePairXform( this->Studylist, true /*invertible*/, this->m_ReferenceVolume->GetMetaInfo( META_FS_PATH ), this->m_FloatingVolume->GetMetaInfo( META_FS_PATH ) );
	  }
	}
      }
    catch ( const ImageXformDB::Exception& ex )
      {
      StdErr << "DB ERROR: " << ex.what() << " on database " << this->m_UpdateDB << "\n";
      }
    }
#endif
}

void
AffineRegistrationCommandLine::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f,
  const int index, const int total )
{
  DebugOutput( 1 ).GetStream().printf( "Entering resolution level %d out of %d...\n", index, total );
  this->Superclass::EnterResolution( v, f, index, total );
}

CallbackResult
AffineRegistrationCommandLine::Register ()
{
  const double baselineTime = Timers::GetTimeProcess();
  CallbackResult Result = Superclass::Register();
  const int elapsed = static_cast<int>( Timers::GetTimeProcess() - baselineTime );

  if ( !this->Time.empty() ) 
    {
    FILE *tfp = fopen( this->Time.c_str(), "w" );
    
    if ( tfp ) 
      {
      fprintf( tfp, "%d\n", elapsed );
      fclose( tfp );
      } 
    else 
      {
      StdErr << "Could not open time file " << this->Time << "\n";
      }
    }
  return Result;
}

} // namespace cmtk

