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

#include <cmtkAffineRegistrationCommandLine.h>

#include <cmtkConsole.h>
#include <cmtkThreads.h>
#include <cmtkTimers.h>
#include <cmtkCommandLine.h>

#include <cmtkTypes.h>
#include <cmtkRegistrationCallback.h>
#include <cmtkProtocolCallback.h>

#include <cmtkMountPoints.h>
#include <cmtkVolumeIO.h>
#include <cmtkAnatomicalOrientation.h>

#include <cmtkClassStream.h>
#include <cmtkClassStreamAffineXform.h>
#include <cmtkCompressedStream.h>
#include <cmtkXformIO.h>

#include <stdio.h>
#include <string.h>

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
( int argc, char* argv[] ) 
{
  this->m_Metric = 0;

  this->m_AutoMultiLevels = 0;
  CoarsestResolution = -1;
  this->m_Exploration = 8;
  this->m_Accuracy = 0.1;
  this->m_Sampling = 1.0;
  OutParametersName = OutMatrixName = Studylist = Protocol = Time = NULL;
  InitXlate = 0;
  this->m_NoSwitch = 0;

  Verbose = 0;

  const char* inStudylist = NULL;
  const char *InitialStudylist = NULL;
  Study1 = Study2 = NULL;

  const char* clArg1; // input studylist or reference image
  const char* clArg2; // empty or floating image

  try 
    {
    CommandLine cl( argc, argv );
    cl.SetProgramInfo( CommandLine::PRG_TITLE, "Rigid and affine registration" );
    cl.SetProgramInfo( CommandLine::PRG_DESCR, "This program performs rigid and affine image registration using multi-resolution optimization of voxel-based image similarity measures." );
    cl.SetProgramInfo( CommandLine::PRG_SYNTX, "[options] [refImage fltImage | initialStudylist]" );
    cl.SetProgramInfo( CommandLine::PRG_CATEG, "CMTK.Image Registration" );

    typedef CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );
    cl.AddSwitch( Key( 'q', "quiet" ), &Verbose, false, "Quiet mode" );

    cl.BeginGroup( "Automation", "Automation Options" );
    cl.AddOption( Key( "auto-multi-levels" ), &this->m_AutoMultiLevels, "Automatic optimization and resolution parameter generation for <n> levels" );

    cl.BeginGroup( "Optimization", "Optimization settings" );
    cl.AddOption( Key( 'e', "exploration" ), &this->m_Exploration, "Exploration [initial optimizer step size]" );
    cl.AddOption( Key( 'a', "accuracy" ), &this->m_Accuracy, "Accuracy [final optimizer step size]" );
    cl.AddOption( Key( 'f', "stepfactor" ), &this->OptimizerStepFactor, "Factor for search step size reduction. Must be > 0.0 and < 1.0 [default: 0.5]" );
    cl.EndGroup();

    cl.BeginGroup( "Resolution", "Image resolution parameters" );
    cl.AddOption( Key( 's', "sampling" ), &this->m_Sampling, "Image sampling (finest resampled image resolution)" );
    cl.AddOption( Key( "coarsest" ), &this->CoarsestResolution, "Upper limit for image sampling in multiresolution hierarchy" );

    cl.AddSwitch( Key( "use-original-data" ), &this->m_UseOriginalData, true, "Use original data in full resolution as final level [default]" );
    cl.AddSwitch( Key( "omit-original-data" ), &this->m_UseOriginalData, false, "Do NOT use original data in full resolution" );
    cl.EndGroup();

    cl.BeginGroup( "Transformation", "Transformation parameters" );
    cl.AddRepeat( Key( "dofs" ), this->NumberDOFs, "Add number of degrees of freedom [can be repeated]" );
    cl.AddRepeat( Key( "dofs-final" ), this->NumberDOFsFinal, "Add number of degrees of freedom for final level only [can be repeated]" );
    
    cl.AddSwitch( Key( 'n', "no-switch" ), &this->m_NoSwitch, 1, "Do not auto-switch reference and floating image for improved computational performance" );
    cl.AddSwitch( Key( 'i', "initxlate" ), &InitXlate, true, "Initialized transformation by translating floating image FOV center onto reference image FOV center" );

    cl.AddOption( Key( "initial" ), &InitialStudylist, "Initialize transformation from given path" );
    cl.AddOption( Key( "initial-inverse" ), &InitialStudylist, "Initialize transformation with inverse from given path", &this->m_InitialXformIsInverse );
    cl.EndGroup();

    cl.BeginGroup( "Image data", "Image data" );
    CommandLine::EnumGroup<int>::SmartPtr
      metricGroup = cl.AddEnum( "registration-metric", &this->m_Metric, "Registration metric for motion estimation by image-to-image registration." );
    metricGroup->AddSwitch( Key( "nmi" ), 0, "Normalized Mutual Information metric" );
    metricGroup->AddSwitch( Key( "mi" ), 1, "Standard Mutual Information metric" );
    metricGroup->AddSwitch( Key( "cr" ), 2, "Correlation Ratio metric" );
    metricGroup->AddSwitch( Key( "msd" ), 4, "Mean Squared Difference metric" );
    metricGroup->AddSwitch( Key( "ncc" ), 5, "Normalized Cross Correlation metric" );

    this->m_PreprocessorRef.AttachToCommandLine( cl );
    this->m_PreprocessorFlt.AttachToCommandLine( cl );

    cl.BeginGroup( "Output", "Output parameters" );
    cl.AddOption( Key( 'o', "outlist" ), &this->Studylist, "Output path for final transformation" );
    cl.AddOption( Key( "out-matrix" ), &this->OutMatrixName, "Output path for final transformation in matrix format" );
    cl.AddOption( Key( "out-parameters" ), &this->OutParametersName, "Output path for final transformation in plain parameter list format" );
    cl.AddOption( Key( 'p', "protocol" ), &this->Protocol, "Optimization protocol output file name" );
    cl.AddOption( Key( 't', "time" ), &this->Time, "Computation time statistics output file name" );
    cl.EndGroup();
    
    cl.Parse();

    clArg1 = cl.GetNext();
    clArg2 = cl.GetNextOptional();
    }
  catch ( CommandLine::Exception ex )
    {
    StdErr << ex << "\n";
    exit( 1 );
    }

  if ( (OptimizerStepFactor <= 0) || (OptimizerStepFactor >= 1) ) 
    {
    StdErr << "ERROR: step factor value " << OptimizerStepFactor << " is invalid. Must be in range (0..1)\n";
    exit( 1 );
    }

  if ( clArg2 ) 
    {
    AffineXform::SmartPtr initialXform( new AffineXform() );
    this->SetInitialXform( initialXform );
    
    Study1 = const_cast<char*>( clArg1 );
    Study2 = const_cast<char*>( clArg2 );
    } 
  else
    {
    inStudylist = clArg1;

    if ( InitialStudylist ) 
      {
      StdErr << "Transformation will be overriden by '--initial' list.\n";
      }
    
    if ( Verbose )
      StdErr << "Reading input studylist " << inStudylist << ".\n";
    
    ClassStream typedStream( MountPoints::Translate(inStudylist), "registration", ClassStream::READ );
    if ( ! typedStream.IsValid() ) 
      {
      StdErr << "Could not open studylist archive " << inStudylist << ".\n";
      exit( 1 );
      }

    typedStream.Seek ( "registration" );
    Study1 = typedStream.ReadString( "reference_study" );
    Study2 = typedStream.ReadString( "floating_study" );
    if ( Study2 )
      {
      AffineXform::SmartPtr affineXform;
      typedStream >> affineXform;
      this->SetInitialXform( affineXform->GetInverse() );
      }
    else
      {
      Study2 = typedStream.ReadString( "model_study" );
      AffineXform::SmartPtr affineXform;
      typedStream >> affineXform;
      this->SetInitialXform( affineXform );
      }

    typedStream.Close();
  }

  if ( InitialStudylist ) 
    {
    Xform::SmartPtr xform( XformIO::Read( InitialStudylist, Verbose ) );
    if ( ! xform ) 
      {
      StdErr << "ERROR: could not read transformation from " << InitialStudylist << "\n";
      exit( 1 );
      }
    
    AffineXform::SmartPtr affine( AffineXform::SmartPtr::DynamicCastFrom( xform ) );
    if ( ! affine )
      {
      StdErr << "ERROR: transformation " << InitialStudylist << " is not affine.\n";
      exit( 1 );
      }
    this->SetInitialXform( affine->GetInverse() );
    }
  
  if ( !Study1 )
    {
    StdErr << "ERROR: reference image path resolved to NULL.\n";
    exit( 1 );
    }
  
  if ( !Study2 )
    {
    StdErr << "ERROR: floating image path resolved to NULL.\n";
    exit( 1 );
    }
  
  UniformVolume::SmartPtr volume( VolumeIO::ReadOriented( Study1, Verbose ) );
  if ( !volume ) throw ConstructorFailed();
  this->SetVolume_1( UniformVolume::SmartPtr( this->m_PreprocessorRef.GetProcessedImage( volume ) ) );

  volume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( Study2, Verbose ) );
  if ( !volume ) throw ConstructorFailed();
  this->SetVolume_2(  UniformVolume::SmartPtr( this->m_PreprocessorFlt.GetProcessedImage( volume ) ) );

  if ( InitXlate ) 
    {
    if ( inStudylist || InitialStudylist ) 
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
  
  if ( Protocol ) 
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
AffineRegistrationCommandLine::OutputResultMatrix( const char* matrixName ) const
{
  Types::Coordinate matrix[4][4];
  this->GetTransformation()->GetMatrix( matrix );

  FILE* mfile = fopen( matrixName, "w" );
  if ( mfile )
    {
    for ( int i = 0; i < 4; ++i )
      {
      fprintf( mfile, "%e\t%e\t%e\t%e\n", matrix[0][i], matrix[1][i], matrix[2][i], matrix[3][i] );
      }
    fclose( mfile );
    }
}

void
AffineRegistrationCommandLine::OutputResultParameters
( const char* paramsName, const CoordinateVector& v ) const
{
  FILE* pfile = fopen( paramsName, "w" );
  if ( pfile )
    {
    for ( unsigned int idx=0; idx < v.Dim; ++idx )
      fprintf( pfile, "#%d: %f\n", idx, v.Elements[idx] );
    fclose( pfile );
    }
}

void
AffineRegistrationCommandLine::OutputResultList( const char* studyList ) const
{
  ClassStream classStream( studyList, "studylist", ClassStream::WRITE );
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
    
  classStream.Open( studyList, "registration", ClassStream::WRITE );
    
  classStream.Begin( "registration" );
  classStream.WriteString( "reference_study", CompressedStream::GetBaseName( Study1 ) );
  classStream.WriteString( "floating_study", CompressedStream::GetBaseName( Study2 ) );
    
  classStream << (*this->GetTransformation()->GetInverse());
    
  classStream.End();
  classStream.Close();
    
  classStream.Open( studyList, "settings", ClassStream::WRITE );
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
    
  classStream.Open( studyList, "statistics", ClassStream::WRITE );
  classStream.WriteDouble( "time", this->GetTotalElapsedTime() );
  classStream.WriteDouble( "walltime", this->GetTotalElapsedWalltime() );
#ifdef CMTK_USE_THREADS
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
AffineRegistrationCommandLine::OutputResult ( const CoordinateVector* v ) const
{
  if ( Verbose ) 
    {
    fprintf( stderr, "\rResulting transformation parameters: \n" );
    for ( unsigned int idx=0; idx<v->Dim; ++idx )
      fprintf( stderr, "#%d: %f\n", idx, v->Elements[idx] );
    }
  
  if ( this->OutMatrixName )
    {
    this->OutputResultMatrix( this->OutMatrixName );
    }

  if ( this->OutParametersName )
    {
    this->OutputResultParameters( this->OutParametersName, *v );
    }

  if ( this->Studylist ) 
    {
    this->OutputResultList( this->Studylist );
    }
}

void
AffineRegistrationCommandLine::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f,
  const int index, const int total )
{
  if ( Verbose )
    fprintf( stderr, "\rEntering resolution level %d out of %d...\n", index, total );
  this->Superclass::EnterResolution( v, f, index, total );
}

CallbackResult
AffineRegistrationCommandLine::Register ()
{
  const double baselineTime = cmtk::Timers::GetTimeProcess();
  CallbackResult Result = Superclass::Register();
  const int elapsed = static_cast<int>( cmtk::Timers::GetTimeProcess() - baselineTime );

  if ( Time ) 
    {
    FILE *tfp = fopen( Time, "w" );
    
    if ( tfp ) 
      {
      fprintf( tfp, "%d\n", elapsed );
      fclose( tfp );
      } 
    else 
      {
      std::cerr << "Could not open time file " << Time << "\n";
      }
    }
  return Result;
}

} // namespace cmtk

