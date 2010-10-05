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

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>
#include <System/cmtkMountPoints.h>

#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkJointHistogram.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkUniformVolumeInterpolatorBase.h>
#include <Base/cmtkUniformVolumeInterpolator.h>
#include <Base/cmtkUniformVolumeInterpolatorPartialVolume.h>

#include <Registration/cmtkVoxelMatchingCrossCorrelation.h>

#include <stdio.h>
#include <string.h>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace similarity
{
#endif
bool Verbose = false;

bool SwapVolumes = false;
bool Absolute = false;

bool Padding0 = false;
bool Padding1 = false;

cmtk::Types::DataItem PaddingValue0 = 0;
cmtk::Types::DataItem PaddingValue1 = 0;

bool LabelMode = false;

cmtk::UniformVolume::SmartPtr Volume0;
cmtk::UniformVolume::SmartPtr Volume1;

bool OutsideBG = false;

const char *Study0 = NULL;
const char *Study1 = NULL;

const char* MaskFileName = NULL;

unsigned int ForceMaxLabel = 0;

const char* HistogramTextFileName = NULL;

bool 
ParseCommandLine( const int argc, const char* argv[] )
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image similarity" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Compute similarity between two images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] image0 image1" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );
    
    cl.AddSwitch( Key( 's', "swap" ), &SwapVolumes, true, "Swap reference and floating image" );
    cl.AddSwitch( Key( "outside-bg" ), &OutsideBG, true, "Assume voxels outside floating image to be background" );
    
    cl.AddOption( Key( "pad0" ), &PaddingValue0, "Define padding value for reference image", &Padding0 );
    cl.AddOption( Key( "pad1" ), &PaddingValue1, "Define padding value for floating image", &Padding1 );
    
    cl.AddSwitch( Key( 'g', "grey" ), &LabelMode, false, "Image pixels are intensities [default]" );
    cl.AddSwitch( Key( 'l', "labels" ), &LabelMode, true, "Image pixels are labels" );
    
    cl.AddOption( Key( "study0" ), &Study0, "Override reference study in studylist" );
    cl.AddOption( Key( "study1" ), &Study1, "Override floating study in studylist" );
    cl.AddOption( Key( 'm', "mask" ), &MaskFileName, "Mask file for optional region-based analysis" );
    
    cl.AddOption( Key( "force-max" ), &ForceMaxLabel, "Force maximum label value" );
    
    cl.AddOption( Key( "histogram-text-file" ), &HistogramTextFileName, "Output file name for histogram (plain text format)." );
    
    if ( !cl.Parse( argc, argv ) ) return false;
    
    Study0 = cl.GetNext();
    Study1 = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  return true;
}

cmtk::JointHistogram<int>*
AnalyseStudies
( unsigned int& voxelCount, double& sumAbs, double& sumSq, double& sumMult, 
  unsigned int& countVoxelsUnequal, cmtk::VoxelMatchingCrossCorrelation& ccMetric, 
  const cmtk::TypedArray* mask = NULL )
{
  sumAbs = sumSq = sumMult = 0;
  voxelCount = countVoxelsUnequal = 0;

  ccMetric.Reset();

  const cmtk::TypedArray *data0 = Volume0->GetData();
  const cmtk::TypedArray *data1 = Volume1->GetData();

  cmtk::JointHistogram<int> *histogram = NULL;
  if ( LabelMode ) 
    {
    const cmtk::Types::DataItemRange range0 = data0->GetRange();
    const cmtk::Types::DataItemRange range1 = data1->GetRange();

    histogram = new cmtk::JointHistogram<int>( 1+static_cast<size_t>( range0.Width() ), 1+static_cast<size_t>( range1.Width() ) );
    histogram->SetRangeCenteredX( range0 );
    histogram->SetRangeCenteredY( range1 );
    }
  else
    {
    histogram = new cmtk::JointHistogram<int>( 256, 256 );
    
    histogram->SetRangeCenteredX( data0->GetRange() );
    histogram->SetRangeCenteredY( data1->GetRange() );
    }
  
  unsigned int recog = 0, subst = 0, reject = 0;
  unsigned int recogFG = 0, substFG = 0, rejectFG = 0, voxelCountFG = 0;
  unsigned int recogM = 0, substM = 0, rejectM = 0, voxelCountM = 0;
  cmtk::Types::DataItem value0, value1, maskValue;

  unsigned int numberOfVoxels = std::min( data0->GetDataSize(), data1->GetDataSize() );
  
  bool insideMask = true;
  for ( size_t r = 0; r < numberOfVoxels; ++r ) 
    {
    bool dataExist0 = data0->Get( value0, r );
    
    bool dataExist1 = data1->Get( value1, r );
    if ( ! dataExist1 && OutsideBG )
      {
      value1 = 0;
      dataExist1 = true;
      }

    if ( dataExist0 && dataExist1 ) 
      {
      insideMask = !mask || (( mask->Get( maskValue, r )  && maskValue ));
      
      if ( insideMask )
	{
	sumAbs += fabs( value0 - value1 );
	sumSq += cmtk::MathUtil::Square( value0 - value1 );
	sumMult += value0 * value1;
	}
      
      voxelCount++;
      if ( value0 != 0 ) ++voxelCountFG;
      if ( insideMask )
	{
	++voxelCountM;       
	histogram->Increment( histogram->ValueToBinX( value0 ), histogram->ValueToBinY( value1 ) );	
	ccMetric.Increment( value0, value1 );
	}
      
      if ( value0 == value1 ) 
	{
	++recog;
	if ( value0 != 0 ) ++recogFG;
	if ( insideMask ) ++recogM;
	} 
      else 
	{
	++subst;
	if ( value0 != 0 ) ++substFG;
	if ( insideMask ) ++substM;
	countVoxelsUnequal++;
	}
      } 
    else
      {
      if ( dataExist0 ) 
	{
	++reject;
	if ( value0 != 0 ) 
	  {
	  ++voxelCountFG;
	  ++rejectFG;
	  }
	if ( insideMask ) 
	  {
	  ++voxelCountM;
	  ++rejectM;
	  }
	}
      if ( dataExist0 || dataExist1 ) 
	countVoxelsUnequal++;
      }
    }
  
  fprintf( stdout, "\nCL\trecog\tsubst\treject\treliab\nCL-all\t%.4f\t%.4f\t%.4f\t%.4f\n",
	   1.0*recog / numberOfVoxels, 1.0*subst / numberOfVoxels, 1.0*reject / numberOfVoxels, 1.0 * recog / (numberOfVoxels - reject) );
  fprintf( stdout, "CL-fg\t%.4f\t%.4f\t%.4f\t%.4f\n",
	   1.0 * recogFG / voxelCountFG, 1.0 * substFG / voxelCountFG, 1.0 * rejectFG / voxelCountFG, 1.0 * recogFG / (voxelCountFG - rejectFG) );
  
  if ( mask ) 
    {
    fprintf( stdout, "CL-mask\t%.4f\t%.4f\t%.4f\t%.4f\n",
	     1.0 * recogM / voxelCountM, 1.0 * substM / voxelCountM, 1.0 * rejectM / voxelCountM, 1.0 * recogM / (voxelCountM - rejectM) );
    }
  fputs( "\n", stdout );
  
  return histogram;
}

int
main ( const int argc, const char* argv[] ) 
{
  if ( ! ParseCommandLine( argc, argv ) ) return 1;

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( Study0, Verbose ) );
  if ( ! volume ) exit( 1 );
  Volume0 = volume;
  if ( Padding0 ) 
    {
    Volume0->GetData()->SetPaddingValue( PaddingValue0 );
    }
  
  cmtk::TypedArray::SmartPtr mask;
  if ( MaskFileName ) 
    {
    cmtk::UniformVolume::SmartPtr maskVolume( cmtk::VolumeIO::ReadOriented( MaskFileName, Verbose ) );
    if ( maskVolume )
      mask = maskVolume->GetData();
    }
  
  volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( Study1, Verbose ) );
  if ( ! volume ) exit( 1 );
  Volume1 = volume;
  if ( Padding1 ) 
    {
    Volume1->GetData()->SetPaddingValue( PaddingValue1 );
    }
  
  unsigned int voxelCount;
  double sumAbs; 
  double sumSq;
  double sumMult;
  unsigned int countVoxelsUnequal;
  
  cmtk::VoxelMatchingCrossCorrelation ccMetric;
  
  cmtk::JointHistogram<int>::SmartPtr histogram( AnalyseStudies( voxelCount, sumAbs, sumSq, sumMult, countVoxelsUnequal, ccMetric, mask ) );
  
  double hX, hY;
  histogram->GetMarginalEntropies( hX, hY );
  const double hXY = histogram->GetJointEntropy();
  
  fprintf( stdout, "STAT\tN\tHX\tHY\nSTATval\t%d\t%.5f\t%.5f\n\n", voxelCount, hX, hY );
  
  fprintf( stdout, "SIM\tDIFF\tMSD\tMAD\tNCC\tHXY\tMI\tNMI\nSIMval\t%d\t%.1f\t%.3f\t%.4f\t%.5f\t%.5f\t%.5f\n",
	   countVoxelsUnequal, sumSq / voxelCount, sumAbs / voxelCount, ccMetric.Get(), hXY, hX + hY - hXY, ( hX + hY ) / hXY );
  
  if ( HistogramTextFileName )
    {
    FILE *file = fopen( HistogramTextFileName, "w" );
    if ( file ) 
      {
      unsigned int binsX = histogram->GetNumBinsX();
      unsigned int binsY = histogram->GetNumBinsY();
      
      for ( unsigned int i = 0; i < binsX; ++i ) 
	{
	for ( unsigned int j = 0; j < binsY; ++j ) 
	  {
	  fprintf( file, "%d\t", histogram->GetBin( i, j ) );
	  }
	fputs( "\n", file );
	}
      fclose( file );
      }
    }
  
  if ( LabelMode ) 
    {
    unsigned int numberLabelsRef = histogram->GetNumBinsX();
    unsigned int numberLabelsFlt = histogram->GetNumBinsY();
    
    // SI : Similarity Index as defined by Dawant et al., TMI 18(10), 1999
    fputs( "\n\t\tTotal\tCorrect\t%\tFalseN\t%\tFalseP\t%\tSI\n", stdout );
    
    unsigned int sumTotal = 0, sumCorrect = 0, countLabels = 0;
    double sumSI = 0;
    double sumSIweighted = 0;
    double avgRecog = 0;
    
    for ( unsigned int i = 0; i < numberLabelsRef; ++i ) 
      {
      unsigned int totalRef = histogram->ProjectToX( i );
      
      if ( totalRef ) 
	{
	const unsigned int correct = ( i < numberLabelsFlt ) ? histogram->GetBin( i, i ) : 0;
	avgRecog += static_cast<double>( (1.0 * correct) / totalRef );
	
	unsigned int missed = 0, wrong = 0;
	for ( unsigned int j = 0; j < numberLabelsFlt; ++j ) 
	  {
	  if ( i != j ) 
	    {
	    if ( (i < numberLabelsFlt) && (j < numberLabelsRef) )
	      {
	      missed += histogram->GetBin( i, j );
	      wrong += histogram->GetBin( j, i );
	      }
	    }
	  }

	const unsigned int totalFlt = ( i < numberLabelsFlt ) ? histogram->ProjectToY( i ) : 0;
	const double SI = static_cast<double>( 2.0 * correct / (totalRef + totalFlt ) );
	fprintf( stdout, "\nLabel #%d:\t%d\t%d\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%.5lf", 
		 i, totalRef, correct, 100.0 * correct / totalRef, missed, 100.0 * missed / totalRef, wrong, 100.0 * wrong / totalRef, SI );
	
	if ( i )
	  { // assume background is label #0 and exclude this
	  ++countLabels;
	  sumTotal += totalRef;
	  sumCorrect += correct;
	  sumSI += SI;
	  sumSIweighted += SI * totalRef;
	  }
	} 
      else 
	{
	// this label does not exist in the reference image
	if ( i <= ForceMaxLabel ) 
	  fprintf( stdout, "\nLabel #%d\t<n/a>", i );
	}
      }
    avgRecog /= (countLabels+1);
    
    fputs( "\n\n\tsumTot\tsumCorr\t%corr\tavgSI\tWmeanSI\n", stdout );
    fprintf( stdout, "Total:\t%d\t%d\t%.2lf\t%.4lf\t%.4lf\n\n", sumTotal, sumCorrect, 100.0 * sumCorrect / sumTotal, sumSI / countLabels, sumSIweighted / sumTotal );
    fprintf( stdout, "AvgRecog:\t%.6f\n", avgRecog );
    }
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace similarity
} // namespace apps
} // namespace cmtk
#endif
