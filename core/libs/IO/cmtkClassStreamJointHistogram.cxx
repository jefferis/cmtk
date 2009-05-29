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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkClassStream.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
ClassStream::operator >> ( JointHistogram<float>*& histogram )
{
  histogram = NULL;
  if ( this->Seek( "JointHistogram" ) != TYPEDSTREAM_OK ) return *this;

  int numBinsX = this->ReadInt( "NumBinsX" );
  Types::DataItem range[2];
  this->ReadItemArray( "RangeX", range, 2 );

  int numBinsY = this->ReadInt( "NumBinsY" );
  histogram = new JointHistogram<float>( numBinsX, numBinsY );
  histogram->SetRangeX( range[0], range[1] );

  this->ReadItemArray( "RangeY", range, 2 );
  histogram->SetRangeY( range[0], range[1] );

  float *bins = Memory::AllocateArray<float>( numBinsX * numBinsY );
  this->ReadFloatArray( "Bins", bins, numBinsX * numBinsY );
  histogram->SetBins( bins );
  delete[] bins;

  this->End();

  return *this;
}

ClassStream& 
ClassStream::operator >> ( JointHistogram<int>*& histogram )
{
  histogram = NULL;
  if ( this->Seek( "JointHistogram" ) != TYPEDSTREAM_OK ) return *this;

  int numBinsX = this->ReadInt( "NumBinsX" );
  Types::DataItem range[2];
  this->ReadItemArray( "RangeX", range, 2 );

  int numBinsY = this->ReadInt( "NumBinsY" );
  histogram = new JointHistogram<int>( numBinsX, numBinsY );
  histogram->SetRangeX( range[0], range[1] );

  this->ReadItemArray( "RangeY", range, 2 );
  histogram->SetRangeY( range[0], range[1] );

  int *bins = Memory::AllocateArray<int>( numBinsX * numBinsY );
  this->ReadIntArray( "Bins", bins, numBinsX * numBinsY );
  histogram->SetBins( bins );
  delete[] bins;

  this->End();

  return *this;
}

ClassStream& 
ClassStream::operator << ( const JointHistogram<float> *histogram )
{
  int numBinsX = histogram->GetNumBinsX();
  int numBinsY = histogram->GetNumBinsY();

  this->Begin( "JointHistogram" );
  this->WriteInt( "NumBinsX", numBinsX);

  Types::DataItem range[2];
  histogram->GetRangeX( range[0], range[1] );
  this->WriteItemArray( "RangeX", range, 2 );

  this->WriteInt( "NumBinsY", numBinsY );

  histogram->GetRangeY( range[0], range[1] );
  this->WriteItemArray( "RangeY", range, 2 );

  float *bins = histogram->GetBins();
  this->WriteFloatArray( "Bins", bins, numBinsX * numBinsY );
  delete[] bins;

  this->End();
  return *this;
}

ClassStream& 
ClassStream::operator << ( const JointHistogram<int> *histogram )
{
  int numBinsX = histogram->GetNumBinsX();
  int numBinsY = histogram->GetNumBinsY();

  this->Begin( "JointHistogram" );
  this->WriteInt( "NumBinsX", numBinsX);

  Types::DataItem range[2];
  histogram->GetRangeX( range[0], range[1] );
  this->WriteItemArray( "RangeX", range, 2 );

  this->WriteInt( "NumBinsY", numBinsY );

  histogram->GetRangeY( range[0], range[1] );
  this->WriteItemArray( "RangeY", range, 2 );

  int *bins = histogram->GetBins();
  this->WriteIntArray( "Bins", bins, numBinsX * numBinsY );
  delete[] bins;

  this->End();
  return *this;
}

} // namespace cmtk
