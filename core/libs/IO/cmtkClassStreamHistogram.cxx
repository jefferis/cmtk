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

#include <cmtkClassStream.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& ClassStream::operator >> ( Histogram<float>*& histogram )
{
  histogram = NULL;
  if ( this->Seek( "Histogram" ) != TYPEDSTREAM_OK ) return *this;

  int numBins = this->ReadInt( "NumBins" );
  Types::DataItem range[2];
  this->ReadItemArray( "Range", range );

  histogram = new Histogram<float>( numBins );
  histogram->SetRange( range[0], range[1] );

  float *bins = Memory::AllocateArray<float>(  numBins  );
  this->ReadFloatArray( "Bins", bins, numBins );
  histogram->SetBins( bins );
  delete[] bins;

  this->End();

  return *this;
}

ClassStream& ClassStream::operator >> ( Histogram<int>*& histogram )
{
  histogram = NULL;
  if ( this->Seek( "Histogram" ) != TYPEDSTREAM_OK ) return *this;

  int numBins = this->ReadInt( "NumBins" );
  Types::DataItem range[2];
  this->ReadItemArray( "Range", range );

  histogram = new Histogram<int>( numBins );
  histogram->SetRange( range[0], range[1] );

  int *bins = Memory::AllocateArray<int>(  numBins  );
  this->ReadIntArray( "Bins", bins, numBins );
  histogram->SetBins( bins );
  delete[] bins;

  this->End();

  return *this;
}

ClassStream& 
ClassStream::operator << ( const Histogram<float> *histogram )
{
  int numBins = histogram->GetNumBins();

  this->Begin( "Histogram" );
  this->WriteInt( "NumBins", numBins );

  Types::DataItem range[2];
  histogram->GetRange( range[0], range[1] );
  this->WriteItemArray( "Range", range, 2 );

  float *bins = histogram->GetBins();
  this->WriteFloatArray( "Bins", bins, numBins );
  delete[] bins;
  
  this->End();
  return *this;
}

ClassStream& 
ClassStream::operator << ( const Histogram<int> *histogram )
{
  int numBins = histogram->GetNumBins();

  this->Begin( "Histogram" );
  this->WriteInt( "NumBins", numBins );

  Types::DataItem range[2];
  histogram->GetRange( range[0], range[1] );
  this->WriteItemArray( "Range", range, 2 );

  int *bins = histogram->GetBins();
  this->WriteIntArray( "Bins", bins, numBins );
  delete[] bins;
  
  this->End();
  return *this;
}

} // namespace cmtk
