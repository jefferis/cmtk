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

#include <cmtkRenderer.h>

#include <cmtkRendererCollection.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

Renderer::Renderer()
{
  Input = NULL;
  Active = true;
  RenderPending = false;
}

Renderer::~Renderer()
{
  if ( Input != NULL ) Input->Delete();
  if ( RendererCollectionInstance != NULL ) {
    // this is a hack to avoid recursive destructor calls:
    this->Reference();
    RendererCollectionInstance->RemoveRenderer( this );
  }
}

void
Renderer::RegisterToCollection()
{
  if ( RendererCollectionInstance == NULL ) {
    RendererCollectionInstance = RendererCollection::New();
  }
  RendererCollectionInstance->AddRenderer( this );
}

void
Renderer::SetInput( ImageRGB *const input )
{
  ReplaceObject( Input, input );
};

long
Renderer::Update()
{
  if ( this->IsActive() )
    this->CheckInputForUpdate( Input );
  return this->Superclass::Update();
}

void 
Renderer::Render()
{
  // Is this renderer being updated already, ie. do we have a recursion here?
  // If no: go through it.
  if ( ! RenderPending ) {
    RenderPending = true;

    // Fake modification to make sure we WILL update ourselves.
    this->UpdateModifiedTime();
    this->Update();

    RenderPending = false;
  }
}

void 
Renderer::WritePPM( const char* filename, const int cmtc, const char** cmtv )
{
  this->Update();
  ImageRGB *capture = this->CaptureDisplay();
  if ( capture == NULL ) return;

  const byte* imageDataRGB = static_cast<const byte*>( capture->GetDataPtr() );
  if ( imageDataRGB == NULL ) return;

  unsigned int width, height;
  capture->GetDims( width, height );

  FILE *fp;
  if ( ( fp = fopen( filename, "wb" ) ) == NULL ) {
    return;
  }
  
  int Result;
  if ( this->GetEffectiveExportColorMode( capture ) == EXPORTCOLOR_RGB ) 
    {
    fprintf( fp, "P6\n" );
    for ( int i = 0; i < cmtc; ++i ) 
      {
      fprintf( fp, "# %s\n", cmtv[i] );
      }
    
    fprintf( fp, "%d %d\n255\n", width, height );
    
    if ( capture->GetAlphaChannel() == IMAGE_RGB )
      Result = fwrite( imageDataRGB, 3, width * height, fp );
    else
      {
      Result = 0;
      for ( size_t px = 0; px < width * height; ++px )
	{
	Result += fwrite( imageDataRGB+4*px, 3, 1, fp );
	}
      }
    } 
  else 
    {
    fprintf( fp, "P5\n" );
    for ( int i = 0; i < cmtc; ++i ) 
      {
      fprintf( fp, "# %s\n", cmtv[i] );
      }
    
    fprintf( fp, "%d %d\n255\n", width, height );
    
    byte *greyData = Memory::AllocateArray<byte>( width*height );
    for ( unsigned int idx=0; idx < width*height; ++idx, imageDataRGB += 3 ) 
      {
      greyData[idx] = MathUtil::Round( 0.30 * imageDataRGB[0] + 0.59 * imageDataRGB[1] + 0.11 * imageDataRGB[2] );
      }
    Result = fwrite( greyData, 1, width * height, fp );
    delete[] greyData;
  }
  
  if ( Result != static_cast<int>(width * height) ) {
    //    InteractionInstance->SetVError( "Could not write image to file %s.",
    //				       filename );
  }
  fclose( fp );

  this->UpdateExecuteTime();
  capture->Delete();
}

ExportColorMode 
Renderer::GetEffectiveExportColorMode( const ImageRGB* capture ) const
{
  if ( (this->m_ExportColorMode == EXPORTCOLOR_AUTO) && capture->IsGreyscale() )
    return EXPORTCOLOR_GREY;

  return EXPORTCOLOR_RGB;
}

} // namespace cmtk
