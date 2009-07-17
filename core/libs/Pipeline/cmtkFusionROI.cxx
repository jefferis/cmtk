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

#include <cmtkFusionROI.h>

#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

FusionROI::FusionROI() 
{ 
  TopImageIndex = 0;
  Mode = 0;

  VerticalWidth = 0.1;
  VerticalOffset = 0;
  HorizontalHeight = 0.1;
  HorizontalOffset = 0;
  CircleX = 0.5;
  CircleY = 0.5;
  CircleDelta = 0.1;
  CircleOffset = 0;
  BoxX = 0.5;
  BoxY = 0.5;
  BoxDeltaX = 0.1;
  BoxDeltaY = 0.1;
  CheckerWidth = 0.1;
  CheckerHeight = 0.1;
  CheckerOffsetX = 0;
  CheckerOffsetY = 0;
};

void FusionROI::Execute()
{
  if ( !Input[0] || !Input[1] ) return;

  ImageRGB *output = this->GetOutput();
  output->CopyStructure( Input[0] );
  output->SetAlphaChannel( IMAGE_RGB );

  unsigned int nx = Input[0]->GetDims()[0];
  unsigned int ny = Input[0]->GetDims()[1];

   if ( ( nx != Input[1]->GetDims()[0] ) || ( ny != Input[1]->GetDims()[1] ) ) 
     {
     return;
     }
   
   RGB* i1 = (RGB*) Input[  TopImageIndex]->GetDataPtr();
   RGB* i2 = (RGB*) Input[1-TopImageIndex]->GetDataPtr();

   RGB* outData = (RGB*) output->GetDataPtr( true /* forceAlloc */ );
   enum mode {
     VERTICAL = 0, 
     HORIZONTAL = 1,
     CIRCLE = 2,
     BOX = 3,
     CHECKERBOARD = 4
   };
   
   switch ( Mode ) 
     {
     case VERTICAL:
       VerticalSlice( outData, i1, i2, nx, ny ); 
       break;
     case HORIZONTAL:
       HorizontalSlice( outData, i1, i2, nx, ny ); 
       break;
     case CIRCLE:
       Circle( outData, i1, i2, nx, ny );
       break;
     case BOX:
       Box( outData, i1, i2, nx, ny );
       break;
     case CHECKERBOARD:
       Checkerboard( outData, i1, i2, nx, ny ); 
       
       break;
     }
   this->UpdateExecuteTime();
}

void
FusionROI::Circles 
( RGB *target, const RGB *source1, const RGB *source2, const int nx, const int ny )
{ 
  const int cX = (int) (nx * CircleX);
  const int cY = (int) (ny * CircleY);
  
  // force a reasonable minimum spacing
  const double rDelta = std::max<double>( CircleDelta * std::max( nx, ny ), 2.0);

  int toggle;
  for ( int y=0; y<ny; ++y, target+=nx, source1+=nx, source2+=nx ) 
    {
    double radius = sqrt( (double) (MathUtil::Square( cX ) + MathUtil::Square( y-cY )) );
    int radiusIdx = (int) ( radius/rDelta + CircleOffset );
    
    toggle = (radiusIdx & 1);
    for ( int x=0; x<nx; toggle ^= 1 ) 
      {
      int delta;
      if ( x < cX ) 
	{
	// left of center: circles get smaller
	const double minRadius = rDelta * ( radiusIdx - CircleOffset );
	
	if ( minRadius > abs( y-cY ) )
	  // intersection with inner circle exists
	  delta = (int) ( ( cX-x ) - sqrt( (double)( MathUtil::Square( minRadius ) - MathUtil::Square( y-cY ) ) ) );
	else
	  // all smaller circles are passed
	  delta = 2 * ( cX-x );
	} 
      else
	{
	// right of center: circles get bigger
	const double maxRadius = rDelta * ( 1 + radiusIdx - CircleOffset );
	delta = (int) sqrt( (double)( MathUtil::Square( maxRadius ) - MathUtil::Square( y-cY ) ) ) + ( cX-x );
	}
      
      delta = std::max( 1, std::min( 1+delta, nx-x ) );
      
      // copy correct field's content
      if ( toggle )
	memcpy( target+x, source2+x, delta*sizeof( RGB ) );
      else
	memcpy( target+x, source1+x, delta*sizeof( RGB ) );
      
      x += delta;
      radius = sqrt( (double)( MathUtil::Square( x-cX ) + MathUtil::Square( y-cY ) ) );
      radiusIdx = (int) ( radius/rDelta + CircleOffset );
      }
    }
}

void
FusionROI::Circle
( RGB *target, const RGB *source1, const RGB *source2, const int nx, const int ny )
{
  const int cX = (int) (nx * CircleX);
  const int cY = (int) (ny * CircleY);

  // force a sensible minimum radius
  const int Radius = std::max( MathUtil::Round( CircleDelta * std::max( nx, ny ) ), 2 );
  const int squareRadius = MathUtil::Square( Radius );

  const int startY = std::max( 0, cY - Radius );
  const int skipStart = nx * startY;
  if ( skipStart > 0 )
    memcpy( target, source2, skipStart * sizeof( RGB ) );
  
  const int endY = std::min( ny, cY + Radius );
  const int skipEnd = nx * ( ny - endY );
  if ( skipEnd > 0 )
    memcpy( target + (cY + Radius)*nx, source2 + (cY + Radius)*nx, skipEnd * sizeof( RGB) );
  
  target  += skipStart; 
  source1 += skipStart;
  source2 += skipStart;

  for ( int y=startY; y<endY; ++y, target+=nx, source1+=nx, source2+=nx ) 
    {
    const int dX = (int) sqrt( (double)(squareRadius - MathUtil::Square( cY-y )) );
    
    const int startLineSkip = std::max( 0, cX-dX );
    const int endLineSkip = std::max( 0, nx - (cX+dX) );
    const int endX = nx - endLineSkip;

    if ( startLineSkip > 0 )
      memcpy( target, source2, startLineSkip * sizeof( RGB ) );

    if ( endX > startLineSkip )
      memcpy( target+startLineSkip, source1+startLineSkip, (endX-startLineSkip) * sizeof( RGB ) );

    if ( endLineSkip > 0 )
      memcpy( target+endX, source2+endX, endLineSkip * sizeof( RGB ) );
    }
};

void
FusionROI::Box 
( RGB *target,const RGB *source1, const RGB *source2, const int nx, const int ny )
{
  const int cX = (int) (nx * BoxX);
  const int cY = (int) (ny * BoxY);

  // force a sensible minimum size
  const int sX = std::max( MathUtil::Round( BoxDeltaX * nx), 2 );
  const int sY = std::max( MathUtil::Round( BoxDeltaY * ny), 2 );

  int startY = std::max( 0, cY - sY );
  int skipStart = nx * startY;
  if ( skipStart > 0 )
    memcpy( target, source2, skipStart * sizeof( RGB ) );
  
  int endY = std::min( ny, cY + sY );
  int skipEnd = nx * ( ny - endY );
  if ( skipEnd > 0 )
    memcpy( target + (cY + sY)*nx, source2 + (cY + sY)*nx, skipEnd * sizeof( RGB) );
  
  target  += skipStart; 
  source1 += skipStart;
  source2 += skipStart;

  const int startLineSkip = std::max( 0, cX-sX );
  const int endLineSkip = std::max( 0, nx - (cX+sX) );
  const int endX = nx - endLineSkip;

  for ( int y=startY; y<endY; ++y, target+=nx, source1+=nx, source2+=nx ) 
    {
    if ( startLineSkip > 0 )
      memcpy( target, source2, startLineSkip * sizeof( RGB ) );
    
    if ( endX > startLineSkip )
      memcpy( target+startLineSkip, source1+startLineSkip, (endX-startLineSkip) * sizeof( RGB ) );
    
    if ( endLineSkip > 0 )
      memcpy( target+endX, source2+endX, endLineSkip * sizeof( RGB ) );
    }
}

void
FusionROI::VerticalSlice
( RGB* target, const RGB *source1, const RGB *source2, const int nx, const int ny, const int xwidth, const int xoffset )
{
  int width;
  if ( xwidth == -1 ) 
    width = (int) (nx * VerticalWidth);
  else
    width = xwidth;

  int offset;
  if ( xoffset == -1 )
    offset = (int) (nx * VerticalWidth * VerticalOffset);
  else
    offset = xoffset;

  const int realWidth = std::min( nx, std::max( 1, width ) );

  for ( int y=0; y<ny; ++y, target+=nx, source1+=nx, source2+=nx ) 
    {
    int toggle = 0;
    
    memcpy( target, source1, offset * sizeof( RGB ) );

    int idx = offset;
    while ( idx < nx ) 
      {
      toggle ^= 1;
      if ( toggle )
	memcpy( target+idx, source2+idx, sizeof( RGB ) * std::min( realWidth, nx-idx ) );
      else
	memcpy( target+idx, source1+idx, sizeof( RGB ) * std::min( realWidth, nx-idx ) );
      idx += realWidth;
      }
    }
}

void 
FusionROI::HorizontalSlice 
( RGB* target, const RGB *source1, const RGB *source2, const int nx, const int ny )
{
  const int height = (int) (ny * HorizontalHeight);
  const int offset = (int) (ny * HorizontalHeight * HorizontalOffset);

  const int realHeight = std::min( ny, std::max( 1, height ) );

  const int rowsize = nx * sizeof( RGB );
  memcpy( target, source1, rowsize*offset );

  const int skipBlock = offset * nx;
  target  += skipBlock;
  source1 += skipBlock;
  source2 += skipBlock;
  
  const int blocksize = realHeight * nx;
  int toggle = 0;

  for ( int idx = offset; idx<ny; idx+=realHeight, target+=blocksize, source1+=blocksize, source2+=blocksize ) 
    {
    toggle ^= 1;
    if ( toggle )
      memcpy( target, source2, rowsize * std::min( realHeight, ( ny-idx ) ) );
    else
      memcpy( target, source1, rowsize * std::min( realHeight, (ny-idx) ) );
    }
}

void
FusionROI::Checkerboard
( RGB* target, const RGB *source1, const RGB *source2, const int nx, const int ny )
{
  const int width = (int) (nx * CheckerWidth);
  const int offsetX = (int) (width * CheckerOffsetX);
  
  const int height = (int) (ny * CheckerHeight);
  const int offsetY = (int) (height * CheckerOffsetY);

  const int realHeight = std::min( ny, std::max(1, height));
  
  VerticalSlice( target, source1, source2, nx, offsetY, width, offsetX );
  
  int skipBlock = offsetY * nx;
  target  += skipBlock;
  source1 += skipBlock;
  source2 += skipBlock;

  const int blocksize = realHeight * nx;
  int toggle = 0;

  for ( int idx = offsetY; idx<ny; idx+=realHeight, target+=blocksize, source1+=blocksize, source2+=blocksize ) 
    {
    toggle ^= 1;
    
    if ( toggle )
      VerticalSlice( target, source2, source1, nx, std::min( realHeight, (ny-idx)), width, offsetX );
    else
      VerticalSlice( target, source1, source2, nx, std::min( realHeight, (ny-idx)), width, offsetX );
  }
}

} // namespace cmtk
