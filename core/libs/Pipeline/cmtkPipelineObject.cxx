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

#include <cmtkPipelineObject.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

PipelineObject::PipelineObject()
{ 
  Owner = NULL;
  // We may use 0 here as the actual time starts at "1". So newly created
  // objects are never marked up-to-date.
  ExecuteTime = 0;
  ExecutePending = 0;
}

int PipelineObject::Register( PipelineObject *const owner )
{
  if ( owner ) 
    {
    Owner = owner;
    if ( PipelineDebugMode ) 
      {
      fprintf( stderr, "Registering %s as owner of %s.\n", owner->GetClassName(), this->GetClassName() );
      }
    }
  this->Object::Reference();
  return this->GetReferenceCount();
}


void PipelineObject::Unregister( PipelineObject *const owner )
{
  // Does the primary owner unregister? Then set primary owner to "none".
  if ( Owner == owner ) Owner = NULL;

  if ( PipelineDebugMode ) 
    {
    fprintf( stderr, "Unregistering %s as owner of %s.\n", owner->GetClassName(), this->GetClassName() );
    }
  
  this->Delete();
}

long PipelineObject::Update()
{
  if ( PipelineDebugMode ) 
    {
    fprintf( stderr, "Entering %s::Update\n", this->GetClassName() );
    }
  this->CheckInputForUpdate( Owner );
  return this->ExecuteIfNecessary();
}

int PipelineObject::CheckInputForUpdate( PipelineObject *const object )
{
  if ( object ) 
    {
    const long ObjectTime = object->Update();
    if ( ObjectTime > ExecuteTime ) 
      {
      ExecutePending = 1;
      return 1;
      }
    }
  return 0;
}

long PipelineObject::ExecuteIfNecessary()
{
  if ( PipelineDebugMode ) 
    {
    fprintf( stderr, "Entering %s::ExecuteIfNecessary\n", this->GetClassName() );
    }
  if ( (this->GetModifiedTime() > ExecuteTime) || ExecutePending ) 
    {
    if ( PipelineDebugMode ) 
      {
      fprintf( stderr, "Calling %s::Execute (MTime %ld > ETime %ld, Pending %d)\n", this->GetClassName(), this->GetModifiedTime(), ExecuteTime, ExecutePending );
      }
    this->Execute();
    this->UpdateExecuteTime();
  } 
  return ExecuteTime;
}

} // namespace cmtk
