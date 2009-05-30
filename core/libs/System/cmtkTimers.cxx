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

#include <cmtkTimers.h>

#include <string.h>
#include <stdio.h>
#include <errno.h>

#ifdef HAVE_SYS_IOCTL_H
#  include <sys/ioctl.h>
#endif

#ifdef HAVE_FCNTL_H
#  include <fcntl.h>
#endif

#ifdef HAVE_SYS_PROCFS_H
#  include <sys/procfs.h>
#endif

#ifdef CMTK_USE_THREADS
#  include <pthread.h>
#endif

double
cmtk::Timers::GetTimeProcess()
{
#ifdef _MSC_VER
  return (double) clock();
#else
  struct tms t; 
  if ( times(&t) )
#ifdef _SC_CLK_TCK
    return (double)( t.tms_utime + t.tms_cutime + t.tms_stime + t.tms_cstime )/sysconf(_SC_CLK_TCK);
#else
    return (double)( t.tms_utime + t.tms_cutime + t.tms_stime + t.tms_cstime )/CLK_TCK;
#endif
  else
    return 0;
#endif // #ifdef _MSC_VER
}

double
cmtk::Timers::GetWalltime()
{
  return time( NULL );
}

double
cmtk::Timers::GetTimeThread()
{
#ifndef _MSC_VER
#ifdef HAVE_SYS_PROCFS_H

  char buffer[80];

  snprintf( buffer, sizeof( buffer ), "/proc/%ld/usage", (long)getpid() );
  FILE *fp = fopen( buffer, "r" );
  if ( fp ) {
#ifdef PRUSAGE_T_IN_SYS_PROCFS_H
    prusage_t usage;
    if ( fread( &usage, sizeof( usage ), 1, fp ) ) {
      fclose( fp );
      //      return ((double) usage.pr_rtime.tv_sec) + 
      //      	((double) usage.pr_rtime.tv_nsec) / NANOSEC;
      return ((double) usage.pr_utime.tv_sec) + ((double) usage.pr_utime.tv_nsec) / NANOSEC +
	((double) usage.pr_stime.tv_sec) + ((double) usage.pr_stime.tv_nsec) / NANOSEC;
    }
#endif
    fclose( fp );
  }
  return 0;

#else // #ifdef HAVE_SYS_PROCFS_H
  struct tms t; 
  if ( times(&t) )
#ifdef _SC_CLK_TCK
    return (double)(t.tms_utime + t.tms_stime)/sysconf(_SC_CLK_TCK);
#else
    return (double)(t.tms_utime + t.tms_stime)/CLK_TCK;
#endif
  else
    return 0;
#endif // #ifdef HAVE_SYS_PROCFS_H
#else // ifndef _MSC_VER
  return (double) clock();
#endif // ifndef _MSC_VER
}
