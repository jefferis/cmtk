#!perl

PREFIX: while ( <STDIN> ) 
  {
    if ( /^\*\// ) 
      {
	last PREFIX;
      }
  }

print "/*
//  Copyright (c) Torsten Rohlfing 1997-2007. All Rights Reserved.
//
//  \$Id\$
//
//  THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE.
//
//  The contents of this file may not be disclosed to third parties,
//  copied or duplicated in any form, in whole or in part, without
//  the prior written permission of the copyright holder.
//
//  \$Revision\$
//
//  \$LastChangedDate\$
//
//  \$LastChangedBy\$
*/
";

BODY: while ( <STDIN> ) 
  {
    if ( /^\/\/ \$Log/ ) 
      {
	last BODY;
      }
    print $_;
  }

print "//  \$Date\$
//
//  \$Author\$
*/\n";
