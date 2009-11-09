##
##  Copyright 1997-2009 Torsten Rohlfing
##  Copyright 2004-2009 SRI International
##
##  This file is part of the Computational Morphometry Toolkit.
##
##  http://www.nitrc.org/projects/cmtk/
##
##  The Computational Morphometry Toolkit is free software: you can
##  redistribute it and/or modify it under the terms of the GNU General Public
##  License as published by the Free Software Foundation, either version 3 of
##  the License, or (at your option) any later version.
##
##  The Computational Morphometry Toolkit is distributed in the hope that it
##  will be useful, but WITHOUT ANY WARRANTY; without even the implied
##  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with the Computational Morphometry Toolkit.  If not, see
##  <http://www.gnu.org/licenses/>.
##
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

import sys
import socket
import gzip
import difflib
import os
import io

BUILDNAME=sys.argv[1]
BINDIR=sys.argv[2]
DATADIR=sys.argv[3]
RUNTEST=sys.argv[4]
VALGRIND=sys.argv[5]

HOSTNAME = socket.gethostname()
tmpdir = os.path.join( BINDIR, '..', 'testing', 'temporary', HOSTNAME, RUNTEST )
if not os.path.isdir(tmpdir):
    os.makedirs(tmpdir)

BASELINE = os.path.join( DATADIR, 'testing', 'baseline', RUNTEST )

os.chdir( os.path.join( DATADIR, 'testing' ,'inputs' ) )

def run(cmd):
    print 'pushd ' + os.getcwd() + '; ' + BINDIR + '/' + cmd + '; popd'
    return os.system( VALGRIND + ' ' + BINDIR + '/' + cmd + " > " + os.path.join( tmpdir, "stdout" ) )

def check_result(name):
    baseline = os.path.join( BASELINE, name )
    result = os.path.join( tmpdir, name )

    if os.path.isfile(result):
        f = io.open(result, 'rb')
        result_content = f.readlines()
        f.close()
    else:
        if os.path.isfile(result + '.gz'):
            f = gzip.open(result + '.gz', 'rb')
            result_content = f.readlines()
            f.close()
        else:
            print 'Results file ' + result + ' does not exist'
            sys.exit(1)

    if os.path.isfile(baseline):
        f = io.open(baseline, 'rb')
        baseline_content = f.readlines()
        f.close()
    else:
        if os.path.isfile(baseline + '.gz'):
            f = gzip.open(baseline + '.gz', 'rb')
            baseline_content = f.readlines()
            f.close()
        else:
            print 'Baseline file ' + baseline + ' does not exist'
            sys.exit(1)

    diff = list( difflib.unified_diff( result_content, baseline_content, fromfile=result, tofile=baseline) )
    if len( diff ) == 0:
        return 0
    sys.stdout.writelines(diff)
    return 1

def check_results(list):
    for r in list:
        if not check_result(r) == 0:
            return 1
    return 0

if RUNTEST=='AffineRegistrationMrMrMSD':
    run("registration -i --dofs 6,9 --msd --match-histograms -o " + tmpdir + " pat001_mr_T1.hdr pat002_mr_T2.hdr")
    exit( check_results(['registration']) )
elif RUNTEST=='xml_film':
    run("film --xml")
    exit( check_results(['stdout']) )
