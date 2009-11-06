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
tmpdir = BINDIR + '/../testing/temporary/' + HOSTNAME + '/' + RUNTEST
if not os.path.isdir(tmpdir):
    os.makedirs(tmpdir)

BASELINE = DATADIR + '/testing/baseline/' + RUNTEST

os.chdir(DATADIR + '/testing/inputs')

def run(cmd):
    print 'pushd ' + os.getcwd() + '; ' + BINDIR + '/' + cmd + '; popd'
    return os.system(VALGRIND + ' ' + BINDIR + '/' + cmd)

def check_result(name):
    baseline = BASELINE + '/' + name
    result = tmpdir + '/' + name

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

    diff = difflib.unified_diff( result_content, baseline_content, fromfile=result, tofile=baseline)
    for line in diff:
        sys.stdout.write(line) 

def check_results(list):
    for r in list:
	check_result(r)

if RUNTEST=='AffineRegistrationMrMrMSD':
    run("registration -i --dofs 6,9 --msd --match-histograms -o " + tmpdir + " pat001_mr_T1.hdr pat002_mr_T2.hdr")
    check_result('registration')
