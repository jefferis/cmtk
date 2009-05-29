#!/bin/sh

HOSTNAME=`hostname`
BUILDS=`ls config/${HOSTNAME}`

setup_build()
{
    local build=$1
    if [ ! -d source/${build} ]; then
	svn co https://neuro.sri.com/repos/igs/trunk/ source/${build}
    fi
}

config_build()
{
    local build=$1

    if [ ! -d builds/${build} ]; then
	mkdir -p builds/${build}
    fi

    pushd builds/${build}
    cmake -C ../../config/${HOSTNAME}/${build} ../../source/${build}
    popd
}

CTEST=/usr/bin/ctest

run_ctest()
{
    local build=$1
    shift
    local protocol="$*"

    pushd builds/${build}
    ${CTEST} -D ${protocol}
    popd
}

for b in ${BUILDS}; do
    if [ ! -d builds/${b} ]; then
	setup_build ${b}
	config_build ${b}
	run_ctest $b Experimental
    else
	config_build ${b}
	run_ctest $b Continuous
    fi
done
