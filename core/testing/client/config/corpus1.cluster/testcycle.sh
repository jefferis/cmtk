#!/bin/sh

export PATH="/bin:/usr/bin:/usr/local/bin:/opt/csw/bin"

lockfile=${HOME}/testcycle.lock
if test -f ${lockfile}; then
        exit 1
fi

touch ${lockfile}
trap "rm -f ${lockfile}; exit" INT TERM EXIT

/opt/csw/bin/svn update
tests=`ls ctest-*.cmake`

Xpid=""
if [ "${DISPLAY}" == "" ]; then
	/usr/X11/bin/Xvfb :1 -screen 0 1024x768x24 -ac &
	Xpid=$!
	export DISPLAY=:1
fi

for t in ${tests}; do
        /usr/local/bin/ctest -S ${t}
done

if [ "${Xpid}" != "" ]; then
  kill ${Xpid}
fi

rm -rf /tmp/tmp.*
rm ${lockfile}
