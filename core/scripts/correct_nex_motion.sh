#!/bin/sh

bindir=@bindir@
registration=${bindir}/registration
reformatx=${bindir}/reformatx
average_images=${bindir}/average_images

tmpdir=`mktemp -d`
if [ "${tmpdir}" == "" ]; then
    echo "ERROR: cannot create tmpdir"
    exit
fi

outfile=$1
reffile=$2
shift 2
fltfiles="$*"

images=${reffile}
for f in ${fltfiles}; do
    name=`basename ${f}`

    echo "Registration: ${name}"
    xform=${tmpdir}/${name}.list
    ${registration} -i -o ${xform} -a 0.01 --dofs 6 -e 2 ${reffile} ${f}

    echo "Reformatting: ${name}"
    nrrd=${tmpdir}/reformat/${name}.nrrd
    ${reformatx} --sinc-cosine --sinc-window-width 5 -o ${nrrd} --floating ${f} ${reffile} ${xform}
    images="${images} ${nrrd}"
done

echo "Averaging"
${average_images} -o ${outfile} ${images}
echo "Done."

##rm -rf ${tmpdir}
