#!/bin/bash

. setup-settings.sh || exit 1

#to locate catdcd
catdcdpath=`/bin/ls -1 -d ${vmd_path}/lib/vmd/plugins/LINUXAMD64/bin/catdcd*` || exit 1
export PATH=$PATH:${catdcdpath}

for lig in ${ligands} ; do

    jobname=${set_basename}-${lig}
    jobdir=${work_dir}/complexes/${jobname}
    cd ${jobdir} || exit 1

    if [ -d r0 ] ; then 
	for repl_dir in r? r?? ; do
	    cd ${jobdir}/${repl_dir} || exit 1
	    python ${scripts_dir}/cleanup.py ${jobname} || exit 1
	done
    fi
done
