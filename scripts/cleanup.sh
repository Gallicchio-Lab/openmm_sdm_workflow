#!/bin/bash

set_basename=hivpr
work_dir=~/Dropbox/src/openmm_sdm_workflow
scripts_dir=${work_dir}/scripts
catdcdpath=~/Dropbox/src/catdcd/LINUXAMD64/bin/catdcd4.0
ligands="1 2 3"


#must have catdcd in PATH
export PATH=$PATH:${catdcdpath}

for lig in ${ligands} ; do

    jobname=${set_basename}-${lig}
    jobdir=${work_dir}/complexes/${jobname}-ilogistic

    cd ${jobdir}
    for repl in `seq 0 15` ; do
	cd r${repl}
	python ${scripts_dir}/cleanup.py ${jobname}
	cd ..
    done
done
