#!/bin/bash
set_basename=mcl1
#work_dir='/home/users/egallicchio/Dropbox/sims/MCL1/mcl1'
#bedam_workflow_dir='/home/users/egallicchio/Dropbox/src/bedam_workflow'
work_dir='/nfs/egallicchio2/Dropbox/sims/MCL1/mcl1'
bedam_workflow_dir='/nfs/egallicchio2/Dropbox/src/bedam_workflow'
scripts_dir=${work_dir}/scripts
impacthome='/home/software/academic-impact-70500'
impacthome_sed='\/home\/software\/academic-impact-70500'
schrodinger=/home/software/schrodinger/Desmond_Maestro_2013.3
ligandcore="name GLOB 'CORE'"

cntltmpl=bedam_workflow_template_openmm.cntl
numcores=4
numthreads=1

receptor=mcl1
agbnpparam=agbnp2.param.agbnp_plugin

cd ${work_dir}

complex_dir=${work_dir}/complexes
mkdir ${complex_dir}

for lig in `cat ${work_dir}/ligands/ligand_list` ; do
#for lig in 23 ; do
    
  jobname=${set_basename}-${lig}
  jobdir=${complex_dir}/${jobname}

  cp ${scripts_dir}/driver-template.py ${jobdir}/${jobname}.py
done
