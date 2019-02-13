#!/bin/bash

set_basename=hivpr
work_dir=/home/users/egallicchio/Dropbox/src/openmm_sdm_workflow
bedam_workflow_dir=/home/users/egallicchio/Dropbox/src/bedam_workflow
scripts_dir=${work_dir}/scripts
impacthome=/home/software/academic-impact-70501
schrodinger=/home/software/schrodinger/Desmond_Maestro_2013.3

cntltmpl=bedam_workflow_template_openmm.cntl
rcpt_cntltmpl=workflow_receptor_template.cntl
numcores=1
numthreads=1

receptor=1hpx
agbnpparam=agbnp2.param.agbnp_plugin


export PYTHONPATH=$PYTHONPATH:${bedam_workflow_dir}
export SCHRODINGER=${schrodinger}
export IMPACTHOME=${impacthome}

cd ${work_dir}

<<EOF
# #setup receptor
jobdir=${work_dir}/receptor 
cp ${scripts_dir}/paramstd.dat ${jobdir}/
cp ${scripts_dir}/add_agbnp2.py ${jobdir}/
cp ${scripts_dir}/${agbnpparam} ${jobdir}/agbnp2.param
cp ${scripts_dir}/add_hct.py ${jobdir}/
cp ${scripts_dir}/hct.param ${jobdir}/
sed "s^<IMPACTHOME>^${impacthome}^;s/<NUMCORES>/${numcores}/;s/<NUM_THREADS>/${numthreads}/;s/<RECEPTOR>/${receptor}.maegz/;s/<JOBNAME>/${jobname}/" < ${scripts_dir}/${rcpt_cntltmpl} > ${jobdir}/${receptor}.cntl
# echo "Running bedam workflow for receptor ..."
jobname=${receptor}
cd ${jobdir}
$SCHRODINGER/run ${scripts_dir}/tempt_asyncre_openmm.py ${jobname}.cntl > ${jobname}_workflow.log
EOF

#set up complexes
complex_dir=${work_dir}/complexes
mkdir ${complex_dir}


for lig in 3 4 ; do
    
  jobname=${set_basename}-${lig}
  jobdir=${complex_dir}/${jobname}
  mkdir -p ${jobdir}

  cp ${scripts_dir}/paramstd.dat ${jobdir}/
  cp ${scripts_dir}/add_agbnp2.py ${jobdir}/
  cp ${scripts_dir}/${agbnpparam} ${jobdir}/agbnp2.param

  cp ${scripts_dir}/add_hct.py ${jobdir}/
  cp ${scripts_dir}/hct.param ${jobdir}/

  cp ${scripts_dir}/runopenmm ${jobdir}/
  cp ${scripts_dir}/nodefile ${jobdir}/

  cp ${scripts_dir}/SDMUtils.py ${jobdir}/
  
#  cp ${work_dir}/receptor/${receptor}.dms ${jobdir}/${jobname}_rcpt.dms
  cp ${scripts_dir}/${receptor}.maegz ${jobdir}/

  cp ${work_dir}/ligands/${lig}.maegz ${jobdir}/
  sed "s^<IMPACTHOME>^${impacthome}^;s/<NUMCORES>/${numcores}/;s/<NUM_THREADS>/${numthreads}/;s/<RECEPTOR>/${receptor}.maegz/;s/<LIGAND>/${lig}.maegz/;s/<JOBNAME>/${jobname}/" < ${scripts_dir}/${cntltmpl} > ${jobdir}/${jobname}.cntl
done


# #run bedam workflow

for lig in 3 4 ; do

    echo "Running bedam workflow for ligand ${lig} ..."
    jobname=${set_basename}-${lig}
    jobdir=${complex_dir}/${jobname}
    cd ${jobdir}
    $SCHRODINGER/run ${scripts_dir}/bedam_workflow_openmm.py ${jobname}.cntl egallicc > ${jobname}_bedam_workflow.log 
    cd ..   
done

<<EOF
for lig in `seq 1 2` ; do
# # ##for lig in 23 ; do

      echo "Running mintherm for ligand ${lig} ..."
      jobname=${set_basename}-${lig}
      jobdir=${complex_dir}/${jobname}
      cd ${jobdir}
      ./runopenmm ${jobname}_mintherm.py > ${jobname}_mintherm.log
      cd ..   
done
EOF

