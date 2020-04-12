#!/bin/bash

#load settings
. setup-settings.sh

export SCHRODINGER=${schrodinger}
export PATH=$PATH:${msys_path}/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${msys_path}/lib
export PYTHONPATH=$PYTHONPATH:${msys_path}/lib/python


des_builder_cmd="
task {
  task = "desmond:auto"
}

build_geometry {
  distil_solute = false
  neutralize_system = false
  rezero_system = false
  solvate_system = false
}

assign_forcefield {
}
"

des_builder_file=des_builder.msj

cd ${work_dir} || exit 1

#setup receptor
echo "Assigning force field parameters to receptor ..."
jobdir=${work_dir}/receptor
cp ${scripts_dir}/add_agbnp2.py ${jobdir}/  && cp ${scripts_dir}/${agbnpparam} ${jobdir}/agbnp2.param  && cd ${jobdir} || exit 1

printf '%s\n' "${des_builder_cmd}" > ${des_builder_file}   && $SCHRODINGER/utilities/multisim -JOBNAME ${set_basename}_receptor -m des_builder.msj ${receptor}.maegz -o ${receptor}.cms -maxjob 1 -WAIT || exit 1

if [ ! -z ${amberhome+x} ]; then
    export AMBERHOME=${amberhome}
    echo "Assigning AMBER parameters to receptor ..."
    bash ${scripts_dir}/ff_opls_to_amber.sh ${receptor}.cms || exit 1
fi

echo "Converting receptor to .dms format ..."
mae2dms ${receptor}.cms ${receptor}.dms                 || exit 1

echo "Assigning AGBNP parameters to receptor ..."
$SCHRODINGER/run add_agbnp2.py ${receptor}.dms          || exit 1

echo "Adding restraints to receptor ..."
python ${scripts_dir}/add_fbh_restraints.py -r "${rest_receptor_sql}" -t "${rest_receptor_tol}" -k "${rest_receptor_kf}"  -f "${receptor}.dms" || exit 1

complex_dir=${work_dir}/complexes
mkdir -p ${complex_dir}                                    || exit 1

#assign force field parameters to ligands and sets up complexes
for lig in ${ligands} ; do
  jobname=${set_basename}-${lig}
  jobdir=${complex_dir}/${jobname}
  mkdir -p ${jobdir}                                    || exit 1

  cp ${scripts_dir}/add_agbnp2.py ${jobdir}/ && cp ${scripts_dir}/${agbnpparam} ${jobdir}/agbnp2.param || exit 1

  cp ${scripts_dir}/runopenmm ${jobdir}/ && cp ${scripts_dir}/nodefile ${jobdir}/                      || exit 1
  
  cp ${work_dir}/receptor/${receptor}.dms ${jobdir}/${jobname}_rcpt.dms && cp ${work_dir}/ligands/${lig}.maegz ${jobdir}/ || exit 1

  cd ${jobdir} || exit 1

  echo "Assigning force field parameters to ligand ${lig} ..."
  echo ${des_builder_cmd} > ${des_builder_file} && $SCHRODINGER/utilities/multisim -JOBNAME ${set_basename}_${jobname} -m des_builder.msj ${lig}.maegz -o ${lig}.cms -maxjob 1 -WAIT || exit 1

  if [ ! -z ${amberhome+x} ]; then
    export AMBERHOME=${amberhome}
    echo "Assigning GAFF parameters to ${lig} ..."
    bash ${scripts_dir}/ff_opls_to_gaff.sh ${lig}.cms || exit 1
  fi
  
  echo "Converting ${lig} to .dms format ..."
  mae2dms ${lig}.cms ${jobname}_lig.dms || exit 1

  echo "Assigning AGBNP parameters to ${lig} ..."
  $SCHRODINGER/run add_agbnp2.py ${jobname}_lig.dms || exit 1
done

#runs workflow on each complex to set up openmm asyncre and openmm calculations
for lig in ${ligands} ; do
    
    echo "Running sdm workflow for ligand ${lig} ..."
    jobname=${set_basename}-${lig}
    jobdir=${complex_dir}/${jobname}
    cd ${jobdir} || exit 1
    
    echo "Writing workflow control file ..."
    sed "s/<JOBNAME>/${jobname}/" < ${scripts_dir}/${cntltmpl} > ${jobdir}/${jobname}.cntl || exit 1

    cp ${scripts_dir}/openmm-driver-template.py ${scripts_dir}/openmm-mintherm-template.py ${scripts_dir}/uwham_analysis-template.R ${jobdir}/ || exit 1
    
    echo "Running workflow ..."
    python ${scripts_dir}/sdm_workflow_openmm.py ${jobname}.cntl || exit 1
    
done

#minimization and thermalization of each complex
for lig in ${ligands} ; do
      echo "Running minimization thermalization for complex with ${lig} ..."
      jobname=${set_basename}-${lig}
      jobdir=${complex_dir}/${jobname}
      cd ${jobdir} || exit 1
      ./runopenmm ${jobname}_mintherm.py > ${jobname}_mintherm.log || exit 1
done
