#!/bin/bash

#usage: bash analyze-equilibration.sh | awk '{print $2, $5, $7}' > t4l-3iodotoluene-ilogistic-b-ratsc-a16-p3-equil.dat

#assumes that analyze.sh has been already executed with discard_samples=0

set_basename=t4l
work_dir=~/Dropbox/sims/bind_phaset/t4l
scripts_dir=${work_dir}/scripts_ilogistic
datafilename=repl.cycle.totE.potE.temp.lambda.ebind.lambda1.lambda2.alpha.u0.w0.dat
#ligands="3iodotoluene-linear-b-tanhsc 3iodotoluene-linear-b-ratsc-a16 3iodotoluene-ilogistic-b-ratsc-a16-p3"
ligands="3iodotoluene-linear-b-ratsc-a16"

for lig in $ligands ; do


for discard_samples in `seq 0 25 1000` ; do


    jobname=${set_basename}-${lig}
    jobdir=${work_dir}/complexes/${jobname}

    cd ${jobdir}
    
    rm -f result.log Rplots.pdf p-l*.dat .RData
    
    cat `/bin/ls -v r*/${datafilename}` | awk "\$2 > ${discard_samples} {print \$0} " > ${datafilename}
    
    R CMD BATCH uwham_analysis.R
    res=`grep 'DGb =' result.log`

    echo "${jobname} ${discard_samples}" ${res}
done

done
