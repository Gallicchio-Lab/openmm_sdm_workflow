#!/bin/bash

. setup-settings.sh || exit 1

datafilename=repl.cycle.totE.potE.temp.lambda.ebind.lambda1.lambda2.alpha.u0.w0.dat

for lig in ${ligands} ; do

    jobname=${set_basename}-${lig}
    jobdir=${work_dir}/complexes/${jobname}

    cd ${jobdir} || exit 1
    rm -f ${datafilename} || exit 1
    if [ -d r0 ] ; then

	echo "free energy analysis for ligand ${lig}"

	for repl_dir in r? r?? ; do
	    
	    repl=${repl_dir#r} || exit 1
	    cd ${jobdir}/${repl_dir} || exit 1
	    rm -f ${datafilename} || exit 1
	    n=0
	    for f in `/bin/ls -v ${jobname}*.out` ; do
		while read -r line ; do
		    words=(${line// / })
		    #lambda,lambda1,lambda2,alpha,u0,w0,potE,bindE
		    # 0       1       2       3   4   5  6     7
 		    echo "${repl} ${n} ${words[4]} ${words[4]} 300.0 ${words[0]} ${words[7]} ${words[1]} ${words[2]} ${words[3]} ${words[4]} ${words[5]} " >> ${datafilename}
		    n=`expr $n + 1 `
		done < $f
	    done
	    
	done

	cd ${jobdir} || exit 1
	
	rm -f result.log Rplots.pdf p-l*.dat || exit 1
    
	( cat `/bin/ls -v r*/${datafilename}` | awk "\$2 > ${discard_samples} {print \$0} " > ${datafilename} ) || exit 1
    
	R CMD BATCH uwham_analysis.R || exit 1
	res=`grep 'DGb =' result.log` || exit 1

	minc=`awk 'BEGIN{min=100000}; $1 ~ /[0-9]+/{if( $11 < min) min = $11};END{print min}' < ${jobname}_stat.txt` || exit 1
	maxc=`awk 'BEGIN{max=-1}; $1 ~ /[0-9]+/{if( $11 > max) max = $11};END{print max}' < ${jobname}_stat.txt` || exit 1

	echo "${jobname} " ${res} " min/max cycles:" $minc $maxc || exit 1
    fi

done
