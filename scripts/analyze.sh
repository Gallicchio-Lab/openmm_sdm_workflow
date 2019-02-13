#!/bin/bash

set_basename=t4l
work_dir=/home/emilio/Dropbox/sims/bind_phaset/t4l
scripts_dir=${work_dir}/scripts_ilogistic
discard_cycles=10
datafilename=repl.cycle.totE.potE.temp.lambda.ebind.lambda1.lambda2.alpha.u0.w0.dat

for lig in 3iodotoluene ; do

    jobname=${set_basename}-${lig}
    jobdir=${work_dir}/complexes/${jobname}-ilogistic

    cd ${jobdir}
    rm -f ${datafilename}
    for repl in `seq 0 15` ; do 
	cd r$repl
	rm -f ${datafilename}
	n=0
	for f in `/bin/ls -v ${jobname}_*.out` ; do 
	    cycle=${f%.out} ; cycle=${cycle#${jobname}_}
	    if [ $cycle -gt ${discard_cycles} ] ; then 
		while read -r line ; do
		    words=(${line// / })
		    #lambda,lambda1,lambda2,alpha,u0,w0,potE,bindE
		    # 0       1       2       3   4   5  6     7
 		    echo "${repl} ${n} ${words[4]} ${words[4]} 300.0 ${words[0]} ${words[7]} ${words[1]} ${words[2]} ${words[3]} ${words[4]} ${words[5]} " >> ${datafilename}
		    n=`expr $n + 1 `
		done < $f
	    fi
	    
	done 
	cd ..
    done

    rm -f result.log Rplots.pdf p-l*.dat
    
    cat `/bin/ls -v r*/${datafilename}` | awk "\$2 > ${discard_cycles} {print \$0} " > ${datafilename}
    
    R CMD BATCH ${scripts_dir}/uwham_analysis.R
    res=`grep 'DGb =' result.log`

    minc=`awk 'BEGIN{min=100000}; $1 ~ /[0-9]+/{if( $11 < min) min = $11};END{print min}' < ${jobname}_stat.txt`
    maxc=`awk 'BEGIN{max=-1}; $1 ~ /[0-9]+/{if( $11 > max) max = $11};END{print max}' < ${jobname}_stat.txt`

    echo "${jobname} " ${res} " min/max cycles:" $minc $maxc
done
