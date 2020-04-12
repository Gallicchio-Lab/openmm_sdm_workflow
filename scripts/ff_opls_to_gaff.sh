#!/bin/bash

#Converts cms file from opls to gaff parameters
#
# usage: bash ff_opls_to_gaff.sh struct.cms
#
#needs SCHRODINGER and AMBERHOME env variables
#write the .maegz file using v1 mae 2012 format

basename=${1%.cms}
cp ${basename}.maegz ${basename}.mae.gz || exit 1
rm -f ${basename}.mae
gunzip ${basename}.mae.gz || exit 1
cp ${basename}.cms ${basename}-opls.cms || exit 1
$SCHRODINGER/utilities/mol2convert -imae ${basename}.mae -omol2 ${basename}.mol2 -n 1:1 || exit 1 
antechamber -fi mol2 -fo mol2 -i ${basename}.mol2 -o ${basename}-p.mol2 -c bcc -at gaff || exit 1
parmchk2 -i ${basename}-p.mol2 -o ${basename}.frcmod -f mol2 || exit 1
$SCHRODINGER/run -FROM desmond ff_amber_to_viparr.py -p $AMBERHOME/dat/leap/parm/gaff.dat -p ${basename}.frcmod -t ${basename}-p.mol2 -f ligand GF1 || exit 1
$SCHRODINGER/run -FROM desmond viparr.py -d GF1 ${basename}.mae ${basename}-out.cms || exit 1
$SCHRODINGER/run -FROM desmond build_constraints.py ${basename}-out.cms ${basename}-gaff.cms || exit 1
cp ${basename}-gaff.cms ${basename}.cms
