#!/bin/bash

#Converts cms file from opls to amber parameters, often protein receptor
#
# usage: bash ff_opls_to_amber.sh struct.cms
#
#needs SCHRODINGER and AMBERHOME env variables
#
#write .maegz file in PDB residue order, cap ends or make sure N and C terminals
#are compatible with AMBER's templates
#

basename=${1%.cms}
cp ${basename}.maegz ${basename}.mae.gz || exit 1
rm -f ${basename}.mae
gunzip ${basename}.mae.gz || exit 1
cp ${basename}.cms ${basename}-opls.cms || exit 1
$SCHRODINGER/run -FROM desmond viparr.py -f amber99SB-ILDN -f tip3p ${basename}.mae ${basename}-out.cms || exit 1
$SCHRODINGER/run -FROM desmond build_constraints.py ${basename}-out.cms ${basename}-amber.cms || exit 1
cp ${basename}-amber.cms ${basename}.cms
