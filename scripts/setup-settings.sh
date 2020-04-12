
#basename for jobs
set_basename=<the name for this project>

#path of working directory
#work_dir=$HOME/${set_basename}
work_dir=<the path to the work directory for this project>

#path of script directory
scripts_dir=${work_dir}/scripts

#schrodinger installation path
#schrodinger=$HOME/schrodinger/Desmond_Maestro_2017.4
schrodinger=<path to schrodinger installation directory>

#ambertools installation path,
#if enabled it will try to assign AMBER/GAFF/AM1-BCC parameters 
#amberhome=/opt/amber

#msys installation path
#msys_path=$HOME/local
msys_path=<path to msys installation directory>

#vmd installation path
#vmd_path=/usr/local
vmd_path=<path to vmd installation directory>

#basename of the receptor (.maegz format expected)
#receptor=t4l99a
receptor=<basename of the .maegz receptor file>

#basenames of the ligands (.maegz format expected)
#ligands="benzene toluene 3iodotoluene"
ligands=<list of ligands>

#control template file
cntltmpl=sdm_workflow_template.cntl

#agbnp parameter file
agbnpparam=agbnp2.param.agbnp_plugin

#positional restraints for the receptor
#rest_receptor_sql=$'name GLOB \'CA\''
#rest_receptor_kf=25.0
#rest_receptor_tol=0.75
rest_receptor_sql=<sqlite atom selection string>
rest_receptor_kf=<force constant in kcal/mol/angstrom^2>
rest_receptor_tol=<tolerance in angstrom>

#how many samples to discard for equilibration
#discard_samples=10
discard_samples=<num samples to discard>
