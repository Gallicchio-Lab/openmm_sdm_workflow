
#basename for jobs
set_basename=t4l

#path of working directory
#work_dir=$HOME/${set_basename}
work_dir=/home/users/rajat/Dropbox/sims/multiple_context_sims

#path of script directory
scripts_dir=${work_dir}/scripts

#schrodinger installation path
#schrodinger=$HOME/schrodinger/Desmond_Maestro_2017.4
schrodinger=/home/software/schrodinger/Desmond_Maestro_2017.4

#msys installation path
#msys_path=$HOME/local
msys_path=/home/users/rajat/Dropbox/md_software/msys_ubuntu16

#vmd installation path
#vmd_path=/usr/local
vmd_path=/usr/local

#basename of the receptor (.maegz format expected)
#receptor=t4l99a
receptor=t4l99a

#basenames of the ligands (.maegz format expected)
#ligands="benzene toluene 3iodotoluene"
ligands=benzene

#control template file
cntltmpl=sdm_workflow_template.cntl

#agbnp parameter file
agbnpparam=agbnp2.param.agbnp_plugin

#positional restraints for the receptor
#rest_receptor_sql=$'name GLOB \'CA\''
#rest_receptor_kf=25.0
#rest_receptor_tol=0.75
rest_receptor_sql=$'name GLOB \'CA\''

#rest_receptor_kf=<force constant in kcal/mol/angstrom^2>
rest_receptor_kf=25.0

#rest_receptor_tol=<tolerance in angstrom>
rest_receptor_tol=0.75

#how many samples to discard for equilibration
#discard_samples=10
discard_samples=0
