# OpenMM SDM Tutorial

1. [Introduction](#introduction)  
2. [Gathering the Software Tools](#gathering-the-software-tools)
3. [Setup the Simulations](#setup-the-simulations)
4. [Run the ASyncRE simulations](#run-the-asyncre-simulations)
5. [Analysis](#analysis)
6. [Appendix](#appendix)


## Introduction

This tutorial walks you over the steps to setup and run absolute binding free energy calculations for a series of ligands of a mutant of T4 lysozyme using the Single Decoupling (SDM) OpenMM's plugin.

SDM is based on an alchemical process in which the ligand is progressively transferred from the solution environment to the receptor binding site. SDM employs an implicit description of the solvent (here we use the AGBNP model) which allows it to avoid the intermediate vacuum state necessary for absolute binding free energy calculations with explicit solvation (Double Decoupling). Hence SDM requires only one free energy calculation as opposed to two with Double Decoupling. Hamiltonian Parallel Replica Exchange with OpenMM is used for conformational sampling.

For a background on the alchemical process and the conformational sampling techniques used with SDM, consult our publications. Here are the key ones:

* [Binding energy distribution analysis method (BEDAM) for estimation of protein-ligand binding affinities](http://www.compmolbiophysbc.org/publications#bedam_2010)
* [Analytic Model of the Free Energy of Alchemical Molecular Binding](http://www.compmolbiophysbc.org/publications#analytical_theory_2018)
* [Recent Theoretical and Computational Advances for Modeling Protein-Ligand Binding Affinities](http://www.compmolbiophysbc.org/publications#pubs_advprot_2011)
* [Theory of binless multi-state free energy estimation with applications to protein-ligand binding](http://www.compmolbiophysbc.org/publications#uwham)
* [Asynchronous Replica Exchange Software for Grid and Heterogeneous Computing](http://www.compmolbiophysbc.org/publications#asyncre_software_2015)
* [Efficient Gaussian Density Formulation of Volume and Surface Areas of Macromolecules on Graphical Processing Units](http://www.compmolbiophysbc.org/publications#gaussvol_2017)
* [AGBNP, an analytic implicit solvent model suitable for molecular dynamics simulations and high-resolution modeling](http://www.compmolbiophysbc.org/publications#AGBNP1)


### Molecular Systems

For this tutorial, we will consider the complexes between benzene, toluene, and 3-iodotoluene and the L99A mutant of the T4 lysozyme receptor (PDB id: 4W53). The systems were prepared using the free academic version of [Maestro/Desmond](https://www.deshawresearch.com/downloads/download_desmond.cgi). Toluene, water molecules and other bound ligands were removed from the 4W53 receptor structure. To reduce the size of the system and speed-up the calculations, residues 1 through 71 were removed. The receptor was then processed with Protein Preparation tool in Maestro to add hydrogen atoms and cap the termini. OpenMM prefers protein structures in which atoms belonging to the same residue are listed consecutively. To do so the processed receptor structure was saved in PDB format and read back into Maestro.

Benzene and 3-iodotoluene were built with academic Maestro using bound toluene as a template. Bond order and formal charges were adjusted and  hydrogen atoms added in Maestro. In most actual applications, the starting bound conformations of ligands are obtained using molecular docking, using Glide for example. The commercial version of the Schrodinger's Suite includes a easy-to-follow tutorial on the steps to prepare protein receptors and ligands for docking. The end result is the receptor structure and a list of prepared and docked ligand structures in the Maestro project table as in the present case.  

The receptor and the ligands were exported from Maestro in separate folders (`ligands/` and `receptor/`) in `.maegz` V1 formatted files (see below).   

## Gathering the Software Tools

SDM calculations require several pieces of software:

1. Conda 2
2. OpenMM 7.2 or later
3. Desmond DMS file reader for OpenMM
4. Desmond (for force field parameter assignment)
5. Msys
6. SDM workflow
7. SDM OpenMM plugin
8. AGBNP OpenMM plugin
9. ASyncRE for OpenMM
10. UWHAM for R
11. VMD for visualization 

We summarize here the steps to install them and configure them. The steps below were tested on a 16.04 Ubuntu system.

### Conda 2

Python bindings for OpenMM and other libraries are best stored in a Conda environment. We recommend the Miniconda 2 environment:

```
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```

In this tutorial we assume that Miniconda 2 is installed under `$HOME/miniconda2`. Adjust as needed.

Activate the conda environment

```
source ~/miniconda2/bin/activate
```

### OpenMM

We recommend installing OpenMM from [source](https://github.com/pandegroup/openmm) since most of the steps required are the same as for building the SDM-related plugins (see below). However binary installations of OpenMM should probably also work. Detailed building instructions for OpenMM are [here](http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code). SDM requires an OpenCL platform with GPUs from NVIDIA (CUDA) or AMD, which we assume are in place.

These are the steps we used to build OpenMM 7.3.1 on an Ubuntu 16.04 system:

```
mkdir $HOME/devel
cd $HOME/devel
wget https://github.com/pandegroup/openmm/archive/7.3.1.tar.gz
tar zxvf 7.3.1.tar.gz
conda install cmake=3.6.3 swig numpy
conda install -c conda-forge doxygen
mkdir build_openmm
cd build_openmm
ccmake -i ../openmm-7.3.1
```

Hit `c` (configure) until all variables are correctly set. Set `CMAKE_INSTALL_PREFIX` to point to the openmm installation directory. Here we assume `$HOME/local/openmm-7.3.1`. If an OpenCL platform is detected (such as from a NVIDIA CUDA installation) it will be enabled automatically. For the present purposes the CUDA platform is optional. Hit `g` to generate the makefiles and `q` to exit `ccmake`, then:

```
make install
make PythonInstall
```

The OpenMM libraries and header files will be installed in `$HOME/local/openmm-7.3.1`  

### Desmond File Reader for OpenMM

SDM uses Desmond DMS-formatted files. OpenMM includes a python library to load molecular files in DMS Desmond format however, the version of the Desmond file reader required for SDM is not yet included in the latest OpenMM sources. To patch the 7.3.1 OpenMM installation above with the latest DMS file reader do the following:

```
cd $HOME/devel
wget https://raw.githubusercontent.com/egallicc/openmm/master/wrappers/python/simtk/openmm/app/desmonddmsfile.py
cp desmonddmsfile.py $HOME/devel/openmm-7.3.1/wrappers/python/simtk/openmm/app/
cd $HOME/devel/build_openmm
make install
make PythonInstall
```

The `sqlitebrowser` application is very useful to inspect DMS files.

### Desmond

In this tutorial we will use academic Desmond as part of the Schrodinger's environment to assign OPLS2005 force field parameters. Download and install [Maestro/Desmond](https://www.deshawresearch.com/downloads/download_desmond.cgi). If you have it, the commercial version of Maestro/Desmond also works, of course.

### Msys

`msys` is a software package developed by DE Shaw Research to manipulate molecular structures. We use it to convert Maestro-formatted files to DMS-formatted files. Installation instructions are in [Appendix A](#msys) below. We assume that `msys` is installed in `$HOME/local` and that executables and libraries there are in the shell search paths. See [Appendix A](#msys).


### SDM Workflow

It's this package

```
cd $HOME/devel
git clone https://github.com/egallicc/openmm_sdm_workflow.git
```

For this tutorial we will use the structures stored under `openmm_sdm_workflow/tutorial` and modify the files under `openmm_sdm_workfklow/tutorial/scripts` to setup the simulations.

### SDM OpenMM plugin

This is a plugin that implements a customized Langevin integrator to compute the alchemical potential in SDM.

```
cd $HOME/devel
git clone https://github.com/rajatkrpal/openmm_sdm_plugin.git
```

then follow the [installation instructions](https://github.com/rajatkrpal/openmm_sdm_plugin/blob/master/README.md) for this package. Point the `CMAKE_INSTALL_PREFIX` and `OPENMM_DIR` to the OpenMM installation directory (`$HOME/local/openmm-7.3.1` in this example).


### AGBNP OpenMM plugin

This is a plugin that implements the AGBNP implicit solvent model in OpenMM

```
cd $HOME/devel
git clone https://github.com/egallicc/openmm_agbnp_plugin.git
```

then follow the [installation instructions](https://github.com/egallicc/openmm_agbnp_plugin/blob/master/README.md) for this package. Point the `CMAKE_INSTALL_PREFIX` and `OPENMM_DIR` to the OpenMM installation directory (`$HOME/local/openmm-7.3.1` in this example).

### ASyncRE for OpenMM

ASyncRE is a package written in python to perform replica exchange simulations in asynchronous mode across a wide range of computational devices and grids. The specific version used by SDM distributes OpenMM jobs across GPU compute servers through passwordless `ssh`. To install it do:

```
conda install numpy configobj paramiko
cd $HOME/devel
git clone https://github.com/baofzhang/async_re-openmm.git
cd async_re-openmm
python setup.py install
```

Here is a [guide](https://www.digitalocean.com/community/tutorials/how-to-set-up-ssh-keys-on-ubuntu-1604) to set up password-less `ssh` across a farm of computers. ASyncRE does not require that the compute servers share a common user space and a shared filesystem.

### UWHAM for R

The unbinned weighted histogram analysis method ([UWHAM](http://www.compmolbiophysbc.org/publications#uwham)) is an implementation of the MBAR thermodynamic reweighting method for R. It is used in this tutorial to compute binding free energies from the output of SDM calculations. The UWHAM package is part of the CRAN archive. To install it do:

```
R
install.packages("UWHAM")
```

### VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is an indispensable tool for visualization and analysis of trajectories. Also, we use the `catdcd` tool distributed with VMD to parse and concatenate .dcd trajectory files. For this tutorial we will assume that VMD is installed in `/usr/local`.


## Setup the Simulations

### Step 1: copy files and scripts

Create a simulation directory. Here we use `t4l` as an example. Then we create three directories to store the receptor and the ligands:

```
mkdir $HOME/t4l
cd $HOME/t4l
```

Copy the receptor and ligand structures for the tutorial into the `ligands` and `receptor` folders. We assume that the SDM workflow has been installed under `$HOME/devel/openmm_sdm_workflow` (see above).

```
cp -r $HOME/devel/openmm_sdm_workflow/tutorial/receptor $HOME/t4l/
cp -r $HOME/devel/openmm_sdm_workflow/tutorial/ligands $HOME/t4l/
```

Now copy the scripts and script templates into the tutorial folder:

```
cp -r $HOME/devel/openmm_sdm_workflow/tutorial/scripts .
cd scripts
```

### Step 2: set the simulation parameters in `setup-settings.sh`

Edit the `setup-settings.sh` file in the newly created `scripts` folder. In each case replace each placeholder with a specific setting using the example as a guide. For example change:

```
receptor=<basename of the .maegz receptor file>
```

with

```
receptor=t4l99a
```

Settings:

* `set_basename`: the name of the project. The name of files, simulation directories, etc. will include this name. For this tutorial set it to `t4l`.
* `work_dir`: the working directory for this project. That is where the `ligands` and `receptor` directories reside. For this tutorial set it to `$HOME/${set_basename}` which will resolve to `$HOME/t4l`
* `scripts_dir`: where the scripts are stored. The default points to `${work_dir}/scripts` which is what we want.
* `schrodinger`: your Schrodinger installation directory. For example `/opt/software/schrodinger/Desmond_Maestro_2017.4`
* `msys_path`: the msys installation directory. If you followed the instructions above it will be `$HOME/local`.
* `vmd_path`: the vmd installation directory, usually `/usr/local`.
* `receptor`: the basename of the receptor file. Set it to `t4l99a` for this tutorial. The setup script will then look for `t4l99a.maegz` in the `receptor` folder.
* `ligands`: a string with the list of basenames of the ligands. Set it to `benzene toluene 3iodotoluene`
* `cntltmpl`: the name of the template control file for the SDM workflow (see below). Leave this unchanged to `sdm_workflow_template.cntl`
* `agbnpparam`: the name of the AGBNP parameter file. Leave this unchanged to `agbnp2.param.agbnp_plugin`.
* `rest_receptor_sql`: a string with a SQLite atom selection specifying which atoms of the receptor will be restrained. SDM uses isotropic, flat-bottom harmonic potentials. For example `name GLOB 'CA'` is used to restrain C-alpha atoms. Notice that the single quotes need to be escaped in the shell. For this tutorial set `rest_receptor_sql` to `$'name GLOB \'CA\''`
* `rest_receptor_kf`: the force constant of the receptor restraints in kcal/mol/angstrom^2. Set it to `25.0`.
* `rest_receptor_tol`: the tolerance of the receptor restraints in angstrom. Set it to `0.75`.
* `discard_samples`: how many samples to discard from the beginning of the trajectory for equilibration. For this tutorial set it to 10 (but make sure to obtain more than 10 samples from each replica).

### Step 3: set the alchemical schedule and ASyncRE settings in the `sdm_workflow_template.cntl` file

Edit the `sdm_workflow_template.cntl` file in the `scripts` folder. In each case replace each placeholder with a specific setting using the example as a guide. For example change:

```
WALL_TIME = <duration of each ASyncRE simulation, in minutes>
```

with

```
WALL_TIME = 480
```

Settings:

* `TEMPERATURE`: the simulation temperature in kelvin. Set it to `300`.
* `TEMPERATURES`: list of replica exchange temperatures. Multi-dimensional replica exchange in temperature and alchemical space is currently not supported. Set this to the single temperature `300`
* `LAMBDAS`: list of alchemical lambda values in comma-separated string. Set it to `' 0.000, 0.057, 0.114, 0.171, 0.229, 0.286, 0.343, 0.400, 0.457, 0.514, 0.571, 0.629, 0.686, 0.743, 0.800, 1.000'`
* `CYCLE_TIME`: period of RE exchanges in seconds. Set it to `30`.
* `WALL_TIME`: wall-clock duration of the RE simulation for each complex in minutes. Set it to `480`.
* `PRODUCTION_STEPS`: number of MD steps for each running cycle of a replica. Set it to `20000`.
* `PRNT_FREQUENCY`: the MD period, in MD steps, for saving binding energy samples. Set it to `5000`.
* `TRJ_FREQUENCY`: the MD period, in MD steps, for saving trajectory frames. Set it to `5000`.
* `REST_LIGAND_CMLIGSQL`: sqlite selection specifying the CM atoms of the ligand. Set it to `'name GLOB 'C?''` to use the CM of the carbon atoms of the ligands.
* `REST_LIGAND_CMRECSQL`: same as above but for the receptor. Set it to `'name GLOB 'CA' AND resid IN (79,84,88,91,96,104,112,113,122,133,150)'` to use the CM determined by the C-alpha atoms of selected residues.
* `REST_LIGAND_CMKF`: force constant of flat bottom CM-CM restraint, in kcal/mol/angstrom^2. Set it to `25.0`
* `REST_LIGAND_CMTOL`: tolerance of flat bottom CM-CM restraint, in angstrom. Set it to `2.5`.
* `SOFT_CORE_UMAX`: max binding energy of soft-core potential in kcal/mol. Set it to `50.0`.
* `SOFT_CORE_ACORE`: exponent of rational soft core potential, dimensionless. Set it to `0.0625` or 1/16.


### Step 4: configure the `runopenmm` launch script

Edit the `runopenmm` launch script to reflect your environment. For example, if OpenMM is installed in `$HOME/local/openmm-7.3.1` and the corresponding python bindings are installed in the Conda environment under `$HOME/miniconda2`, the `runopenmm` should look like:

```
#!/bin/bash
export LD_LIBRARY_PATH=$HOME/local/openmm-7.3.1/lib:$HOME/local/openmm-7.3.1/lib/plugins:$LD_LIBRARY_PATH
$HOME/miniconda2/bin/python $1
```

Technically, the `runopenmm` launch script should reflect the OpenMM installation environment on the GPU compute machines (see below). In this tutorial we assume that the local machine and the remote compute machines (which could include the local machine) execute OpenMM in the same way.

### Step 5: prepare the `nodefile`

The ASyncRE system dispatches jobs to GPU computing devices on the machines listed in the `nodefile`. The format of each line is as follows:

```
<machine name> , <OpenCL platform>:<OpenCL GPU device id> , <num GPUs> , <platform name> , <remote username>, <remote temp directory>
```

For example, when using 4 GPUs on 2 servers named `compute-server-1` and `compute-server-2`, we will set the nodefile as:

```
compute-server-1,0:0,1,OpenCL,,/tmp
compute-server-1,0:1,1,OpenCL,,/tmp
compute-server-2,1:0,1,OpenCL,,/tmp
compute-server-2,1:1,1,OpenCL,,/tmp 
```

This assumes that on the first machine the OpenCL platform for the GPUs is the first platform (0), and for the second server it is the second OpenCL platform (1). If you are running this tutorial on your desktop with only one GPU, a `nodefile` such as the following will probably work: 

```
locahost,0:0,1,OpenCL,,/tmp
```

The `<num GPUs>` setting is for MD threads running on multiple GPUs, that are not currently supported. The only allowed value is 1.

In this tutorial, the RE simulation of each complex employs 16 replicas. ASyncRE assumes that there are more replicas than computing devices. In this case, it is not recommended to use more than 8 GPUs.

### Step 6: run the workflow

```
cd $HOME/t4l/scripts
bash ./setup-sdm.sh
```

The workflow will first set up the receptor, and then each complex in turn.

At the end, the script will run a minimization and thermalization cycle for each complex using the local GPU. If the setup is performed on a machine without GPUs, comment out the corresponding section and run the minimization thermalization for each complex manually. For example:

```
cd $HOME/t4l/complexes/t4l-toluene
./runopenmm t4l-toluene_mintherm.py
```

## Run the ASyncRE simulations

Go to the simulation directories of each complex and launch the ASyncRE simulations. For example, assuming ASyncRE is installed under `$HOME/devel/async_re-openmm`:

```
cd $HOME/t4l/complexes/t4l-toluene
python $HOME/devel/async_re-openmm/bedamtempt_async_re.py t4l-toluene_asyncre.cntl
```

Each RE simulation is set to run for 8 hours (see `WALL_TIME` above). The amount of samples collected will depend on the number of GPUs utilized. The more samples, the better converged will be the free energy estimate.

Monitor the progress of each ASyncRE simulation by inspecting the `*_stat.txt` file. For example,

```
cd $HOME/t4l/complexes/t4l-toluene
cat t4l-toluene_stat.txt
```

will show something like:

```
Replica  State  Lambda Lambda1 Lambda2 Alpha U0 W0coeff Temperature Status  Cycle 
     0       1   0.057  0.057  0.057 0.000 0.000 0.000 300     W      2 
     1       0   0.000  0.000  0.000 0.000 0.000 0.000 300     W      2 
     2       3   0.171  0.171  0.171 0.000 0.000 0.000 300     W      2 
     3       6   0.343  0.343  0.343 0.000 0.000 0.000 300     W      2 
     4       4   0.229  0.229  0.229 0.000 0.000 0.000 300     W      2 
     5       2   0.114  0.114  0.114 0.000 0.000 0.000 300     W      2 
     6       5   0.286  0.286  0.286 0.000 0.000 0.000 300     W      2 
     7       7   0.400  0.400  0.400 0.000 0.000 0.000 300     W      2 
     8       8   0.457  0.457  0.457 0.000 0.000 0.000 300     R      1 
     9       9   0.514  0.514  0.514 0.000 0.000 0.000 300     W      1 
    10      10   0.571  0.571  0.571 0.000 0.000 0.000 300     W      1 
    11      11   0.629  0.629  0.629 0.000 0.000 0.000 300     W      1 
    12      12   0.686  0.686  0.686 0.000 0.000 0.000 300     W      1 
    13      13   0.743  0.743  0.743 0.000 0.000 0.000 300     W      1 
    14      14   0.800  0.800  0.800 0.000 0.000 0.000 300     W      1 
    15      15   1.000  1.000  1.000 0.000 0.000 0.000 300     W      1 
Running = 1
Waiting = 15
```

In this case it shows that replicas 0 through 7 have completed one cycle and are waiting to perform a second. Replica 8 is currently running. Replicas 9 through 15 are waiting to run their first cycle. Because in this case `PRODUCTION_STEPS=20000` and `PRNT_FREQUENCY=5000`, each cycle is 20,000 MD steps and produces 4 binding energy samples. The `t4l-toluene_stat.txt` also shows the current lambda-state assigned to each replica. We see for example that replica 0 and replica 1 have exchanged states. 

Each replica runs in a separate subdirectory named `r0`, `r1`, etc. Each cycle generates .out, .dms, .log, .pdb, and .dcd files. For example the binding energy samples for the first cycle of replica 2 are stored in the file `r2/t4l-toluene_1.out` (last column). The .dms files are used to start the following cycle. The .pdb and .dcd files are used for trajectory visualization. See below.  

The ASyncRE process can be killed at any time with `^C` and optionally restarted. However, replicas currently running on remote machines are likely to keep running and may have to be killed before ASyncRE can be restarted. To start from scratch (that is from the first cycle) remove the replicas directories by doing `rm -r r? r??`. 

## Analysis

### Cleanup

At any time while the simulations are running, do:

```
cd $HOME/t4l/scripts
bash ./cleanup.sh
```

to clean up the replica directories. The script deletes the .log and .err files and other unnecessary files, except for the last 3 cycles. The script also concatenates the .out and .dcd files. This reduces drastically the number of files and simplifies the inspection of trajectory files.

### Free energy analysis

At any time while the simulations are running, do:

```
cd $HOME/t4l/scripts
bash ./analyze.sh
```

to obtain the binding free energies of the complexes. The first `discard_samples` samples (set in the `setup-settings.sh` file, see above) are discarded. In this tutorial each cycle generates 4 samples.

The script produces an output such as:

```
free energy analysis for ligand benzene
t4l-benzene  DGb = -5.993421 +- 0.3762201 DE = -16.66318 +- 0.1983421  min/max cycles: 4 5
free energy analysis for ligand toluene
t4l-toluene  DGb = -7.617383 +- 0.2985353 DE = -20.26433 +- 0.2746842  min/max cycles: 6 7
```

`DGb` is the binding free energ in kcal/mol. `DE` is the average binding energy in the coupled ensemble (lambda=1).

The output of the R program is in each complex directory. For toluene, for example, it will be in `$HOME/t4l/complexes/t4l-toluene/uwham_analysis.Rout`. Look for errors in this output file if the free energies are not printed. The program also produces plots of the free energy profiles and of the binding energy distributions. For toluene, for example, the plots will be in `$HOME/t4l/complexes/t4l-toluene/Rplots.pdf`.

### Visualization

It is always a good idea to inspect trajectory to see what is going on. After cleanup (see above) the trajectory of a replica can be loaded as follows (for example):

```
cd $HOME/t4l/complexes/t4l-toluene/r4
vmd -f t4l-toluene_6.pdb t4l-toluene.dcd
```

the pdb file is used to define the topology and any present in the replica directory would do. Once in vmd, delete frame 0 to remove it. 

## Appendix

### <a name="msys"></a> A. Installation of `msys`

The `msys` package provides the `mae2dms` utility to convert Maestro files into the DMS format.

We assume an Ubuntu 16.04 system and a Conda 2 environment. Below we assume the Conda environment is located in `$HOME/miniconda2`. Change it to match your environment:

```
source $HOME/miniconda2/bin/activate
```

`msys`'s primary requirements are `boost`, `scons` and `sqlite3`:

```
sudo apt install libboost-dev libboost-all-dev
sudo apt install scons
sudo apt install sqlite3 libsqlite3-dev
```

Retrieve the `msys` sources from the DE Shaw Research github archive. In the commands below we assume this is done in the `devel` folder of the user home directory. Change it to match your build environment. 

```
cd $HOME/devel
git clone https://github.com/DEShawResearch/msys.git
```

Consult the installation section of the README file that comes with `msys` if the steps below do not work for you

```
cd $HOME/devel/msys
export PYTHONPATH=$HOME/devel/msys/external:$PYTHONPATH
scons -j4
scons -j4 PYTHONVER=27 MSYS_BOOST_PYTHON_SUFFIX=-py
```

For this tutorial we assume that `msys` will be installed in `$HOME/local`. Change to match your choice.

```
scons -j4 PYTHONVER=27 MSYS_BOOST_PYTHON_SUFFIX=-py install PREFIX=$HOME/local
```

make sure that `$HOME/local` is in your search paths:

```
export PATH=$PATH:$HOME/local/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/lib
export PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python
```

Test the installation

```
cd $HOME
dms-info $HOME/devel/msys/tests/smarts_tests/PC32699.mae
```
