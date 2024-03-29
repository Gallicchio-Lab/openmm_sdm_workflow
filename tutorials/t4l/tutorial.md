# OpenMM SDM Tutorial for a set of T4 Lysozyme complexes

1. [Introduction](#introduction)  
2. [Gathering the Software Tools](#gathering-the-software-tools)
3. [Setup the Simulations](#setup-the-simulations)
4. [Run the ASyncRE simulations](#run-the-asyncre-simulations)
5. [Analysis](#analysis)

## Introduction

This tutorial is a step-by-step guide to setup and run alchemical absolute binding free energy calculations for a series of complexes of a T4 lysozyme mutant using the OpenMM's Single Decoupling Method (SDM) plugin. 

SDM is based on an alchemical process in which the ligand is progressively transferred from the solution environment to the receptor binding site. SDM employs an implicit description of the solvent (the AGBNP model used here), which allows it to avoid the intermediate vacuum state necessary for absolute binding free energy calculations with explicit solvation, such as in the Double Decoupling Method (DDM). Hence, SDM requires only one free energy calculation, as opposed to two with DDM. Hamiltonian Parallel Replica Exchange with OpenMM is used for conformational sampling.

For a background on the alchemical process and the conformational sampling techniques used with SDM, consult our publications. Here are the key ones:

* [Binding energy distribution analysis method (BEDAM) for estimation of protein-ligand binding affinities](http://www.compmolbiophysbc.org/publications#bedam_2010)
* [Analytic Model of the Free Energy of Alchemical Molecular Binding](http://www.compmolbiophysbc.org/publications#analytical_theory_2018)
* [Perturbation Potentials to Overcome Order/Disorder Transitions in Alchemical Binding Free Energy Calculations](http://www.compmolbiophysbc.org/project-updates/manuscriptonorderdisordertransitionsinalchemicalbindingfreeenergycalculations)
* [Recent Theoretical and Computational Advances for Modeling Protein-Ligand Binding Affinities](http://www.compmolbiophysbc.org/publications#pubs_advprot_2011)
* [Theory of binless multi-state free energy estimation with applications to protein-ligand binding](http://www.compmolbiophysbc.org/publications#uwham)
* [Asynchronous Replica Exchange Software for Grid and Heterogeneous Computing](http://www.compmolbiophysbc.org/publications#asyncre_software_2015)
* [Efficient Gaussian Density Formulation of Volume and Surface Areas of Macromolecules on Graphical Processing Units](http://www.compmolbiophysbc.org/publications#gaussvol_2017)
* [AGBNP, an analytic implicit solvent model suitable for molecular dynamics simulations and high-resolution modeling](http://www.compmolbiophysbc.org/publications#AGBNP1)


### Molecular Systems

In this tutorial, we consider the complexes between benzene, toluene, and 3-iodotoluene with the L99A mutant of the T4 lysozyme receptor (PDB id: 4W53). The systems were prepared using the free academic version of [Maestro/Desmond](https://www.deshawresearch.com/downloads/download_desmond.cgi). Toluene, water molecules and other bound ligands were removed from the 4W53 receptor structure. To reduce the size of the system and speed-up the calculations, residues 1 through 71 on the receptor were removed. The receptor was then processed with Protein Preparation tool in Maestro to add hydrogen atoms and cap the termini. OpenMM prefers protein structures in which atoms belonging to the same residue are listed consecutively. To do so the processed receptor structure was saved in PDB format and read back into Maestro.

Benzene and 3-iodotoluene were built with academic Maestro using bound toluene as a template. Bond order and formal charges were adjusted and hydrogen atoms were added in Maestro. In most actual applications, the starting bound conformations of ligands are obtained using molecular docking, using Glide for example. The commercial version of the Schrodinger's Suite includes an easy-to-follow tutorial on the steps to prepare protein receptors and ligands for docking. The end result is the receptor structure and a list of prepared and docked ligand structures in the Maestro project table as in the present case.  

The receptor and the ligands were exported from Maestro in separate folders (`ligands/` and `receptor/`) in `.maegz` V1 formatted files (see below).   

## Setup on Docker

Although the system can be setup manually, all the software tools are gathered in Docker to expedite setup. First, download and install Docker accordingly. Then proceed with obtaining our Docker container with:

```
docker pull egallicchio/centos610-sdmsetup:latest
```
Once the pull has been completed, run the container:

```
docker run -it --rm egallicchio/centos610-sdmsetup:latest
```

Copy the SDM directory to Docker: 

```
docker cp openmm_sdm_workflow [docker id]:/tmp/
```

The files can be modified before or after the copy, however the setup must be run in `$HOME` on Docker. Files can be read in `$HOME` but not written. The `tmp/` directory bypasses this complication.


### Step 1: set the simulation parameters in `setup-settings.sh`

Edit the `setup-settings.sh` file in the `scripts` folder. For each setting, replace every placeholder designated as `< >` with a specific input using the instructions below.
For example, change the receptor name from:
```
receptor=<basename of the .maegz receptor file>
```
to:

```
receptor=t4l99a
```

Settings:

* `set_basename`: the name of the project, that labels files, simulation directories, etc. For this tutorial set it to `t4l`.
* `main_dir`: the working directory for this project. For this tutorial set it to `$HOME/openmm_sdm_workflow`.
* `work_dir`: where the `ligands` and `receptor` directories are stored. For this tutorial set it to `${main_dir}/tutorials/${set_basename}`.
* `scripts_dir`: where the scripts are stored. For this tutorial set it to `${main_dir}/scripts`.
* `schrodinger`: your Schrodinger installation directory. For example `/opt/Desmond_Maestro_2018.4`
* `msys_path`: the msys installation directory. If you followed the instructions above it will be `$HOME/local`.
* `vmd_path`: the vmd installation directory, usually `/usr/local`.
* `receptor`: the basename of the receptor file. `t4l99a` for this tutorial. The setup script will then look for `t4l99a.maegz` in the `receptor` folder.
* `ligands`: a string with the list of basenames of the ligands. Set it to `benzene toluene 3iodotoluene` here.
* `cntltmpl`: the name of the template control file for the SDM workflow (see below). Leave this unchanged to `sdm_workflow_template.cntl`
* `agbnpparam`: the name of the AGBNP parameter file for the solvent contex. Leave this unchanged to `agbnp2.param.agbnp_plugin`.
* `rest_receptor_sql`: a string with a SQLite atom selection specifying which atoms of the receptor will be restrained. SDM uses isotropic, flat-bottom harmonic potentials. For example `name GLOB 'CA'` is used to restrain C-alpha atoms. Notice that the single quotes need to be escaped in the shell. For this tutorial set `rest_receptor_sql` to `$'name GLOB \'CA\''`
* `rest_receptor_kf`: the force constant of the receptor restraints in kcal/mol/angstrom^2. Set it to `25.0`.
* `rest_receptor_tol`: the tolerance of the receptor restraints in angstrom. Set it to `0.75`.
* `discard_samples`: how many samples to discard from the beginning of the trajectory for equilibration. Set it to 0. (For any other setting, make sure to obtain more than the number of samples discarded from each replica). 

### Step 2: set the alchemical schedule and ASyncRE settings in the `sdm_workflow_template.cntl` file

Edit the `sdm_workflow_template.cntl` file in the `scripts` folder. For each setting, replace every placeholder designated as `< >` with a specific input using the instructions below.
For example, chang ethe `WALL_TIME` from:

```
WALL_TIME = <duration of each ASyncRE simulation, in minutes>
```

to:

```
WALL_TIME = 480
```

Settings:

* `JOB_TRANSPORT`: job transport mechanism. Set this setting to 'LOCAL_OPENMM', which enables ASyncRE communication via memory. See below, for more information about this setting.
* `TEMPERATURES`: list of replica exchange temperatures. Set this to the single temperature `300`.
* `LAMBDAS`: list of alchemical lambda values in comma-separated string. Set it to `' 0.000, 0.057, 0.114, 0.171, 0.229, 0.286, 0.343, 0.400, 0.457, 0.514, 0.571, 0.629, 0.686, 0.743, 0.800, 1.000'`
* `CYCLE_TIME`: period of RE exchanges in seconds. Set it to `30`.
* `WALL_TIME`: wall-clock duration of the RE simulation for each complex in minutes. Set it to `480`.
* `PRODUCTION_STEPS`: number of MD steps for each running cycle of a replica before checkpointing. Set it to `1000`.
* `PRNT_FREQUENCY`: the MD period, in MD steps, for saving binding energy samples. Set it to `5000`.
* `TRJ_FREQUENCY`: the MD period, in MD steps, for saving trajectory frames. Set it to `5000`.
* `REST_LIGAND_CMLIGSQL`: sqlite selection specifying the CM atoms of the ligand. Set it to `'name GLOB 'C?''` to use the CM of the carbon atoms of the ligands.
* `REST_LIGAND_CMRECSQL`: same as above but for the receptor. Set it to `'name GLOB 'CA' AND resid IN (79,84,88,91,96,104,112,113,122,133,150)'` to use the CM determined by the C-alpha atoms of selected residues.
* `REST_LIGAND_CMKF`: force constant of flat bottom CM-CM restraint, in kcal/mol/angstrom^2. Set it to `25.0`.
* `REST_LIGAND_CMTOL`: tolerance of flat bottom CM-CM restraint, in angstrom. Set it to `2.5`.
* `SOFT_CORE_UMAX`: max binding energy of soft-core potential in kcal/mol. Set it to `50.0`.
* `SOFT_CORE_ACORE`: exponent of rational soft core potential, dimensionless. Set it to `0.0625` or 1/16.


### Step 3: configure the `runopenmm` launch script

Edit the `runopenmm` launch script to reflect your environment. For example, if OpenMM is installed in `$HOME/local/openmm-7.3.1` and the corresponding python bindings are installed in the Conda environment under `$HOME/miniconda2`, the `runopenmm` should look like:

```
#!/bin/bash
openmm_dir=$HOME/local/openmm-7.3.1
pythondir=$HOME/miniconda2
export LD_LIBRARY_PATH=${openmm_dir}/lib:${openmm_dir}/lib/plugins:$LD_LIBRARY_PATH
${pythondir}/bin/python "$@"
```

Technically, the `runopenmm` launch script should reflect the OpenMM installation environment on the GPU compute machines (see below). In this tutorial we assume that the local machine and the remote compute machines (which may include the local machine) execute OpenMM similarly.

### Step 4: prepare the `nodefile`

The ASyncRE system dispatches jobs to GPU computing devices on the machines listed in the `nodefile`. The format of each line is as follows:

```
<machine name> , <OpenCL platform id>:<OpenCL GPU device id> , <num GPUs> , <platform name> , <remote username>, <remote temp directory>
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

The `<num GPUs>` setting is for MD threads running on multiple GPUs, which is a setting that is not currently supported. The only allowed value is 1.

In this tutorial, the RE simulation of each complex employs 16 replicas. ASyncRE assumes that there are more replicas than computing devices. In this case, it is not recommended to use more than 8 GPUs.

### Step 5: run the setup workflow on Docker

Copy the directory from `tmp/` to `$HOME` on Docker. 

```
cd $HOME
cp /tmp/openmm_sdm_workflow . 
cd /openmm_sdm_workflow/scripts 
bash ./setup-sdm.sh
```

The workflow will first set up the receptor, and then each complex in turn.

Minimization and thermalization must be conducted in separately for each complex prior to conducting the simulation.  

### Step 6: minimization and thermalization 

Minimization and thermalization can be done in the local computer, as well as on a supercomputer with the appropriate procedure. This tutorial will assume local implementation. For the complex of t4l and benzene:

```
./runopenmm t4l-benzene_mintherm.py
```

This step will produce the t4l-benzene_lig_0.dms and t4l-benzene_rcpt_0.dms files needed to initiate ASyncRE.


## Run the ASyncRE simulations

Go to the simulation directories of each complex and launch the ASyncRE simulations. For example, assuming ASyncRE is installed under `$HOME/devel/async_re-openmm`:

```
cd $HOME/t4l/complexes/t4l-benzene
./runopenmm $HOME/devel/async_re-openmm/bedamtempt_async_re.py t4l-benzene_asyncre.cntl
```

Each RE simulation is set to run for 8 hours (corresponding to 480 minutes as indicated in `WALL_TIME` above). The amount of samples collected will depend on the number of GPUs utilized, which will also allow for more checkpointing. The more samples that are collected in a simulation, the better converged the binding free energy will be.

Monitor the progress of each ASyncRE simulation by inspecting the `[complex]_stat.txt` file. For example,

```
cd $HOME/t4l/complexes/t4l-benzene
cat t4l-benzene_stat.txt
```

will demonstrate information similar to below:

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

In this case it shows that replicas 0 through 7 have completed one cycle and are waiting to perform a second. Replica 8 is currently running. Replicas 9 through 15 are waiting to run their first cycle. In this case `PRODUCTION_STEPS=1000` and `PRNT_FREQUENCY=5000`, meaning that 1/5 of the total cycles is the total number of samples that are collected. For example, if the simulation produces 50 cycles, at most 10 samples are collected in total. The `t4l-benzene_stat.txt` also shows the current lambda-state assigned to each replica. We see for example that replica 0 and replica 1 have exchanged states. 

Each replica runs in a separate subdirectory named `r0`, `r1`, etc. Each cycle generates `*.out`, `*.dms`, `*.log`, `*.pdb`, and `*.dcd files`. For example the binding energy samples for the first cycle of replica 2 are stored in the file `r2/t4l-benzene_1.out` (last column). The `*.dms` files are used to begin the next cycle. The `*.pdb` and `*.dcd` files are used for trajectory visualization. See below.  

The ASyncRE process can be killed at any time with `^C` and optionally restarted. However, replicas currently running on remote machines are likely to keep running and may have to be killed before ASyncRE can be restarted. To start from scratch (that is from the first cycle) remove the replicas directories by the command `rm -r r? r??`. 

## Job Transport Mechanism: LOCAL_OPENMM vs. SSH

ASyncRE supports two mechanisms to run replicas on GPU devices. The first is the file-based SSH job transport mechanism that is suitable for a farm of GPU computing servers accessible through ssh. This transport mode implements files of OpenMM directly and enables ASyncRE children to communicate via folder. 

The `LOCAL_OPENMM` transport is more suitable for a set of local GPUs residing on the machine running the ASyncRE manager. In this case, the LOCAL_OPENMM transport yields faster job turn-around and more frequent replica exchanges. 

To run this tutorial examples with the SSH transport:

1. In the `sdm_workflow_template.cntl` file set:

* `JOB_TRANSPORT` to `SSH` 
* `PRODUCTION_STEPS` to `20000`.
* `PRNT_FREQUENCY` to `5000`.
* `TRJ_FREQUENCY` to `5000`.

to perform (asynchronous) exchanges every `20000` steps and record data every `5000` steps. As a result, this setting will collect 4 samples from the cycles that are produced.

2. Set the `nodefile` to point to a set of GPUs on the local machine. For example:

```
locahost,0:0,1,OpenCL,,/tmp
locahost,0:1,1,OpenCL,,/tmp
locahost,0:2,1,OpenCL,,/tmp
locahost,0:3,1,OpenCL,,/tmp
```

for a compute server with 4 GPUs when the OpenCL platform is the first platform (0).

The 'LOCAL_OPENMM' transport mechanism stores data in compacted form as the simulation proceeds and does not require the cleanup step before analysis (see below).

## Analysis

### Cleanup

If running with the `SSH` job transport mechanism, at any time while the simulations are running, do:

```
cd $HOME/t4l/scripts
bash ./cleanup.sh
```

to clean up the replica directories. The script deletes the .log and .err files and other unnecessary files, with the exception of the last 3 cycles. The script also concatenates the `*.out` and `*.dcd files`. This significantly reduces the number of files and simplifies the inspection of trajectory files. This step is not necessary when using the `LOCAL_OPENMM` job transport.

### Free energy analysis

The free energy estimate can be obtained at any time during the simulations, with:

```
cd $HOME/t4l/scripts
bash ./analyze.sh
```

to obtain the binding free energies of the complexes. The first `discard_samples` samples (set in the `setup-settings.sh` file, see above) are discarded. In this tutorial each cycle generates 1/5 samples.

The script produces an output such as:

```
free energy analysis for ligand benzene
t4l-benzene  DGb = 0.1430045 +- 0.4394011 DE = -10.66937 +- 0.4186568  min/max cycles: 46 47
```

`DGb` is the binding free energ in kcal/mol. `DE` is the average binding energy in the coupled ensemble (lambda=1).

The output of the R program is in each complex directory. For benzene, for example, it will be in `$HOME/t4l/complexes/t4l-benzene/uwham_analysis.Rout`. Inspect for errors in this output file if the free energies are not printed. The program also produces plots of the free energy profiles and of the binding energy distributions. For benzene, for example, the plots will be in `$HOME/t4l/complexes/t4l-benzene/Rplots.pdf`.

### Visualization

Best practice involves inspecting MD trajectories to visualize the simulation. After cleanup (see above), the trajectory of a replica can be loaded as follows:

```
cd $HOME/t4l/complexes/t4l-benzene/r4
vmd -f t4l-benzene.pdb t4l-toluene.dcd
```

The pdb file is used to define the topology and it is conventionally found in the main complex directory. On VMD, remove frame 0. 

Note that the alchemical lambda value changes with time. Do not be surprised if the ligand overlaps with the receptor during the alchemical trajectory. This is expected when lambda is small. The trajectory will hopefully show multiple binding and unbinding events. These are an indication of optimal convergence.
