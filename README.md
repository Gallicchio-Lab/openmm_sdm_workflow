# OpenMM SDM workflow

VERSION 0.2.0

A workflow to setup the Single Decoupling Method (SDM) for alchemical binding free energy calculations in implicit solvent with OpenMM. A similar method on explicit solvent is upcoming.

## Contributors

Emilio Gallicchio egallicchio@brooklyn.cuny.edu

Rajat Pal rajatfor2014@gmail.com

Sheenam Sheenam ssheenam@gradcenter.cuny.edu

Solmaz Azimi sazimi@gradcenter.cuny.edu

## License

This software is released under the LGPL license.

## Credits

This software is maintained by the Gallicchio Laboratory at the Department of Chemistry in Brooklyn College of CUNY ([compmolbiophysbc.org](http://compmolbiophysbc.org)). Development and maintenance of this software is supported in part from a grant from the National Science Foundation ([CAREER 1750511](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1750511&HistoricalAwards=false)).

## Documentation and Tutorials

The documentation is structured in the

[T4 Lysozyme tutorial](tutorials/t4l/tutorial.md),

which prepares SDM calculations for a farm of GPUs, whether local or distributed, on multiple servers.   

## References

For a background on the alchemical process and the conformational sampling techniques used with SDM, consult our publications. Here are some key ones:

* [Binding energy distribution analysis method (BEDAM) for estimation of protein-ligand binding affinities](http://www.compmolbiophysbc.org/publications#bedam_2010)
* [Analytic Model of the Free Energy of Alchemical Molecular Binding](http://www.compmolbiophysbc.org/publications#analytical_theory_2018)
* [Perturbation Potentials to Overcome Order/Disorder Transitions in Alchemical Binding Free Energy Calculations](http://www.compmolbiophysbc.org/project-updates/manuscriptonorderdisordertransitionsinalchemicalbindingfreeenergycalculations)
* [Recent Theoretical and Computational Advances for Modeling Protein-Ligand Binding Affinities](http://www.compmolbiophysbc.org/publications#pubs_advprot_2011)
* [Theory of binless multi-state free energy estimation with applications to protein-ligand binding](http://www.compmolbiophysbc.org/publications#uwham)
* [Asynchronous Replica Exchange Software for Grid and Heterogeneous Computing](http://www.compmolbiophysbc.org/publications#asyncre_software_2015)
* [Efficient Gaussian Density Formulation of Volume and Surface Areas of Macromolecules on Graphical Processing Units](http://www.compmolbiophysbc.org/publications#gaussvol_2017)
* [AGBNP, an analytic implicit solvent model suitable for molecular dynamics simulations and high-resolution modeling](http://www.compmolbiophysbc.org/publications#AGBNP1)


## Installation

The initial system setup is conducted on our [OpenMM/SDM Docker image](http://www.compmolbiophysbc.org/research/research-blog/acentosdockerimageforopenmmsdmdevelopment). It includes all the software necessary for running SDM calculations, including: 

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

Below, a summary of the packages are provided that include instructions on their installation and configuration. In the interest of ease and efficiency, the Docker Image is recommended. 

Note the steps below were tested on a 16.04 Ubuntu system.

### Conda 2

Python bindings for OpenMM and other libraries are best stored in a Conda environment. We recommend the Miniconda 2 environment:

```
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```

In this tutorial we assume that Miniconda 2 is installed under `$HOME/miniconda2`. Adjust as needed.

Activate the conda environment:

```
source ~/miniconda2/bin/activate
```

### OpenMM

We recommend installing OpenMM from [source](https://github.com/pandegroup/openmm) as several necessary steps are similar to that of building SDM-related plugins (see below). Binary installations of OpenMM can also work. Detailed building instructions for OpenMM are [here](http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code). SDM requires an OpenCL platform with GPUs from NVIDIA (CUDA) or AMD, which we assume are in place.

The steps we used to build OpenMM 7.3.1 on an Ubuntu 16.04 system involve:

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

Hit `c` (configure) until all variables are correctly set. Set `CMAKE_INSTALL_PREFIX` to point to the openmm installation directory. Here we assume `$HOME/local/openmm-7.3.1`. If an OpenCL platform is detected (such as from a NVIDIA CUDA installation) it will be enabled automatically. In the context of this tutorial, the CUDA platform is optional. Hit `g` to generate the makefiles and `q` to exit `ccmake`, continue with:

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

The `sqlitebrowser` application is very useful for inspecting DMS files.

### Desmond

In this tutorial, in order to assign OPLS2005 force field parameters, we will use academic Desmond as part of the Schrodinger's environment. Download and install [Maestro/Desmond](https://www.deshawresearch.com/downloads/download_desmond.cgi). Naturally, the commercial version of Maestro/Desmond is also functional.

### Msys

`msys` is a software package developed by DE Shaw Research to manipulate molecular structures. We use it to convert Maestro-formatted files to DMS-formatted files. Installation instructions are in [Appendix A](#msys) below. We assume that `msys` is installed in `$HOME/local` and that executables and libraries there are in the shell search paths. See [Appendix A](#msys).


### SDM Workflow

The package includes:

```
cd $HOME/devel
git clone https://github.com/egallicc/openmm_sdm_workflow.git
```

In order to set up the simulations for this tutorual,  we will use the structures stored under `openmm_sdm_workflow/tutorial` and modify the files under `openmm_sdm_workfklow/tutorial/scripts`. 

### SDM OpenMM plugin

This plugin implements a customized Langevin integrator to compute the alchemical potential in SDM.

```
cd $HOME/devel
git clone https://github.com/rajatkrpal/openmm_sdm_plugin.git
```

Follow the [installation instructions](https://github.com/rajatkrpal/openmm_sdm_plugin/blob/master/README.md) for this package. Point the `CMAKE_INSTALL_PREFIX` and `OPENMM_DIR` to the OpenMM installation directory (`$HOME/local/openmm-7.3.1` in this example).


### AGBNP OpenMM plugin

This tutorial is conducted in the context of an implicit solvent. This plugin implements the AGBNP implicit solvent model in OpenMM. A tutorial with explicit solvent is in progress.

```
cd $HOME/devel
git clone https://github.com/egallicc/openmm_agbnp_plugin.git
```

Follow the [installation instructions](https://github.com/egallicc/openmm_agbnp_plugin/blob/master/README.md) for this package. Point the `CMAKE_INSTALL_PREFIX` and `OPENMM_DIR` to the OpenMM installation directory (`$HOME/local/openmm-7.3.1` in this example).

### ASyncRE for OpenMM

ASyncRE performs replica exchange simulations in asynchronous mode across a wide range of computational devices and grids. The specific version used by SDM distributes OpenMM jobs across GPU compute servers through passwordless `ssh`. For installation, proceed as:

```
conda install numpy configobj paramiko
cd $HOME/devel
git clone https://github.com/egallicc/async_re-openmm.git
cd async_re-openmm
python setup.py install
```

To set up password-less `ssh` across a farm of computers, refer to [guide](https://www.digitalocean.com/community/tutorials/how-to-set-up-ssh-keys-on-ubuntu-1604). ASyncRE does not require that the compute servers share a common user space and a shared filesystem.

### UWHAM for R

The unbinned weighted histogram analysis method ([UWHAM](http://www.compmolbiophysbc.org/publications#uwham)) is an implementation of the MBAR thermodynamic reweighting method for R. In this tutorial, UWHAM computes the binding free energies from the output of SDM calculations. The UWHAM package is part of the CRAN archive. For installation, proceed as:

```
R
install.packages("UWHAM")
```

### VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is an indispensable tool for visualization and analysis of trajectories. We use the `catdcd` tool distributed with VMD to parse and concatenate .dcd trajectory files. For this tutorial we will assume that VMD is installed in `/usr/local`.

## Appendix

### <a name="msys"></a> A. Installation of `msys`

The `msys` package provides the `mae2dms` utility to convert Maestro files into the DMS format.

We assume an Ubuntu 16.04 system and a Conda 2 environment. Below we assume the Conda environment is located in `$HOME/miniconda2`.Modify it to match your environment accordingly:

```
source $HOME/miniconda2/bin/activate
```

The primary requirement for`msys`are `boost`, `scons` and `sqlite3`:

```
sudo apt install libboost-dev libboost-all-dev
sudo apt install scons
sudo apt install sqlite3 libsqlite3-dev
```

Retrieve the `msys` sources from the DE Shaw Research github archive. In the commands below, we assume this is done in the `devel` folder of the user home directory. Modify it accordingly to match your build environment. 

```
cd $HOME/devel
git clone https://github.com/DEShawResearch/msys.git
```

Should the steps below cause issues, consult the installation section of the README file that comes with `msys`.

```
cd $HOME/devel/msys
export PYTHONPATH=$HOME/devel/msys/external:$PYTHONPATH
scons -j4
scons -j4 PYTHONVER=27 MSYS_BOOST_PYTHON_SUFFIX=-py
```

For this tutorial we assume that `msys` will be installed in `$HOME/local`. Modify it accordingly to match another set-up.

```
scons -j4 PYTHONVER=27 MSYS_BOOST_PYTHON_SUFFIX=-py install PREFIX=$HOME/local
```

Ensure that `$HOME/local` is in your search paths:

```
export PATH=$PATH:$HOME/local/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/lib
export PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python
```

Test the installation:

```
cd $HOME
dms-info $HOME/devel/msys/tests/smarts_tests/PC32699.mae
```
