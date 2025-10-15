# Mesoscale Network Models

This repository contains various mesoscale network models for studying mechanical properties of elastomeric materials.

All projects utilize the opensource [LAMMPS](https://www.lammps.org/#gsc.tab=0) MD codebase, under the custom [meso-network](https://github.com/SoftMatterLab24/lammps/tree/meso-network) branch. Please see the Installing LAMMPS section for more guidance.

## Branches
Code for specific models can be found within the following branches

- **Poly** -- Elastic networks with polydispersity
- **Visco** -- Viscoelastic networks with bimodal chain length distribution

## Installing & Running LAMMPS
While LAMMPS offers both Windows and Mac compatability, we recommend using a Linux environment. For Windows users the recommended option is to install Windows Subsystem for Linux (WSL). A usefull guide can be found [here](https://docs.lammps.org/Howto_wsl.html). For a Linux distribution try Ubuntu 20.04.6 LTS.

### Installation Guide
**Clean installation**
1. Open a WSL terminal and navigate to the home directory or location you wish to install LAMMPS
```
cd ~
```
2. Clone the repository (SoftMatterLab24)
```
git clone https://github.com/SoftMatterLab24/lammps.git
```
3. Switch to the meso-network branch
```
cd lammps
git checkout meso-network
```
4. Make a build directory
```
mkdir build
cd build
```
5. Run cmake with proper packages (more packages can be apppended as necessary)
```
cmake -D PKG_BPM=yes -D PKG_EXTRA-FIX=yes -D PKG_GRANULAR=yes -D PKG_MISC=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_TNT=yes -D PKG_EXTRA-MOLECULE=yes -D PKG_COLLOID=yes -D PKG_MC=yes -D PKG_BROWNIAN=yes ../cmake
```
6. Compile lammps (enter the number of processors you have, if you don’t know choose 4)
```
make -j 4
```

**To add another package after LAMMPS is built**
1. Enter your lammps build folder
2. Append the new package (you must have the period at the end!)
```
cmake -D PKG_<NAME>=on .
```
3. Rebuild LAMMPS
```
make -j 4
```

### Running LAMMPS
For a more in-depth guide please see [Basics of running LAMMPS](https://docs.lammps.org/Run_basics.html). LAMMPS is run from the command-line, reading an input script via the ```-in``` command-line flag. The input script may use any name or file extention, however, the ```file.in``` format is preferred. LAMMPS is normally run in the directory where your input script is located. This is also where output files are produced by default.

***DO NOT*** run lammps within the build folder, as this can significatly clog your LAMMPS installation. Additionally, data may be lost when attempting to freshly reinstall LAMMPS. It is therefore recommended to run LAMMPS within an entirely seperate directory from the ```lammps``` installation.

**Serial Execution**
To run LAMMPS in serial, (i.e. one processor):
```
/path/to/lammps/src/build/lmp -in file.in
```
where ```/path/to/lammps/``` is the directory you installed LAMMPS, typically ```/home/usr/lammps/```.

**Parallel Execution**
To run LAMMPS in parallel, if it was built with more than 1 processor:
```
mpirun -np 4 /path/to/lammps/src/build/lmp -in file.in
```
where the ```-np``` command-line flag specifies the number of processors to use.

### Troubleshooting
Below is a list of commonly encountered issues when trying to install lammps. Before trying one of the specific solutions below, please ensure that the requisite packages are installed:
```
sudo apt install -y cmake build-essential ccache gfortran openmpi-bin libopenmpi-dev \
                    libfftw3-dev libjpeg-dev libpng-dev python3-dev python3-pip \
                    python3-virtualenv libblas-dev liblapack-dev libhdf5-serial-dev \
                    hdf5-tools
```

**Selecting the correct C++ Compiler**

If LAMMPS fails to build following ``` make -j 8 ``` giving an error about the compiler, CMAKE may not be using the correct compiler when attempting to build LAMMPS. After running ```cmake``` (Step 5) check that ```C++ Compiler:``` is set to ``` /usr/bin/c++ ``` in the build configuration, although the exact location may differ. If this is not the case:
1. Check to see if the C++ compiler is installed by running
```
which c++
```
This should return the directory location of the compiler, i.e. ``` /usr/bin/c++ ```. If it returns blank, this means C++ is either not installed or cannot be found.

2. Install the ```build-essentials``` package, which includes various compilers and development tools:
```
sudo apt install build-essential
```
Then check to see if LAMMPS compiles. Make sure to rebuild LAMMPS with ```cmake```. If this does not work then: 

3. Manually set the compiler by including the following CMAKE build option
```
-D CMAKE_CXX_COMPILER=c++
```
This should set CMAKE to use the correct compiler during configuration. 
