# Mesoscale Network Models

This repository contains various mesoscale network models for studying mechanical properties of elastomeric materials.

All projects utilize the opensource [LAMMPS](https://www.lammps.org/#gsc.tab=0) MD codebase, under the custom [meso-network](https://github.com/SoftMatterLab24/lammps/tree/meso-network) branch. Please see the Installing LAMMPS section for more guidance.

## Branches (OUTDATED)
Code for specific models can be found within the following branches

- **Poly** -- Elastic networks with polydispersity
- **Visco** -- Viscoelastic networks with bimodal chain length distribution

## Using the Network Generator
The network generator is a multi-purpose tool that can be used to generate many different types of networks (hexagonal lattice, random, double) whose properties (i.e. contour length, connectivity, etc.) seek to follow physically realistic statistics. The network generator requires MATLAB Version 2016b or later. There are two ways to use the Generator: 1. (Open and run the ```/NetworkGen/NetworkGenDriver.m``` directly) OR 2. (Create a "network" object from the ```/NetworkGen/network.m``` class).

### 1. Running the Generator from the Driver
To execute the generator from the driver, first download the source code either manually or by cloning the repository. Once installed open ```/NetworkGen/NetworkGenDriver.m```

...to be expanded upon

### 2. Generating a network from a network object
An alternative approach is to create a network object, modify the desired properties, and then pass the object to the ```generateNetwork``` function. This may be beneficial if many networks each with different properties need to be generated, which could otherwise become cumbersome to make with the driver. To use this method:
1. Open a Linux terminal and navigate to the home directory or location you wish to install the generator. Then clone this repository:
```
git clone https://github.com/SoftMatterLab24/MesoNetwork.git
```
2. Once installed open a powershell terminal and in the root of the cloned directory run
```powershell
.\addPath.bat
```
This will save the network class to your MATLAB PATH, so that it can be accessed from any open editor. A defender window will prompt you to allow MATLAB to make these changes. Once this is done exit powershell.

3. Now in a MATLAB script you can create a network object from the network class, modify its properties, and generate a network as:
```matlab
clear all

% Create a network object
n = network;  

% Modify network properties
n.Lx = 100;               % Change the x-dim boundary length
n.dist_type = 'bimodal';  % Change the Kuhn segment distribution type

% Generate a network
[Domain, Atoms, Bonds, Nvec, order] = generateNetwork(n);
```

### 3. Using the networks with LAMMPS
Unless saving is disabled the generator will automatically create the required data and table files needed to run LAMMPS simulations. By default four files are saved:
1. data file (.dat) which contains the atoms and bonds data to be used by the [read_data](https://docs.lammps.org/read_data.html) command.
2. bond table file (bond.table) which stores extra bond data to be used by various bond styles in the [meso-network](https://github.com/SoftMatterLab24/lammps/tree/meso-network) branch
3. local density table file (.localdensity.table) which stores tabulated data for the [local_density](https://docs.lammps.org/pair_local_density.html) pairstyle
4. A log file (.log) with summary information about the generated network, which is also used by Simulation Builder.

## Using the Simulation Builder

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
cmake -D PKG_BPM=yes -D PKG_EXTRA-FIX=yes -D PKG_GRANULAR=yes -D PKG_MISC=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_TNT=yes -D PKG_EXTRA-MOLECULE=yes -D PKG_COLLOID=yes -D PKG_MC=yes -D PKG_BROWNIAN=yes -D PKG_MANYBODY=yes ../cmake
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
