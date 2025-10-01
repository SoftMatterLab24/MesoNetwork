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
