# NetworkGen

A MATLAB-based mesoscale polymer network generator that produces simulation-ready input files for [LAMMPS](https://www.lammps.org/). NetworkGen gives researchers fine-grained control over network topology, strand length distributions, and defect structures.

**Documentation:** [soft-matter-lab.github.io/networkgen](https://soft-matter-lab.github.io/networkgen)

---

## Features

- Random and hexagonal lattice node placement
- Monodisperse, polydisperse, and bimodal strand length distributions
- Double network generation
- Void defects with controllable size, shape, and spatial distribution
- Geometric and topological lattice disorder
- Multi-type atom and bond support
- Local density potential table generation
- Structured log files with network statistics
- Config Builder web tool for generating scripts without writing code
- Python support via [Octave](https://octave.org/) and [oct2py](https://github.com/blink1073/oct2py)

---

## Getting Started

### MATLAB

**Requirements:** MATLAB R2020a or later.

**Option 1 — Download a release**

Download the latest release from the [Releases page](https://github.com/soft-matter-lab/networkgen/releases). On Windows a one-click installer is available. Otherwise unzip and add to your MATLAB path:

```matlab
addpath(genpath('/path/to/NetworkGen'));
```

**Option 2 — Clone**

```bash
git clone https://github.com/soft-matter-lab/networkgen.git
```

```matlab
addpath(genpath('/path/to/networkgen/NetworkGen'));
```

**Verify:**

```matlab
net = network();
disp(net)
```

### Python (via Octave)

**Requirements:** Octave 7.0+, Python 3.10+, `oct2py`, `numpy`

```python
from oct2py import octave

octave.addpath(octave.genpath('/path/to/NetworkGen'))

octave.eval("""
    net = network();
    net.domain.Lx = 10;
    net.domain.Ly = 10;
    net.flags.iplot = false;
    net.generateNetwork();
""")
```

---

## Usage

Create a `network` object, configure its properties, and call `generateNetwork()`:

```matlab
net = network();

% Domain
net.domain.Lx = 150;
net.domain.Ly = 150;
net.domain.boundary = 'fixed';
net.domain.write_location = './output';

% Architecture
net.architecture.geometry = 'random';
net.architecture.strand_typology.mode = 'bimodal';
net.architecture.strand_typology.bimodal.mean_1 = 10;
net.architecture.strand_typology.bimodal.mean_2 = 40;

% Flags
net.flags.isave = true;
net.flags.iplot = false;

net.generateNetwork();
```

Use the **[Config Builder](https://soft-matter-lab.github.io/networkgen/config-builder)** on the documentation site to generate MATLAB or Python scripts interactively without writing code.

---

## Output Files

| File | Description |
|------|-------------|
| `.dat` | LAMMPS data file — atom positions, bond connectivity, type assignments |
| `.lammps` | LAMMPS visualization file for OVITO or VMD |
| `bond.table` | Bond potential lookup table |
| `console_*.log` | Command window transcript |
| `network_*.log` | Network statistics and generation parameters |

---

## Installing & Running LAMMPS

NetworkGen output files are designed for use with the custom [meso-network](https://github.com/SoftMatterLab24/lammps/tree/meso-network) branch of LAMMPS. NetworkGen itself does not require LAMMPS — it is only needed to run the generated simulation files.

While LAMMPS offers both Windows and Mac compatibility, we recommend using a Linux environment. For Windows users the recommended option is to install Windows Subsystem for Linux (WSL). A useful guide can be found [here](https://docs.lammps.org/Howto_wsl.html). For a Linux distribution try Ubuntu 20.04.6 LTS.

### Installation Guide

**Clean installation**

1. Open a WSL terminal and navigate to the home directory or location you wish to install LAMMPS
```bash
cd ~
```

2. Clone the repository
```bash
git clone https://github.com/SoftMatterLab24/lammps.git
```

3. Switch to the meso-network branch
```bash
cd lammps
git checkout meso-network
```

4. Make a build directory
```bash
mkdir build
cd build
```

5. Run cmake with required packages
```bash
cmake -D PKG_BPM=yes -D PKG_EXTRA-FIX=yes -D PKG_GRANULAR=yes -D PKG_MISC=yes \
      -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_TNT=yes \
      -D PKG_EXTRA-MOLECULE=yes -D PKG_COLLOID=yes -D PKG_MC=yes \
      -D PKG_BROWNIAN=yes -D PKG_MANYBODY=yes -D PKG_POLY-NET=yes ../cmake
```

6. Compile (replace `4` with your processor count)
```bash
make -j 4
```

**Adding a package after LAMMPS is built**

```bash
cmake -D PKG_<NAME>=on .
make -j 4
```

### Running LAMMPS

**Serial:**
```bash
/path/to/lammps/build/lmp -in file.in
```

**Parallel:**
```bash
mpirun -np 4 /path/to/lammps/build/lmp -in file.in
```

> **Note:** Do not run LAMMPS from within the build folder.

### Troubleshooting

Install prerequisite packages first:

```bash
sudo apt install -y cmake build-essential ccache gfortran openmpi-bin libopenmpi-dev \
                    libfftw3-dev libjpeg-dev libpng-dev python3-dev python3-pip \
                    python3-virtualenv libblas-dev liblapack-dev libhdf5-serial-dev \
                    hdf5-tools
```

**C++ compiler errors**

If the build fails with a compiler error, check that CMake found the right compiler:

```bash
which c++
```

If the output is blank, install build tools:

```bash
sudo apt install build-essential
```

Or set the compiler manually in the cmake step:

```bash
-D CMAKE_CXX_COMPILER=c++
```

---

## License

See [LICENSE](LICENSE) for details.
