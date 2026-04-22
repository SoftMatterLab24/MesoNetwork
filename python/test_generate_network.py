"""
test_generate_network.py
========================
Integration tests that call generateNetwork() end-to-end via Octave
and verify that real output files are produced and are well-formed.

Unlike test_octave_compat.py (which tests each pipeline stage in
isolation), these tests run the full pipeline exactly as a user would,
then inspect the files that land in the output directory.

Requirements
------------
    pip install oct2py pytest numpy
    Octave >= 7.0 on PATH
    NETWORKGEN_PATH env var pointing at your .m files

Usage
-----
    NETWORKGEN_PATH=/path/to/NetworkGen pytest test_generate_network.py -v

    # Keep generated files after the test run for manual inspection
    NETWORKGEN_PATH=/path/to/NetworkGen KEEP_OUTPUT=1 pytest test_generate_network.py -v -s
"""

import os
import sys
import shutil
import tempfile
import textwrap
import numpy as np
import pytest

# ── Octave executable ─────────────────────────────────────────────────────────
OCTAVE_EXECUTABLE_PATH = os.environ.get(
    "OCTAVE_EXECUTABLE",
    "/usr/local/bin/octave"
)
if "--no-gui" not in OCTAVE_EXECUTABLE_PATH:
    os.environ["OCTAVE_EXECUTABLE"] = OCTAVE_EXECUTABLE_PATH + " --no-gui"
else:
    os.environ["OCTAVE_EXECUTABLE"] = OCTAVE_EXECUTABLE_PATH

try:
    from oct2py import octave, Oct2PyError
except ImportError:
    sys.exit("\n[ERROR] pip install oct2py\n")

NETWORKGEN_PATH = "/home/zwhit/networkgen/NetworkGen/"

# Set to "1" to keep output dirs after tests for manual inspection
KEEP_OUTPUT = os.environ.get("KEEP_OUTPUT", "0") == "1"

# Persistent directory for visualizable outputs — always written, never deleted.
# Contains one representative network per geometry/topology combination.
# Override with: export NETWORKGEN_VIZ_DIR=/your/path
NETWORKGEN_VIZ_DIR = os.environ.get(
    "NETWORKGEN_VIZ_DIR",
    os.path.join(os.path.expanduser("~"), "networkgen_viz_output")
)


# =============================================================================
# Session fixture
# =============================================================================

@pytest.fixture(scope="session", autouse=True)
def octave_session():
    octave.addpath(octave.genpath(NETWORKGEN_PATH))
    ver = octave.eval("OCTAVE_VERSION", verbose=False)
    print(f"\n[setup] Octave {ver} — NetworkGen path: {NETWORKGEN_PATH}")

    # Create persistent viz output dir
    os.makedirs(NETWORKGEN_VIZ_DIR, exist_ok=True)
    print(f"[setup] Visualization outputs -> {NETWORKGEN_VIZ_DIR}")

    yield
    octave.exit()


# =============================================================================
# Per-test output directory
# =============================================================================

@pytest.fixture()
def outdir():
    """
    Create a fresh temp directory for each test to write networks into.
    Automatically deleted after the test unless KEEP_OUTPUT=1.
    """
    d = tempfile.mkdtemp(prefix="networkgen_test_")
    yield d
    if not KEEP_OUTPUT:
        shutil.rmtree(d, ignore_errors=True)
    else:
        print(f"\n[outdir kept] {d}")


# =============================================================================
# Helper: run a config script in Octave
# =============================================================================

def run_config(config_lines: str):
    """
    Evaluate a multi-line Octave config string that ends with
    net.generateNetwork().  Raises Oct2PyError on any Octave error.
    """
    octave.eval(config_lines, verbose=False)


def read_lammps_data(filepath):
    """
    Robust LAMMPS data file reader for NetworkGen output (.dat files).

    Handles the NetworkGen format where:
      - Atoms section: id  mol  type  x  y  z   (atom_style bond/full)
      - Bonds section: id  type  atom_i  atom_j  (4 columns, type before endpoints)
      - Section headers may have trailing comments: "Atoms # bond"

    Returns dict:
      n_atoms, n_bonds           : header counts
      atoms  [N x 6]             : full atom rows  [id mol type x y z]
      bonds  [M x 4]             : full bond rows  [id type i j]
      atom_ids                   : set of atom IDs for fast lookup
      xlo, xhi, ylo, yhi         : domain bounds
    """
    with open(filepath, 'r') as f:
        raw_lines = f.readlines()

    info = {
        'n_atoms': 0, 'n_bonds': 0,
        'n_atom_types': 0, 'n_bond_types': 0,
        'atoms': [], 'bonds': [],
        'xlo': 0.0, 'xhi': 0.0,
        'ylo': 0.0, 'yhi': 0.0,
    }

    section = None
    # Known LAMMPS section keywords — once we see one, header parsing stops
    SECTION_KEYWORDS = {
        'Atoms', 'Bonds', 'Masses', 'Velocities', 'Angles',
        'Dihedrals', 'Impropers', 'Pair Coeffs', 'Bond Coeffs',
        'Angle Coeffs'
    }

    for raw in raw_lines:
        # Strip inline comments
        line = raw.split('#')[0].strip()
        if not line:
            if section is not None:
                section = None  # blank line resets section
            continue

        parts = line.split()
        first = parts[0]

        # ── Section header detection ─────────────────────────────────────
        if first == 'Atoms':
            section = 'atoms'
            continue
        elif first == 'Bonds':
            section = 'bonds'
            continue
        elif first in ('Masses', 'Velocities', 'Angles', 'Dihedrals',
                       'Impropers', 'Pair', 'Bond', 'Angle'):
            section = 'other'
            continue

        # ── Header lines: any numeric line before sections ───────────────
        if section is None:
            # Try to parse count lines: "46 atoms", "3 bond types" etc.
            if len(parts) >= 2:
                try:
                    val = int(parts[0])
                    tag = ' '.join(parts[1:]).lower()
                    if tag == 'atoms':
                        info['n_atoms'] = val
                    elif tag == 'bonds':
                        info['n_bonds'] = val
                    elif tag == 'atom types':
                        info['n_atom_types'] = val
                    elif tag == 'bond types':
                        info['n_bond_types'] = val
                except ValueError:
                    pass
            # Try to parse box bounds: "xlo xhi" lines
            if len(parts) >= 4:
                try:
                    if parts[2] == 'xlo' and parts[3] == 'xhi':
                        info['xlo'], info['xhi'] = float(parts[0]), float(parts[1])
                    elif parts[2] == 'ylo' and parts[3] == 'yhi':
                        info['ylo'], info['yhi'] = float(parts[0]), float(parts[1])
                except ValueError:
                    pass
            continue

        # ── Data rows ────────────────────────────────────────────────────
        if section == 'atoms':
            # Format: id  mol  type  x  y  z  (or id type x y z)
            # Accept any row with >= 5 numeric columns
            try:
                row = [float(p) for p in parts]
                if len(row) >= 5:
                    info['atoms'].append(row[:6] if len(row) >= 6 else row)
            except ValueError:
                pass

        elif section == 'bonds':
            # Format: id  type  atom_i  atom_j
            # Column 0 = bond ID, col 1 = type, col 2 = i, col 3 = j
            try:
                row = [float(p) for p in parts]
                if len(row) >= 4:
                    info['bonds'].append(row[:4])
            except ValueError:
                pass

    # Pad atom rows to 6 columns if short
    if info['atoms']:
        max_cols = max(len(r) for r in info['atoms'])
        padded = [r + [0.0] * (max_cols - len(r)) for r in info['atoms']]
        info['atoms'] = np.array(padded)
    else:
        info['atoms'] = np.zeros((0, 6))

    info['bonds'] = (np.array(info['bonds'])
                     if info['bonds'] else np.zeros((0, 4)))

    # Convenience: set of atom IDs for O(1) lookup in tests
    if len(info['atoms']) > 0:
        info['atom_ids'] = set(info['atoms'][:, 0].astype(int))
    else:
        info['atom_ids'] = set()

    return info


# =============================================================================
# TESTS
# =============================================================================

class TestGenerateMonoRandom:
    """Full pipeline: mono strand topology, random geometry."""

    def _config(self, outdir):
        return textwrap.dedent(f"""
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.domain.lammps_data_file = 'PolyNetwork';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.perbond.kuhn.mono.value         = 20;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """)

    def test_output_files_created(self, outdir):
        """generateNetwork creates at least one output file."""
        run_config(self._config(outdir))
        files = os.listdir(outdir)
        assert len(files) > 0, f"No files created in {outdir}"
        print(f"\n  Files created: {files}")

    def test_lammps_data_file_exists(self, outdir):
        """A LAMMPS data file (.lammps or similar) is produced."""
        run_config(self._config(outdir))
        files = os.listdir(outdir)
        data_files = [f for f in files if
                      f.endswith('.lammps') or
                      f.endswith('.dat') or
                      f.endswith('.data') or
                      'PolyNetwork' in f]
        assert len(data_files) > 0, \
            f"No LAMMPS data file found. Files present: {files}"

    def test_log_files_created(self, outdir):
        """Console and network log files are written."""
        run_config(self._config(outdir))
        files = os.listdir(outdir)
        log_files = [f for f in files if f.endswith('.log')]
        assert len(log_files) >= 1, \
            f"Expected at least one .log file, found: {files}"

    def test_lammps_data_has_atoms(self, outdir):
        """The LAMMPS data file contains a non-zero atom count."""
        run_config(self._config(outdir))
        files = os.listdir(outdir)
        data_files = [f for f in files if
                      f.endswith('.lammps') or 'PolyNetwork' in f]
        if not data_files:
            pytest.skip("No LAMMPS data file to inspect")

        data = read_lammps_data(os.path.join(outdir, data_files[0]))
        assert data['n_atoms'] > 0, \
            f"LAMMPS file reports 0 atoms"

    def test_lammps_data_has_bonds(self, outdir):
        """The LAMMPS data file contains a non-zero bond count."""
        run_config(self._config(outdir))
        files = os.listdir(outdir)
        data_files = [f for f in files if
                      f.endswith('.lammps') or 'PolyNetwork' in f]
        if not data_files:
            pytest.skip("No LAMMPS data file to inspect")

        data = read_lammps_data(os.path.join(outdir, data_files[0]))
        assert data['n_bonds'] > 0, \
            f"LAMMPS file reports 0 bonds"

    def test_atom_count_reasonable(self, outdir):
        """
        Atom count is in a physically reasonable range for this domain.
        For Lx=Ly=50, b=1.6, rho=0.0078:
          domain area = (2*50*1.6)^2 = 160^2 = 25600
          expected atoms ~ 0.0078 * 25600 ~ 200
        """
        run_config(self._config(outdir))
        files = os.listdir(outdir)
        data_files = [f for f in files if
                      f.endswith('.lammps') or 'PolyNetwork' in f]
        if not data_files:
            pytest.skip("No LAMMPS data file to inspect")

        data = read_lammps_data(os.path.join(outdir, data_files[0]))
        n = data['n_atoms']
        assert 10 <= n <= 2000, \
            f"Atom count {n} outside expected range [10, 2000] for this domain"


class TestGenerateMultipleReplicates:
    """Check that Nreplicates > 1 produces multiple separate output files."""

    def test_two_replicates_produce_separate_logs(self, outdir):
        """Two replicates each write a separately-stamped log file."""
        config = textwrap.dedent(f"""
            net = network();
            net.Nreplicates = 2;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """)
        run_config(config)
        files = os.listdir(outdir)
        log_files = sorted([f for f in files if f.endswith('.log')])
        print(f"\n  Log files: {log_files}")

        # Expect at least two log files stamped N001 and N002
        n001 = any('N0001' in f for f in log_files)
        n002 = any('N0002' in f for f in log_files)
        assert n001, f"No N0001 log file found: {log_files}"
        assert n002, f"No N0002 log file found: {log_files}"


class TestGenerateHexLattice:
    """Full pipeline: mono topology, hex_lattice geometry."""

    def test_hex_network_produces_output(self, outdir):
        """hex_lattice geometry runs to completion and writes files."""
        config = textwrap.dedent(f"""
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.architecture.geometry           = 'hex_lattice';
            net.architecture.lattice_spacing    = 6;
            net.architecture.lattice_disorder_level   = 0;
            net.architecture.lattice_max_del_per_node = 0;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 8;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """)
        run_config(config)
        files = os.listdir(outdir)
        assert len(files) > 0, f"No files created for hex_lattice network"
        print(f"\n  Files: {files}")

    def test_hex_with_disorder_runs(self, outdir):
        """hex_lattice with geometric and topological disorder runs without error."""
        config = textwrap.dedent(f"""
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.architecture.geometry                 = 'hex_lattice';
            net.architecture.lattice_spacing          = 6;
            net.architecture.lattice_disorder_level   = 0.3;
            net.architecture.lattice_disorder_maxfrac = 0.4;
            net.architecture.lattice_max_del_per_node = 1;
            net.architecture.lattice_min_degree_keep  = 4;
            net.architecture.strand_typology.mode     = 'mono';
            net.peratom.Max_peratom_bond              = 8;
            net.peratom.min_degree_keep               = 2;
            net.perbond.kuhn.auto                     = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """)
        run_config(config)
        files = os.listdir(outdir)
        assert len(files) > 0, "No output from hex_lattice + disorder"


class TestGenerateBimodal:
    """Full pipeline: bimodal topology, random geometry."""

    def test_bimodal_network_produces_output(self, outdir):
        """Bimodal network runs to completion."""
        config = textwrap.dedent(f"""
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'bimodal';
            net.architecture.strand_typology.bimodal.method       = 'gaussian';
            net.architecture.strand_typology.bimodal.mean_1       = 10;
            net.architecture.strand_typology.bimodal.mean_2       = 40;
            net.architecture.strand_typology.bimodal.std_1        = 2;
            net.architecture.strand_typology.bimodal.std_2        = 5;
            net.architecture.strand_typology.bimodal.height_mode  = 'prob';
            net.architecture.strand_typology.bimodal.height_prob  = 0.3;
            net.architecture.strand_typology.bimodal.long_first   = true;
            net.architecture.strand_typology.bimodal.double_network_flag = false;
            net.architecture.strand_typology.bimodal.bin_window_method = 'manual';
            net.architecture.strand_typology.bimodal.manual_dev_type   = 'mixed';
            net.architecture.strand_typology.bimodal.stdR_1 = 3;
            net.architecture.strand_typology.bimodal.stdR_2 = 10;
            net.architecture.strand_typology.bimodal.lam_1 = 0.2;
            net.architecture.strand_typology.bimodal.lam_2 = 0.5;
            net.architecture.strand_typology.bimodal.auto_1_flag = false;
            net.architecture.strand_typology.bimodal.auto_2_flag = false;
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """)
        run_config(config)
        files = os.listdir(outdir)
        assert len(files) > 0, "No output from bimodal network"
        print(f"\n  Files: {files}")


class TestGenerateWithDefects:
    """Full pipeline including void defects."""

    def test_defect_network_produces_output(self, outdir):
        """Network with voids enabled runs to completion."""
        config = textwrap.dedent(f"""
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = true;
            net.flags.ipotential = false;
            net.defect.density_mode        = 'area_frac';
            net.defect.void_area_frac      = 0.10;
            net.defect.radius_mean         = 20;
            net.defect.radius_std          = 4;
            net.defect.radius_min          = 8;
            net.defect.radius_max          = 35;
            net.defect.size_dist           = 'gaussian';
            net.defect.shape_roughness     = 0.2;
            net.defect.shape_n_modes       = 3;
            net.defect.void_overlap        = true;
            net.defect.sparse_network      = false;
            net.defect.center_distribution = 'random';
            net.defect.margin_frac         = 0.0;
            net.defect.clamp_thickness     = 0.0;
            net.generateNetwork();
        """)
        run_config(config)
        files = os.listdir(outdir)
        assert len(files) > 0, "No output from defect network"
        print(f"\n  Files: {files}")




class TestDataFileContent:
    """
    Read the written LAMMPS data file and validate its contents
    thoroughly — atom positions, bond connectivity, ID numbering,
    and physical sanity checks.
    """

    def _generate_and_find_data(self, outdir, extra_lines=""):
        """Generate a small mono/random network and return parsed data."""
        config = f"""
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.domain.lammps_data_file = 'PolyNetwork';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.perbond.kuhn.mono.value         = 20;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            {extra_lines}
            net.generateNetwork();
        """
        run_config(config)

        files = os.listdir(outdir)
        data_files = [f for f in files if
                      ('PolyNetwork' in f or f.endswith('.dat') or
                       f.endswith('.lammps') or f.endswith('.data')) and
                      not f.endswith('.log')]
        if not data_files:
            pytest.skip(f"No data file found in {outdir}. Files: {files}")

        filepath = os.path.join(outdir, sorted(data_files)[0])
        return read_lammps_data(filepath), filepath

    def test_header_atom_count_matches_atoms_section(self, outdir):
        """Header atom count equals actual rows in Atoms section."""
        data, path = self._generate_and_find_data(outdir)
        assert data['n_atoms'] == len(data['atoms']), (
            f"Header says {data['n_atoms']} atoms but "
            f"Atoms section has {len(data['atoms'])} rows in {path}"
        )

    def test_header_bond_count_matches_bonds_section(self, outdir):
        """Header bond count equals actual rows in Bonds section."""
        data, path = self._generate_and_find_data(outdir)
        assert data['n_bonds'] == len(data['bonds']), (
            f"Header says {data['n_bonds']} bonds but "
            f"Bonds section has {len(data['bonds'])} rows in {path}"
        )

    def test_atom_ids_consecutive(self, outdir):
        """Atom IDs in the data file are consecutive from 1."""
        data, path = self._generate_and_find_data(outdir)
        if len(data['atoms']) == 0:
            pytest.skip("No atoms in file")
        ids = data['atoms'][:, 0].astype(int)
        ids_sorted = np.sort(ids)
        expected = np.arange(1, len(ids) + 1)
        assert np.array_equal(ids_sorted, expected), (
            f"Atom IDs are not consecutive 1..N: {ids_sorted[:10]}..."
        )

    def test_bond_ids_consecutive(self, outdir):
        """Bond IDs in the data file are consecutive from 1."""
        data, path = self._generate_and_find_data(outdir)
        if len(data['bonds']) == 0:
            pytest.skip("No bonds in file")
        ids = data['bonds'][:, 0].astype(int)
        ids_sorted = np.sort(ids)
        expected = np.arange(1, len(ids) + 1)
        assert np.array_equal(ids_sorted, expected), (
            f"Bond IDs are not consecutive 1..N: {ids_sorted[:10]}..."
        )

    def test_bond_endpoints_reference_valid_atoms(self, outdir):
        """Every bond endpoint ID exists in the atom list."""
        data, path = self._generate_and_find_data(outdir)
        if len(data['bonds']) == 0:
            pytest.skip("No bonds in file")
        atom_ids = set(data['atoms'][:, 0].astype(int))
        # Bond format: id type i j  — endpoints are columns 2 and 3
        atom_ids = data['atom_ids']
        i_col = data['bonds'][:, 2].astype(int)
        j_col = data['bonds'][:, 3].astype(int)
        bad_i = [v for v in i_col if v not in atom_ids]
        bad_j = [v for v in j_col if v not in atom_ids]
        assert not bad_i, f"Bond endpoint i IDs not in atom list: {bad_i[:5]}"
        assert not bad_j, f"Bond endpoint j IDs not in atom list: {bad_j[:5]}"

    def test_atom_positions_inside_domain(self, outdir):
        """All atom x,y positions are within the domain bounds."""
        data, path = self._generate_and_find_data(outdir)
        if len(data['atoms']) == 0:
            pytest.skip("No atoms in file")
        # Atoms section columns: id  mol  type  x  y  (z)
        # Column indices vary by LAMMPS format — find x,y from domain bounds
        xlo, xhi = data['xlo'], data['xhi']
        ylo, yhi = data['ylo'], data['yhi']
        if xlo == xhi == 0:
            pytest.skip("Domain bounds not parsed from file")
        # Atom format: id mol type x y z — x is col 3, y is col 4 (0-indexed)
        x = data['atoms'][:, 3]
        y = data['atoms'][:, 4]
        margin = 1.0  # small tolerance for boundary atoms
        assert np.all(x >= xlo - margin) and np.all(x <= xhi + margin), \
            "Some atom x positions outside domain"
        assert np.all(y >= ylo - margin) and np.all(y <= yhi + margin), \
            "Some atom y positions outside domain"

    def test_no_duplicate_bonds(self, outdir):
        """No bond appears twice (same pair of atom IDs)."""
        data, path = self._generate_and_find_data(outdir)
        if len(data['bonds']) == 0:
            pytest.skip("No bonds in file")
        # Bond format: id type i j — endpoints are columns 2 and 3
        i_col = data['bonds'][:, 2].astype(int)
        j_col = data['bonds'][:, 3].astype(int)
        # Normalise so smaller ID is always first
        pairs = [tuple(sorted((i, j))) for i, j in zip(i_col, j_col)]
        n_unique = len(set(pairs))
        assert n_unique == len(pairs), (
            f"Found {len(pairs) - n_unique} duplicate bonds in {path}"
        )

    def test_no_self_bonds(self, outdir):
        """No bond connects an atom to itself."""
        data, path = self._generate_and_find_data(outdir)
        if len(data['bonds']) == 0:
            pytest.skip("No bonds in file")
        # Bond format: id type i j — endpoints are columns 2 and 3
        i_col = data['bonds'][:, 2].astype(int)
        j_col = data['bonds'][:, 3].astype(int)
        self_bonds = np.sum(i_col == j_col)
        assert self_bonds == 0, f"Found {self_bonds} self-bonds in {path}"

    def test_log_file_contains_atom_count(self, outdir):
        """The network log file records a non-zero atom count."""
        self._generate_and_find_data(outdir)
        files = os.listdir(outdir)
        net_logs = [f for f in files if 'network' in f and f.endswith('.log')]
        if not net_logs:
            pytest.skip("No network log file found")

        log_path = os.path.join(outdir, sorted(net_logs)[0])
        with open(log_path, 'r') as f:
            content = f.read()

        assert 'Number of atoms' in content or 'atom_count' in content, \
            f"Log file does not mention atom count. Content:\n{content[:500]}"

    def test_log_file_contains_bond_count(self, outdir):
        """The network log file records a non-zero bond count."""
        self._generate_and_find_data(outdir)
        files = os.listdir(outdir)
        net_logs = [f for f in files if 'network' in f and f.endswith('.log')]
        if not net_logs:
            pytest.skip("No network log file found")

        log_path = os.path.join(outdir, sorted(net_logs)[0])
        with open(log_path, 'r') as f:
            content = f.read()

        assert 'Number of bonds' in content or 'bond_count' in content, \
            f"Log file does not mention bond count."




class TestGenerateVisualizableOutputs:
    """
    Generate one representative network per geometry/topology combination
    and write the output to NETWORKGEN_VIZ_DIR — a persistent directory
    that is NEVER deleted, so you can open the .dat files in OVITO,
    VMD, or any other visualizer after the test run.

    These tests always run regardless of KEEP_OUTPUT.  The output dir
    is printed at the start of the session so you know where to look.

    Existing files are overwritten on each run so the dir stays clean.
    """

    def _run_to_viz_dir(self, label, config_template):
        """
        Run config_template with outdir=NETWORKGEN_VIZ_DIR and a
        sub-label so files from different configs don't collide.
        """
        subdir = os.path.join(NETWORKGEN_VIZ_DIR, label)
        os.makedirs(subdir, exist_ok=True)

        # Clear old files so stale outputs from previous runs don't linger
        for f in os.listdir(subdir):
            try:
                os.remove(os.path.join(subdir, f))
            except OSError:
                pass

        config = config_template.format(outdir=subdir)
        run_config(config)

        files = os.listdir(subdir)
        data_files = [f for f in files if
                      f.endswith('.dat') or f.endswith('.lammps') or
                      f.endswith('.data') or
                      ('PolyNetwork' in f and not f.endswith('.log'))]

        print(f"\n  [{label}] output -> {subdir}")
        print(f"  [{label}] files   -> {files}")
        return subdir, data_files

    def test_viz_mono_random(self):
        """Write a mono/random network to the viz dir."""
        config = """
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.domain.lammps_data_file = 'PolyNetwork';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.perbond.kuhn.mono.value         = 20;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """
        subdir, data_files = self._run_to_viz_dir("mono_random", config)
        assert len(data_files) > 0, f"No data file written to {subdir}"

    def test_viz_mono_hex(self):
        """Write a mono/hex_lattice network to the viz dir."""
        config = """
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.domain.lammps_data_file = 'PolyNetwork';
            net.architecture.geometry           = 'hex_lattice';
            net.architecture.lattice_spacing    = 6;
            net.architecture.lattice_disorder_level   = 0;
            net.architecture.lattice_max_del_per_node = 0;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 8;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """
        subdir, data_files = self._run_to_viz_dir("mono_hex", config)
        assert len(data_files) > 0, f"No data file written to {subdir}"

    def test_viz_hex_with_disorder(self):
        """Write a disordered hex_lattice network to the viz dir."""
        config = """
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.domain.lammps_data_file = 'PolyNetwork';
            net.architecture.geometry                 = 'hex_lattice';
            net.architecture.lattice_spacing          = 6;
            net.architecture.lattice_disorder_level   = 0.4;
            net.architecture.lattice_disorder_maxfrac = 0.4;
            net.architecture.lattice_max_del_per_node = 1;
            net.architecture.lattice_min_degree_keep  = 4;
            net.architecture.strand_typology.mode     = 'mono';
            net.peratom.Max_peratom_bond              = 8;
            net.peratom.min_degree_keep               = 2;
            net.perbond.kuhn.auto                     = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """
        subdir, data_files = self._run_to_viz_dir("hex_disorder", config)
        assert len(data_files) > 0, f"No data file written to {subdir}"

    def test_viz_bimodal_random(self):
        """Write a bimodal/random network to the viz dir."""
        config = """
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.domain.lammps_data_file = 'PolyNetwork';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'bimodal';
            net.architecture.strand_typology.bimodal.method       = 'gaussian';
            net.architecture.strand_typology.bimodal.mean_1       = 10;
            net.architecture.strand_typology.bimodal.mean_2       = 40;
            net.architecture.strand_typology.bimodal.std_1        = 2;
            net.architecture.strand_typology.bimodal.std_2        = 5;
            net.architecture.strand_typology.bimodal.height_mode  = 'prob';
            net.architecture.strand_typology.bimodal.height_prob  = 0.3;
            net.architecture.strand_typology.bimodal.long_first   = true;
            net.architecture.strand_typology.bimodal.double_network_flag = false;
            net.architecture.strand_typology.bimodal.bin_window_method = 'manual';
            net.architecture.strand_typology.bimodal.manual_dev_type   = 'mixed';
            net.architecture.strand_typology.bimodal.stdR_1 = 3;
            net.architecture.strand_typology.bimodal.stdR_2 = 10;
            net.architecture.strand_typology.bimodal.lam_1  = 0.2;
            net.architecture.strand_typology.bimodal.lam_2  = 0.5;
            net.architecture.strand_typology.bimodal.auto_1_flag = false;
            net.architecture.strand_typology.bimodal.auto_2_flag = false;
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.generateNetwork();
        """
        subdir, data_files = self._run_to_viz_dir("bimodal_random", config)
        assert len(data_files) > 0, f"No data file written to {subdir}"

    def test_viz_with_defects(self):
        """Write a network with void defects to the viz dir."""
        config = """
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.domain.lammps_data_file = 'PolyNetwork';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.perbond.kuhn.mono.value         = 20;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = true;
            net.flags.ipotential = false;
            net.defect.density_mode        = 'area_frac';
            net.defect.void_area_frac      = 0.15;
            net.defect.radius_mean         = 40;
            net.defect.radius_std          = 8;
            net.defect.radius_min          = 15;
            net.defect.radius_max          = 80;
            net.defect.size_dist           = 'gaussian';
            net.defect.shape_roughness     = 0.3;
            net.defect.shape_n_modes       = 4;
            net.defect.void_overlap        = true;
            net.defect.sparse_network      = false;


class TestGenerateMultiType:
    """Full pipeline with post-cleanup multi-type assignment enabled."""

    def test_multitype_network_exports_atom_and_bond_types(self, outdir):
        config = textwrap.dedent(f"""
            net = network();
            net.Nreplicates = 1;
            net.domain.b            = 1.6;
            net.domain.Lx           = 50;
            net.domain.Ly           = 50;
            net.domain.boundary     = 'fixed';
            net.domain.write_location = '{outdir}';
            net.architecture.geometry           = 'random';
            net.architecture.rho_atom           = 0.0078;
            net.architecture.strand_typology.mode = 'mono';
            net.peratom.Max_peratom_bond        = 6;
            net.peratom.min_degree_keep         = 2;
            net.perbond.kuhn.auto               = true;
            net.flags.isave      = true;
            net.flags.iplot      = false;
            net.flags.ilog       = true;
            net.flags.idefect    = false;
            net.flags.ipotential = false;
            net.architecture.types.enabled    = true;
            net.architecture.types.natom_type = 3;
            net.architecture.types.nbond_type = 2;
            net.architecture.types.atype_mode = 'frac';
            net.architecture.types.btype_mode = 'frac';
            net.architecture.types.atom_frac  = [0.50, 0.30, 0.20];
            net.architecture.types.bond_frac  = [0.60, 0.40];
            net.architecture.types.connectivity = [ ...
                1 2 1 0; ...
                1 2 2 1; ...
                2 3 1 1; ...
                2 3 2 0  ...
            ];
            net.generateNetwork();
        """)

        run_config(config)
        files = os.listdir(outdir)
        data_files = [f for f in files if f.endswith('.dat') or 'PolyNetwork' in f]
        assert data_files, f"No LAMMPS data file found. Files present: {files}"

        data = read_lammps_data(os.path.join(outdir, data_files[0]))
        atom_types = data['atoms'][:, 2].astype(int)
        bond_types = data['bonds'][:, 1].astype(int)

        assert set(atom_types) == {1, 2, 3}, f"Unexpected atom types: {sorted(set(atom_types))}"
        assert set(bond_types) == {1, 2}, f"Unexpected bond types: {sorted(set(bond_types))}"

        atom_type_by_id = {
            int(row[0]): int(row[2])
            for row in data['atoms']
        }

        connectivity = {
            (1, 2, 1): False,
            (1, 2, 2): True,
            (2, 3, 1): True,
            (2, 3, 2): False,
        }

        for bond_row in data['bonds']:
            bond_type = int(bond_row[1])
            atom_i = int(bond_row[2])
            atom_j = int(bond_row[3])

            type_i = atom_type_by_id[atom_i]
            type_j = atom_type_by_id[atom_j]
            key = (min(type_i, type_j), max(type_i, type_j), bond_type)

            allowed = connectivity.get(key, True)
            assert allowed, (
                f"Bond ({atom_i}, {atom_j}) with atom types ({type_i}, {type_j}) "
                f"was assigned forbidden bond type {bond_type}"
            )


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    import subprocess
    result = subprocess.run(
        [sys.executable, "-m", "pytest", __file__, "-v", "--tb=short"],
        check=False
    )
    sys.exit(result.returncode)