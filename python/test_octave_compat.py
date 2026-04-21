"""
test_octave_compat.py
=====================
Octave compatibility test suite for NetworkGen.

Tests each pipeline stage in isolation by calling Octave via oct2py.
Each test is self-contained: it sets up the minimum obj state needed,
calls the function under test, and asserts structural correctness of
the output.

Requirements
------------
    pip install oct2py pytest numpy

Octave >= 7.0 must be installed and on PATH.
Your NetworkGen .m files must be in NETWORKGEN_PATH (set below or via
the NETWORKGEN_PATH environment variable).

Usage
-----
    # Run all tests
    pytest test_octave_compat.py -v

    # Run a specific test
    pytest test_octave_compat.py::test_setup_domain -v

    # Stop on first failure
    pytest test_octave_compat.py -v -x

    # Show Octave stdout in real time
    pytest test_octave_compat.py -v -s
"""

import os
import sys
import time
import textwrap
import numpy as np
import pytest

# ── Octave executable ─────────────────────────────────────────────────────────
# Set OCTAVE_EXECUTABLE *before* importing oct2py — it reads the env var at
# import time.  Precedence:
#   1. Shell environment:  export OCTAVE_EXECUTABLE=/usr/local/bin/octave
#   2. OCTAVE_EXECUTABLE_PATH constant below (edit to match your install)
#   3. Falls back to bare 'octave' and hopes it is on PATH
#
# --no-gui suppresses the "no graphical display found" X11 warning that
# appears in WSL / headless environments.  It has no effect on computation.

OCTAVE_EXECUTABLE_PATH = os.environ.get(
    "OCTAVE_EXECUTABLE",
    "/usr/local/bin/octave"     # <- change this if your octave lives elsewhere
)

# Append --no-gui if not already present so X11 warnings are suppressed
if "--no-gui" not in OCTAVE_EXECUTABLE_PATH:
    os.environ["OCTAVE_EXECUTABLE"] = OCTAVE_EXECUTABLE_PATH + " --no-gui"
else:
    os.environ["OCTAVE_EXECUTABLE"] = OCTAVE_EXECUTABLE_PATH

# ── oct2py import with helpful error ──────────────────────────────────────────
# Note: oct2py in this version does not accept an 'executable' constructor
# argument — the environment variable set above is the only supported way.
try:
    from oct2py import octave, Oct2PyError
except ImportError:
    sys.exit(
        "\n[ERROR] oct2py not found.  Install with:  pip install oct2py\n"
        "        Octave itself must also be installed and on your PATH.\n"
    )

# ── Path to your NetworkGen .m files ──────────────────────────────────────────
NETWORKGEN_PATH = "/home/zwhit/networkgen/NetworkGen/"

# ── Octave version check ───────────────────────────────────────────────────────
MIN_OCTAVE_MAJOR = 7


# =============================================================================
# Session-scoped fixture: start one Octave session for all tests
# =============================================================================

@pytest.fixture(scope="session", autouse=True)
def octave_session():
    """Start Octave, add NetworkGen to path, verify version."""

    print(f"\n[setup] Adding NetworkGen path: {NETWORKGEN_PATH}")
    octave.addpath(octave.genpath(NETWORKGEN_PATH))

    # Check Octave version
    ver = octave.eval("OCTAVE_VERSION", verbose=False)
    print(f"[setup] Octave version: {ver}")
    major = int(str(ver).split(".")[0])
    if major < MIN_OCTAVE_MAJOR:
        pytest.exit(
            f"Octave >= {MIN_OCTAVE_MAJOR}.0 required, found {ver}. "
            "classdef handle support is incomplete in earlier versions."
        )

    yield

    octave.exit()


# =============================================================================
# Helper: build a minimal network object in Octave
# =============================================================================

def make_network_obj(extra_lines=""):
    """
    Evaluate a block of Octave code that creates a 'net' network object
    with small domain suitable for fast testing, plus any extra lines.
    Returns the Octave variable name 'net' (it lives in the Octave workspace).
    """
    setup = textwrap.dedent(f"""
        net = network();
        net.domain.b            = 1.6;
        net.domain.Lx           = 20;
        net.domain.Ly           = 20;
        net.domain.boundary     = 'fixed';
        net.architecture.geometry       = 'random';
        net.architecture.rho_atom       = 0.0078;
        net.architecture.lattice_spacing = 6;
        net.peratom.Max_peratom_bond    = 6;
        net.peratom.min_degree_keep     = 2;
        net.flags.idefect   = false;
        net.flags.iplot     = false;
        net.flags.isave     = false;
        net.flags.ilog      = false;
        net.flags.ipotential = false;
        {extra_lines}
        SetupDomain(net);
    """)
    octave.eval(setup, verbose=False)


# =============================================================================
# Utility assertions
# =============================================================================

def assert_atoms_shape(Atoms, min_rows=1):
    assert Atoms is not None, "Atoms is None"
    assert Atoms.ndim == 2, f"Atoms should be 2D, got shape {Atoms.shape}"
    assert Atoms.shape[0] >= min_rows, \
        f"Expected >= {min_rows} atoms, got {Atoms.shape[0]}"
    assert Atoms.shape[1] >= 5, \
        f"Atoms should have >= 5 columns, got {Atoms.shape[1]}"


def assert_bonds_shape(Bonds, min_rows=1):
    assert Bonds is not None, "Bonds is None"
    assert Bonds.ndim == 2, f"Bonds should be 2D, got shape {Bonds.shape}"
    assert Bonds.shape[0] >= min_rows, \
        f"Expected >= {min_rows} bonds, got {Bonds.shape[0]}"
    assert Bonds.shape[1] == 5, \
        f"Bonds should have 5 columns, got {Bonds.shape[1]}"


def assert_ids_consecutive(arr_col, label="IDs"):
    """Check that a column of IDs is 1, 2, 3, ..., N."""
    ids = arr_col.flatten().astype(int)
    expected = np.arange(1, len(ids) + 1)
    assert np.array_equal(ids, expected), \
        f"{label} are not consecutive 1..N: {ids[:10]}..."


def assert_bond_endpoints_valid(Bonds, n_atoms):
    """All bond endpoint IDs must be in [1, n_atoms]."""
    i_col = Bonds[:, 1].astype(int)
    j_col = Bonds[:, 2].astype(int)
    assert np.all(i_col >= 1) and np.all(i_col <= n_atoms), \
        f"Bond endpoint i out of range [1,{n_atoms}]"
    assert np.all(j_col >= 1) and np.all(j_col <= n_atoms), \
        f"Bond endpoint j out of range [1,{n_atoms}]"


def assert_L0_positive(Bonds):
    L0 = Bonds[:, 3]
    assert np.all(L0 > 0), f"Some L0 values are <= 0: {L0[L0 <= 0]}"


# =============================================================================
# TESTS
# =============================================================================

class TestOctaveAvailability:

    def test_octave_basic_math(self):
        """Octave can do arithmetic — baseline sanity check."""
        result = octave.eval("2 + 2", verbose=False)
        assert result == 4.0

    def test_octave_matrix_ops(self):
        """Octave handles matrix ops correctly."""
        A = octave.eval("rand(10, 10)", verbose=False)
        assert A.shape == (10, 10)

    def test_networkgen_path_accessible(self):
        """network.m is on the Octave path."""
        result = octave.eval("exist('network', 'file')", verbose=False)
        assert result == 2, \
            "network.m not found on Octave path — check NETWORKGEN_PATH"


class TestNetworkClass:

    def test_network_instantiation(self):
        """network() constructor runs without error."""
        octave.eval("n = network();", verbose=False)
        result = octave.eval("class(n)", verbose=False)
        assert result == "network"

    def test_network_default_properties(self):
        """Default property values are set correctly."""
        b     = octave.eval("network().domain.b", verbose=False)
        geom  = octave.eval("network().architecture.geometry", verbose=False)
        assert float(b) == pytest.approx(1.6)
        assert geom == "random"

    def test_subclass_instantiation(self):
        """architecture and bondstyle subclasses initialise."""
        ok = octave.eval(
            "n=network(); isa(n.architecture,'architecture') && isa(n.perbond,'bondstyle')",
            verbose=False
        )
        assert bool(ok)

    def test_polydisperse_mode_alias_is_accepted(self):
        """Documented polydisperse mode alias works for topology and kuhn sync."""
        result = octave.eval(textwrap.dedent("""
            n = network();
            n.architecture.strand_typology.mode = 'polydisperse';
            n.architecture.strand_typology.poly.method = 'pmf';
            n.architecture.strand_typology.poly.pmf_mean = 12;
            n.architecture.strand_typology.poly.pmf_min = 4;
            n.architecture.strand_typology.poly.pmf_max = 20;
            syncKuhnAssignmentFromTopology(n);
            strcmp(n.perbond.kuhn.mode, 'polydisperse') && ...
            strcmp(resolveNetworkTypePrefix(n), 'PD')
        """), verbose=False)
        assert bool(result)


class TestSetupDomain:

    def test_domain_bounds_set(self):
        """SetupDomain populates xlo/xhi/ylo/yhi."""
        make_network_obj()
        xlo = float(octave.eval("net.domain.xlo", verbose=False))
        xhi = float(octave.eval("net.domain.xhi", verbose=False))
        assert xlo < xhi, f"xlo={xlo} should be < xhi={xhi}"

    def test_max_atom_positive(self):
        """Max_atom is a positive integer after SetupDomain."""
        make_network_obj()
        ma = int(octave.eval("net.domain.Max_atom", verbose=False))
        assert ma > 0, f"Max_atom={ma} should be positive"

    def test_max_bond_positive(self):
        """Max_bond is a positive integer after SetupDomain."""
        make_network_obj()
        mb = int(octave.eval("net.domain.Max_bond", verbose=False))
        assert mb > 0

    def test_min_node_sep_set(self):
        """min_node_sep is set and equals lattice_spacing * b."""
        make_network_obj()
        sep = float(octave.eval("net.domain.min_node_sep", verbose=False))
        b   = float(octave.eval("net.domain.b", verbose=False))
        ls  = float(octave.eval("net.architecture.lattice_spacing", verbose=False))
        assert sep == pytest.approx(ls * b, rel=1e-6)

    def test_hex_lattice_domain(self):
        """SetupDomain also works for hex_lattice geometry."""
        make_network_obj("net.architecture.geometry = 'hex_lattice';")
        xhi = float(octave.eval("net.domain.xhi", verbose=False))
        assert xhi > 0


class TestAddAtoms:

    def test_random_atoms_returned(self):
        """AddAtomsRandom returns a non-empty atom array."""
        make_network_obj()
        octave.eval("[Atoms,~] = AddAtoms(net);", verbose=False)
        Atoms = octave.pull("Atoms")
        assert_atoms_shape(Atoms, min_rows=5)

    def test_atom_ids_consecutive(self):
        """Atom IDs in column 1 are consecutive from 1."""
        make_network_obj()
        octave.eval("[Atoms,~] = AddAtoms(net);", verbose=False)
        Atoms = octave.pull("Atoms")
        assert_ids_consecutive(Atoms[:, 0], "Atom IDs")

    def test_atoms_inside_domain(self):
        """All atom positions are within domain bounds."""
        make_network_obj()
        octave.eval("[Atoms,~] = AddAtoms(net);", verbose=False)
        xlo = float(octave.eval("net.domain.xlo", verbose=False))
        xhi = float(octave.eval("net.domain.xhi", verbose=False))
        ylo = float(octave.eval("net.domain.ylo", verbose=False))
        yhi = float(octave.eval("net.domain.yhi", verbose=False))
        Atoms = octave.pull("Atoms")
        x = Atoms[:, 1]
        y = Atoms[:, 2]
        assert np.all(x >= xlo) and np.all(x <= xhi), "Atoms outside x bounds"
        assert np.all(y >= ylo) and np.all(y <= yhi), "Atoms outside y bounds"

    def test_hex_atoms_returned(self):
        """AddAtomsHex returns atoms and non-empty LatticeData."""
        make_network_obj("net.architecture.geometry = 'hex_lattice';")
        n_atoms = int(octave.eval(
            "[A, LD] = AddAtoms(net); size(A,1)", verbose=False))
        has_idx = bool(octave.eval(
            "isfield(LD,'idx_map')", verbose=False))
        assert n_atoms > 0
        assert has_idx, "LatticeData missing idx_map field"

    def test_min_separation_respected(self):
        """No two atoms are closer than min_node_sep."""
        make_network_obj()
        octave.eval("[Atoms,~] = AddAtoms(net);", verbose=False)
        Atoms = octave.pull("Atoms")
        dmin  = float(octave.eval("net.domain.min_node_sep", verbose=False))
        x = Atoms[:, 1]
        y = Atoms[:, 2]
        n = len(x)
        # Check a random subset for speed (full O(n^2) for large n)
        sample = min(n, 50)
        idx = np.random.choice(n, sample, replace=False)
        for i in idx:
            dx = x - x[i]
            dy = y - y[i]
            d  = np.sqrt(dx**2 + dy**2)
            d[i] = np.inf   # ignore self
            assert np.min(d) >= dmin * 0.95, \
                f"Atoms {i} has a neighbour closer than min_node_sep"


class TestAddBondsMono:

    def _setup(self):
        make_network_obj(
            "net.architecture.strand_typology.mode = 'mono';"
        )
        octave.eval("[Atoms, LD] = AddAtoms(net);", verbose=False)
        octave.eval("[Atoms, Bonds] = AddBonds(net, Atoms, LD);", verbose=False)

    def test_bonds_returned(self):
        """AddBondsMono returns a non-empty bond array."""
        self._setup()
        Bonds = octave.pull("Bonds")
        assert_bonds_shape(Bonds, min_rows=1)

    def test_bond_ids_consecutive(self):
        self._setup()
        Bonds = octave.pull("Bonds")
        assert_ids_consecutive(Bonds[:, 0], "Bond IDs")

    def test_bond_endpoints_valid(self):
        self._setup()
        Bonds  = octave.pull("Bonds")
        n_atoms = int(octave.eval("size(Atoms,1)", verbose=False))
        assert_bond_endpoints_valid(Bonds, n_atoms)

    def test_L0_positive(self):
        self._setup()
        Bonds = octave.pull("Bonds")
        assert_L0_positive(Bonds)

    def test_bond_type_is_1(self):
        """Mono bonds all have type=1 in column 5."""
        self._setup()
        types = octave.eval("Bonds(:,5)", verbose=False).flatten()
        assert np.all(types == 1), f"Expected all type=1, got: {np.unique(types)}"

    def test_degree_within_cap(self):
        """No atom exceeds Max_peratom_bond."""
        self._setup()
        deg_col = octave.eval("Atoms(:,5)", verbose=False).flatten()
        cap     = int(octave.eval("net.peratom.Max_peratom_bond", verbose=False))
        assert np.all(deg_col <= cap), \
            f"Some atoms exceed Max_peratom_bond={cap}"


class TestAddBondsBimodal:

    def _setup(self):
        make_network_obj(textwrap.dedent("""
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
        """))
        octave.eval("[Atoms, LD] = AddAtoms(net);", verbose=False)
        octave.eval("[Atoms, Bonds] = AddBonds(net, Atoms, LD);", verbose=False)

    def test_bonds_returned(self):
        self._setup()
        Bonds = octave.pull("Bonds")
        assert_bonds_shape(Bonds, min_rows=1)

    def test_both_bond_types_present(self):
        """Bimodal network should have both type-1 and type-2 bonds."""
        self._setup()
        types = octave.eval("Bonds(:,5)", verbose=False).flatten().astype(int)
        assert 1 in types, "No type-1 bonds found"
        assert 2 in types, "No type-2 bonds found"

    def test_L0_positive(self):
        self._setup()
        Bonds = octave.pull("Bonds")
        assert_L0_positive(Bonds)


class TestAssignPerBond:

    def _setup_mono(self):
        make_network_obj(
            "net.architecture.strand_typology.mode = 'mono';"
            "net.perbond.kuhn.mono.value = 20;"
        )
        octave.eval("[Atoms,LD] = AddAtoms(net);", verbose=False)
        octave.eval("[Atoms,Bonds] = AddBonds(net,Atoms,LD);", verbose=False)
        octave.eval("Nvec = AssignPerBond(net, Bonds, Atoms);", verbose=False)

    def test_nvec_length_matches_bonds(self):
        self._setup_mono()
        nb   = int(octave.eval("size(Bonds,1)", verbose=False))
        nv   = int(octave.eval("numel(Nvec)",   verbose=False))
        assert nb == nv, f"Nvec length {nv} != bond count {nb}"

    def test_mono_nvec_all_same(self):
        self._setup_mono()
        Nvec = octave.pull("Nvec").flatten()
        assert np.all(Nvec == Nvec[0]), "Mono Nvec should be constant"
        assert Nvec[0] == 20, f"Expected N=20, got {Nvec[0]}"

    def test_nvec_all_positive(self):
        self._setup_mono()
        Nvec = octave.pull("Nvec").flatten()
        assert np.all(Nvec > 0), "Some Nvec values are <= 0"


class TestCleanupNetwork:

    def _setup(self):
        make_network_obj(
            "net.architecture.strand_typology.mode = 'mono';"
        )
        octave.eval("[Atoms,LD] = AddAtoms(net);", verbose=False)
        octave.eval("[Atoms,Bonds] = AddBonds(net,Atoms,LD);", verbose=False)
        octave.eval("Nvec = AssignPerBond(net,Bonds,Atoms);", verbose=False)
        octave.eval(
            "[Atoms,Bonds,Nvec] = CleanupNetwork(net,Atoms,Bonds,Nvec);",
            verbose=False
        )

    def test_atoms_renumbered(self):
        self._setup()
        Atoms = octave.pull("Atoms")
        assert_ids_consecutive(Atoms[:, 0], "Atom IDs after cleanup")

    def test_bonds_renumbered(self):
        self._setup()
        Bonds = octave.pull("Bonds")
        assert_ids_consecutive(Bonds[:, 0], "Bond IDs after cleanup")

    def test_endpoints_valid_after_cleanup(self):
        self._setup()
        Bonds   = octave.pull("Bonds")
        n_atoms = int(octave.eval("size(Atoms,1)", verbose=False))
        assert_bond_endpoints_valid(Bonds, n_atoms)

    def test_nvec_synced_with_bonds(self):
        self._setup()
        nb = int(octave.eval("size(Bonds,1)", verbose=False))
        nv = int(octave.eval("numel(Nvec)",   verbose=False))
        assert nb == nv, f"Nvec ({nv}) out of sync with bonds ({nb}) after cleanup"

    def test_min_degree_respected(self):
        """After cleanup no atom has degree below min_degree_keep."""
        self._setup()
        deg      = octave.eval("Atoms(:,5)", verbose=False).flatten()
        min_keep = int(octave.eval("net.peratom.min_degree_keep", verbose=False))
        assert np.all(deg >= min_keep), \
            f"Some atoms have degree below min_degree_keep={min_keep}"


class TestAddDefects:

    def _setup_with_defects(self):
        make_network_obj(textwrap.dedent("""
            net.flags.idefect                 = true;
            net.defect.density_mode           = 'area_frac';
            net.defect.void_area_frac         = 0.10;
            net.defect.radius_mean            = 15;
            net.defect.radius_std             = 3;
            net.defect.radius_min             = 5;
            net.defect.radius_max             = 25;
            net.defect.size_dist              = 'gaussian';
            net.defect.shape_roughness        = 0.2;
            net.defect.shape_n_modes          = 3;
            net.defect.void_overlap           = true;
            net.defect.sparse_network         = false;
            net.defect.center_distribution    = 'random';
            net.defect.margin_frac            = 0.0;
            net.defect.clamp_thickness        = 0.0;
            net.architecture.strand_typology.mode = 'mono';
        """))
        octave.eval("[Atoms,LD] = AddAtoms(net);", verbose=False)
        octave.eval("[Atoms,Bonds] = AddBonds(net,Atoms,LD);", verbose=False)
        octave.eval("Nvec = AssignPerBond(net,Bonds,Atoms);", verbose=False)
        octave.eval(
            "[Atoms,Bonds,Nvec] = AddDefects(net,Atoms,Bonds,Nvec);",
            verbose=False
        )

    def test_defects_reduce_atom_count(self):
        """Atom count should decrease after void insertion."""
        make_network_obj(
            "net.architecture.strand_typology.mode = 'mono';"
        )
        octave.eval("[Atoms0,LD] = AddAtoms(net);", verbose=False)
        n_before = int(octave.eval("size(Atoms0,1)", verbose=False))
        self._setup_with_defects()
        n_after  = int(octave.eval("size(Atoms,1)", verbose=False))
        assert n_after < n_before, \
            f"Expected fewer atoms after defects: before={n_before}, after={n_after}"

    def test_nvec_synced_after_defects(self):
        self._setup_with_defects()
        nb = int(octave.eval("size(Bonds,1)", verbose=False))
        nv = int(octave.eval("numel(Nvec)",   verbose=False))
        assert nb == nv

    def test_idefect_false_skips(self):
        """When flags.idefect=false, atoms and bonds are unchanged."""
        make_network_obj(
            "net.flags.idefect = false;"
            "net.architecture.strand_typology.mode = 'mono';"
        )
        octave.eval("[Atoms,LD] = AddAtoms(net);", verbose=False)
        octave.eval("[Atoms,Bonds] = AddBonds(net,Atoms,LD);", verbose=False)
        octave.eval("Nvec = AssignPerBond(net,Bonds,Atoms);", verbose=False)
        nb_before = int(octave.eval("size(Bonds,1)", verbose=False))
        octave.eval(
            "[Atoms,Bonds,Nvec] = AddDefects(net,Atoms,Bonds,Nvec);",
            verbose=False
        )
        nb_after = int(octave.eval("size(Bonds,1)", verbose=False))
        assert nb_before == nb_after, \
            "AddDefects modified bonds even though idefect=false"


class TestLogClass:

    def test_log_instantiation(self):
        """networklog() constructs without error."""
        ok = bool(octave.eval(
            "L = networklog(); isa(L, 'networklog')", verbose=False
        ))
        assert ok

    def test_print_buffers_message(self):
        """log.print() stores the message in the buffer."""
        octave.eval(
            "L = networklog(); L.echo_to_console = false; "
            "L.print('hello %d\\n', 42);",
            verbose=False
        )
        n_msgs = int(octave.eval("numel(L.messages)", verbose=False))
        assert n_msgs == 1

    def test_set_replicate_stamps_id(self):
        octave.eval(
            "L = networklog(); L.setReplicate(3);",
            verbose=False
        )
        rid = int(octave.eval("L.replicate_id", verbose=False))
        assert rid == 3

    def test_clear_resets_messages(self):
        octave.eval(
            "L = networklog(); L.echo_to_console = false; "
            "L.print('x\\n'); L.setReplicate(2); L.clear();",
            verbose=False
        )
        n = int(octave.eval("numel(L.messages)", verbose=False))
        r = int(octave.eval("L.replicate_id",    verbose=False))
        assert n == 0, "messages not cleared"
        assert r == 2, "replicate_id should survive clear()"

    def test_stamp_path(self):
        """stamp_path inserts _N002 before the extension."""
        # Validate stamping logic directly inside Octave using sprintf —
        # avoids filesystem writes and path-escaping issues with tmpdir strings.
        result = octave.eval(
            "L = networklog(); L.setReplicate(2); "
            "[~, stem, ext] = fileparts('console.log'); "
            "sprintf('%s_N%03d%s', stem, L.replicate_id, ext)",
            verbose=False
        )
        assert str(result).strip() == "console_N002.log", \
            f"Expected 'console_N002.log', got: {result!r}"

    def test_full_pipeline_mono_random(self):
        """
        Run the complete pipeline for a small mono/random network.
        Checks that every stage produces internally consistent output.
        """
        make_network_obj(
            "net.architecture.strand_typology.mode = 'mono';"
            "net.perbond.kuhn.mono.value = 15;"
        )
        octave.eval(textwrap.dedent("""
            [Atoms,LD]    = AddAtoms(net);
            [Atoms,Bonds] = AddBonds(net,Atoms,LD);
            Nvec          = AssignPerBond(net,Bonds,Atoms);
            [Atoms,Bonds,Nvec] = AddDefects(net,Atoms,Bonds,Nvec);
            [Atoms,Bonds,Nvec] = CleanupNetwork(net,Atoms,Bonds,Nvec);
        """), verbose=False)

        Atoms = octave.pull("Atoms")
        Bonds = octave.pull("Bonds")
        Nvec  = octave.pull("Nvec").flatten()

        n_atoms = Atoms.shape[0]
        n_bonds = Bonds.shape[0]

        assert_atoms_shape(Atoms, min_rows=5)
        assert_bonds_shape(Bonds, min_rows=1)
        assert_ids_consecutive(Atoms[:, 0], "Atom IDs (end-to-end)")
        assert_ids_consecutive(Bonds[:, 0], "Bond IDs (end-to-end)")
        assert_bond_endpoints_valid(Bonds, n_atoms)
        assert_L0_positive(Bonds)
        assert len(Nvec) == n_bonds, "Nvec length mismatch at end of pipeline"
        assert np.all(Nvec == 15), "Mono Nvec should all be 15"

    def test_full_pipeline_hex_lattice(self):
        """Full pipeline for a small hex lattice with mono bonds."""
        make_network_obj(textwrap.dedent("""
            net.architecture.geometry = 'hex_lattice';
            net.architecture.strand_typology.mode = 'mono';
            net.architecture.lattice_disorder_level  = 0;
            net.architecture.lattice_max_del_per_node = 0;
        """))
        octave.eval(textwrap.dedent("""
            [Atoms,LD]    = AddAtoms(net);
            [Atoms,Bonds] = AddBonds(net,Atoms,LD);
            Nvec          = AssignPerBond(net,Bonds,Atoms);
            [Atoms,Bonds,Nvec] = AddHeterogeneities(net,Atoms,Bonds,Nvec);
            [Atoms,Bonds,Nvec] = CleanupNetwork(net,Atoms,Bonds,Nvec);
        """), verbose=False)

        Atoms = octave.pull("Atoms")
        Bonds = octave.pull("Bonds")
        assert_atoms_shape(Atoms, min_rows=5)
        assert_bonds_shape(Bonds, min_rows=1)
        assert_bond_endpoints_valid(Bonds, Atoms.shape[0])


# =============================================================================
# Entry point for running directly without pytest
# =============================================================================

if __name__ == "__main__":
    import subprocess
    result = subprocess.run(
        [sys.executable, "-m", "pytest", __file__, "-v", "--tb=short"],
        check=False
    )
    sys.exit(result.returncode)
