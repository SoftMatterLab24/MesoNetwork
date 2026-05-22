%% Test script for manual defect coordinate feature
% Tests the new manual center_distribution mode for defects

clear all; close all;

%% Test 1: Manual defect placement with valid coordinates
fprintf('\n=== Test 1: Manual defect placement ===\n');

net = network();
net.flags.idefect = true;
net.flags.iplot = false;
net.flags.ilog = true;

net.domain.Lx = 100;
net.domain.Ly = 100;
net.domain.b = 1.0;

net.architecture.geometry = 'random';
net.architecture.rho_atom = 0.01;
net.architecture.strand_typology.mode = 'mono';

net.peratom.Max_peratom_bond = 6;
net.peratom.min_degree_keep = 2;

net.perbond.kuhn.auto = true;
net.perbond.kuhn.mono.value = 20;

% Manual defect settings
net.defect.density_mode = 'count';
net.defect.n_voids = 3;
net.defect.size_dist = 'gaussian';
net.defect.radius_mean = 5.0;
net.defect.radius_std = 1.0;
net.defect.radius_min = 2.0;
net.defect.radius_max = 8.0;
net.defect.center_distribution = 'manual';
net.defect.manual_centers = [25.0, 30.0; 50.0, 50.0; 75.0, 35.0];
net.defect.void_overlap = false;

% Generate atoms and bonds
[Atoms, LD] = AddAtoms(net);
fprintf('Created %d atoms\n', size(Atoms, 1));

[Atoms, Bonds] = AddBonds(net, Atoms, LD);
fprintf('Created %d bonds\n', size(Bonds, 1));

Nvec = AssignPerBond(net, Bonds, Atoms);

% Add defects
n_atoms_before = size(Atoms, 1);
n_bonds_before = size(Bonds, 1);

[Atoms, Bonds, Nvec] = AddDefects(net, Atoms, Bonds, Nvec);

n_atoms_after = size(Atoms, 1);
n_bonds_after = size(Bonds, 1);

fprintf('Atoms: %d -> %d (removed %d)\n', n_atoms_before, n_atoms_after, n_atoms_before - n_atoms_after);
fprintf('Bonds: %d -> %d (removed %d)\n', n_bonds_before, n_bonds_after, n_bonds_before - n_bonds_after);

assert(n_atoms_after < n_atoms_before, 'ERROR: Atoms should be removed by defects');
assert(n_bonds_after < n_bonds_before, 'ERROR: Bonds should be removed by defects');
fprintf('✓ Test 1 PASSED\n');

%% Test 2: Error handling - manual mode without density_mode='count'
fprintf('\n=== Test 2: Error handling - manual mode with area_frac ===\n');

net2 = network();
net2.flags.idefect = true;
net2.flags.ilog = true;

net2.domain.Lx = 100;
net2.domain.Ly = 100;
net2.domain.b = 1.0;

net2.architecture.geometry = 'random';
net2.architecture.rho_atom = 0.01;
net2.architecture.strand_typology.mode = 'mono';

net2.peratom.Max_peratom_bond = 6;
net2.peratom.min_degree_keep = 2;

net2.perbond.kuhn.auto = true;
net2.perbond.kuhn.mono.value = 20;

% Invalid: area_frac with manual distribution
net2.defect.density_mode = 'area_frac';
net2.defect.void_area_frac = 0.1;
net2.defect.center_distribution = 'manual';
net2.defect.manual_centers = [25.0, 30.0; 50.0, 50.0];

[Atoms, LD] = AddAtoms(net2);
[Atoms, Bonds] = AddBonds(net2, Atoms, LD);
Nvec = AssignPerBond(net2, Bonds, Atoms);

try
    [Atoms, Bonds, Nvec] = AddDefects(net2, Atoms, Bonds, Nvec);
    fprintf('✗ Test 2 FAILED: Should have thrown an error\n');
catch ME
    fprintf('✓ Test 2 PASSED: Error correctly thrown\n');
    fprintf('  Error message: %s\n', ME.message);
end

%% Test 3: Error handling - manual mode without manual_centers
fprintf('\n=== Test 3: Error handling - manual mode without manual_centers ===\n');

net3 = network();
net3.flags.idefect = true;
net3.flags.ilog = true;

net3.domain.Lx = 100;
net3.domain.Ly = 100;
net3.domain.b = 1.0;

net3.architecture.geometry = 'random';
net3.architecture.rho_atom = 0.01;
net3.architecture.strand_typology.mode = 'mono';

net3.peratom.Max_peratom_bond = 6;
net3.peratom.min_degree_keep = 2;

net3.perbond.kuhn.auto = true;
net3.perbond.kuhn.mono.value = 20;

% Invalid: manual mode but no coordinates
net3.defect.density_mode = 'count';
net3.defect.n_voids = 2;
net3.defect.center_distribution = 'manual';
net3.defect.manual_centers = [];

[Atoms, LD] = AddAtoms(net3);
[Atoms, Bonds] = AddBonds(net3, Atoms, LD);
Nvec = AssignPerBond(net3, Bonds, Atoms);

try
    [Atoms, Bonds, Nvec] = AddDefects(net3, Atoms, Bonds, Nvec);
    fprintf('✗ Test 3 FAILED: Should have thrown an error\n');
catch ME
    fprintf('✓ Test 3 PASSED: Error correctly thrown\n');
    fprintf('  Error message: %s\n', ME.message);
end

%% Test 4: Error handling - wrong dimensions for manual_centers
fprintf('\n=== Test 4: Error handling - wrong dimensions for manual_centers ===\n');

net4 = network();
net4.flags.idefect = true;
net4.flags.ilog = true;

net4.domain.Lx = 100;
net4.domain.Ly = 100;
net4.domain.b = 1.0;

net4.architecture.geometry = 'random';
net4.architecture.rho_atom = 0.01;
net4.architecture.strand_typology.mode = 'mono';

net4.peratom.Max_peratom_bond = 6;
net4.peratom.min_degree_keep = 2;

net4.perbond.kuhn.auto = true;
net4.perbond.kuhn.mono.value = 20;

% Invalid: coordinates count doesn't match n_voids
net4.defect.density_mode = 'count';
net4.defect.n_voids = 3;
net4.defect.center_distribution = 'manual';
net4.defect.manual_centers = [25.0, 30.0; 50.0, 50.0];  % Only 2 points, need 3

[Atoms, LD] = AddAtoms(net4);
[Atoms, Bonds] = AddBonds(net4, Atoms, LD);
Nvec = AssignPerBond(net4, Bonds, Atoms);

try
    [Atoms, Bonds, Nvec] = AddDefects(net4, Atoms, Bonds, Nvec);
    fprintf('✗ Test 4 FAILED: Should have thrown an error\n');
catch ME
    fprintf('✓ Test 4 PASSED: Error correctly thrown\n');
    fprintf('  Error message: %s\n', ME.message);
end

fprintf('\n=== All tests completed ===\n');
