---
custom_edit_url: null
sidebar_position: 98
---

# Full Settings Reference

A complete flat reference of all NetworkGen settings. For detailed descriptions and examples see the individual setting pages linked in each section header.

## [Domain & Boundary](./network-setup/domain-boundary)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.domain.b` | `double` | (0, ∞) | `1.6` |
| `net.domain.Lx` | `double` | (0, ∞) | `150` |
| `net.domain.Ly` | `double` | (0, ∞) | `150` |
| `net.domain.scale` | `double` | (0, ∞) | `1` |
| `net.domain.boundary` | `string` | `'fixed'` \| `'periodic'` | `'fixed'` |
| `net.domain.seed` | `int` | [1, ∞) | `12345` |
| `net.domain.write_location` | `string` | any valid path | `'./networks'` |
| `net.domain.lammps_data_file` | `string` | any string | `'PolyNetwork'` |
| `net.domain.lammps_viz_file` | `string` | any string | `'PolyVisual'` |
| `net.domain.bond_table_file` | `string` | any string | `'bond'` |
| `net.domain.smp_number` | `int` | [1, ∞) | `1` |

## [Flags & Output](./network-setup/flags-output)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.flags.isave` | `boolean` | `true` \| `false` | `true` |
| `net.flags.iplot` | `boolean` | `true` \| `false` | `true` |
| `net.flags.ilog` | `boolean` | `true` \| `false` | `true` |
| `net.flags.savemode` | `boolean` | `true` \| `false` | `true` |
| `net.flags.imanualseed` | `boolean` | `true` \| `false` | `false` |
| `net.flags.idefect` | `boolean` | `true` \| `false` | `false` |
| `net.flags.ipotential` | `boolean` | `true` \| `false` | `true` |
| `net.flags.idumpsettings` | `boolean` | `true` \| `false` | `false` |
| `net.flags.iverbose_settings` | `boolean` | `true` \| `false` | `false` |

## [Geometry](./architecture/geometry)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.architecture.geometry` | `string` | `'random'` \| `'hex_lattice'` | `'random'` |
| `net.architecture.rho_atom` | `double` | (0, ∞) | `0.0078` |
| `net.architecture.lattice_spacing` | `double` | (0, ∞) | `6` |
| `net.architecture.spacing_multiplier_mode` | `string` | `'auto'` \| `'manual'` | `'auto'` |
| `net.architecture.spacing_multiplier` | `double` | [0, ∞) | `1.2` |
| `net.architecture.lattice_disorder_level` | `double` | [0, 1] | `1` |
| `net.architecture.lattice_disorder_maxfrac` | `double` | [0, 1] | `0.4` |
| `net.architecture.lattice_max_del_per_node` | `int` | [0, ∞) | `1` |
| `net.architecture.lattice_min_degree_keep` | `int` | [1, ∞) | `5` |

## [Peratom](./architecture/peratom)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.peratom.Max_peratom_bond` | `int` | [3, ∞) | — |
| `net.peratom.min_degree_keep` | `int` | [1, ∞) | `2` |

## [Strand Typology](./architecture/strand-typology)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.architecture.strand_typology.mode` | `string` | `'mono'` \| `'poly'` \| `'polydisperse'` \| `'bimodal'` | `'mono'` |
| `net.architecture.strand_typology.auto` | `boolean` | `true` \| `false` | `true` |

## [Assignment Modes](./assignment-modes/overview)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.architecture.strand_typology.mono.value` | `double` | (0, ∞) | `20` |
| `net.architecture.strand_typology.uniform.min_value` | `double` | (0, ∞) | `5` |
| `net.architecture.strand_typology.uniform.max_value` | `double` | (0, ∞) | `40` |
| `net.architecture.strand_typology.poly.method` | `string` | `'geom'` \| `'range'` \| `'pmf'` | `'pmf'` |
| `net.architecture.strand_typology.poly.min_value` | `double` | (0, ∞) | `1` |
| `net.architecture.strand_typology.poly.align_to_length` | `string` | `'ascend'` \| `'none'` | `'ascend'` |
| `net.architecture.strand_typology.poly.rounding` | `string` | `'round'` \| `'ceil'` \| `'floor'` | `'round'` |
| `net.architecture.strand_typology.poly.pmf_mean` | `double` | (0, ∞) | `40` |
| `net.architecture.strand_typology.poly.pmf_min` | `double` | (0, ∞) | `20` |
| `net.architecture.strand_typology.poly.pmf_max` | `double` | (0, ∞) | `120` |
| `net.architecture.strand_typology.poly.target_min` | `double` | (0, ∞) | `5` |
| `net.architecture.strand_typology.poly.target_max` | `double` | (0, ∞) | `120` |
| `net.architecture.strand_typology.bimodal.method` | `string` | `'single'` \| `'geom'` \| `'gaussian'` | `'gaussian'` |
| `net.architecture.strand_typology.bimodal.mean_1` | `double` | (0, `mean_2`) | `35` |
| `net.architecture.strand_typology.bimodal.mean_2` | `double` | (`mean_1`, ∞) | `60` |
| `net.architecture.strand_typology.bimodal.std_1` | `double` | (0, ∞) | `10` |
| `net.architecture.strand_typology.bimodal.std_2` | `double` | (0, ∞) | `5` |
| `net.architecture.strand_typology.bimodal.height_mode` | `string` | `'prob'` \| `'count'` | `'prob'` |
| `net.architecture.strand_typology.bimodal.height_prob` | `double` | [0, 1] | `1.0` |
| `net.architecture.strand_typology.bimodal.height_count` | `int` | [1, ∞) | `2` |
| `net.architecture.strand_typology.bimodal.long_first` | `boolean` | `true` \| `false` | `true` |
| `net.architecture.strand_typology.bimodal.double_network_flag` | `boolean` | `true` \| `false` | `true` |
| `net.architecture.strand_typology.bimodal.alpha` | `double` | (0, ∞) | `3.0` |
| `net.architecture.strand_typology.bimodal.auto_1_flag` | `boolean` | `true` \| `false` | `false` |
| `net.architecture.strand_typology.bimodal.auto_2_flag` | `boolean` | `true` \| `false` | `true` |
| `net.architecture.strand_typology.bimodal.bin_window_method` | `string` | `'manual'` \| `'adaptive'` | `'manual'` |
| `net.architecture.strand_typology.bimodal.manual_dev_type` | `string` | `'mixed'` \| `'kuhn'` \| `'both'` | `'mixed'` |
| `net.architecture.strand_typology.bimodal.stdR_1` | `double` | (0, ∞) | `3` |
| `net.architecture.strand_typology.bimodal.stdR_2` | `double` | (0, ∞) | `10` |
| `net.architecture.strand_typology.bimodal.lam_1` | `double` | (0, 1) | `0.2` |
| `net.architecture.strand_typology.bimodal.lam_2` | `double` | (0, 1) | `0.5` |
| `net.architecture.strand_typology.bimodal.min_value` | `int` | [1, ∞) | `1` |

## [Perbond & Potentials](./advanced-features/potentials-perbond)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.perbond.kuhn.auto` | `boolean` | `true` \| `false` | `true` |
| `net.perbond.kuhn.mode` | `string` | `'mono'` \| `'poly'` \| `'polydisperse'` \| `'bimodal'` | `'mono'` |
| `net.pot.k_LD` | `double` | (0, ∞) | `0.414` |
| `net.pot.N_rho` | `int` | [1, ∞) | `100000` |
| `net.pot.rho_min` | `double` | (0, `rho_max`) | `0.0` |
| `net.pot.rho_max` | `double` | (`rho_min`, ∞) | `500` |

## [Defects](./advanced-features/defects)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.defect.density_mode` | `string` | `'count'` \| `'area_frac'` | `'area_frac'` |
| `net.defect.n_voids` | `int` | [0, ∞) | `0` |
| `net.defect.void_area_frac` | `double` | [0, 1] | `0.75` |
| `net.defect.size_dist` | `string` | `'fixed'` \| `'gaussian'` \| `'exponential'` | `'gaussian'` |
| `net.defect.radius_mean` | `double` | (0, ∞) | `12` |
| `net.defect.radius_std` | `double` | (0, ∞) | `4` |
| `net.defect.radius_min` | `double` | (0, `radius_max`) | `2` |
| `net.defect.radius_max` | `double` | (`radius_min`, ∞) | `30` |
| `net.defect.shape_roughness` | `double` | [0, 1] | `0.3` |
| `net.defect.shape_n_modes` | `int` | [1, ∞) | `2` |
| `net.defect.void_overlap` | `boolean` | `true` \| `false` | `false` |
| `net.defect.center_distribution` | `string` | `'random'` \| `'uniform'` \| `'clustered'` \| `'manual'` | `'clustered'` |
| `net.defect.manual_centers` | `double array` | [N × 2] | `[]` |
| `net.defect.n_cluster_parents` | `int` | [1, ∞) | `2` |
| `net.defect.cluster_spread` | `double` | (0, ∞) | `10` |
| `net.defect.margin_frac` | `double` | [0, 0.5) | `0.15` |
| `net.defect.prune_isolated` | `boolean` | `true` \| `false` | `true` |
| `net.defect.sparse_network` | `boolean` | `true` \| `false` | `true` |
| `net.defect.bridge_width` | `double` | (0, ∞) | `1` |
| `net.defect.wall_thickness` | `double` | (0, ∞) | `18` |
| `net.defect.clamp_thickness` | `double` | (0, ∞) | `0.12` |

## [Multi-type Networks](./advanced-features/multi-type)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `net.architecture.types.enabled` | `boolean` | `true` \| `false` | `false` |
| `net.architecture.types.natom_type` | `int` | [1, ∞) | `1` |
| `net.architecture.types.nbond_type` | `int` | [1, ∞) | `1` |
| `net.architecture.types.atype_mode` | `string` | `'fixed'` \| `'frac'` | `'frac'` |
| `net.architecture.types.btype_mode` | `string` | `'fixed'` \| `'frac'` | `'frac'` |
| `net.architecture.types.atom_count` | `int array` | size [1 x `natom_type`] | `0` |
| `net.architecture.types.bond_count` | `int array` | size [1 x `nbond_type`] | `0` |
| `net.architecture.types.atom_frac` | `double array` | size [1 x `natom_type`], sum=1 | `1` |
| `net.architecture.types.bond_frac` | `double array` | size [1 x `nbond_type`], sum=1 | `1` |
| `net.architecture.types.connectivity` | `int array` | size [N x 4] | `[]` |
| `net.architecture.types.atype_sel_method` | `string` | `'random'` | `'random'` |
| `net.architecture.types.btype_sel_method` | `string` | `'random'` | `'random'` |
