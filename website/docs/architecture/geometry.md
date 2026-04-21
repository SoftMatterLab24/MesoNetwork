---
custom_edit_url: null
sidebar_position: 1
---

# Geometry

These settings control the spatial arrangement of nodes in the network.

---

### `geometry`

| Type | Args | Default |
|------|------|---------|
| `string` | `'random'` \| `'hex'` | `'random'` |

Defines the base pattern used to place nodes in the domain.

- **random** — nodes are placed stochastically, producing an amorphous network topology typical of real polymer networks
- **hex** — nodes are placed on a hexagonal lattice, producing a regular network. Disorder can be introduced via `lattice_disorder_level`

```matlab
net.architecture.geometry = 'hex_lattice';
```

---

### `lattice_spacing`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `6` |

The nominal spacing between nodes on the lattice, in units of `b`. Previously referred to as `min_node_sep`. Only relevant when `geometry = 'hex_lattice'`.

```matlab
net.architecture.lattice_spacing = 1.5;
```

---

### `spacing_multiplier_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'auto'` \| `'manual'` | `'auto'` |

Controls whether the spacing multiplier is computed automatically or set manually. Only relevant when `geometry = 'hex_lattice'`.

```matlab
net.architecture.spacing_multiplier_mode = 'manual';
```

---

### `spacing_multiplier`

| Type | Args | Default |
|------|------|---------|
| `double` | [0, ∞) | `1.2` |

A manual scaling factor applied to `lattice_spacing`. Only used when `spacing_multiplier_mode = 'manual'`.

```matlab
net.architecture.spacing_multiplier = 1.2;
```

---

### `lattice_disorder_level`

| Type | Args | Default |
|------|------|---------|
| `double` | [0, 1] | `1` |

Controls the magnitude of random positional perturbations applied to lattice nodes. `0` = perfect lattice, `1` = maximum disorder. Only relevant when `geometry = 'hex_lattice'`.

```matlab
net.architecture.lattice_disorder_level = 0.3;
```

---

### `lattice_disorder_maxfrac`

| Type | Args | Default |
|------|------|---------|
| `double` | [0, 1] | `0.4` |

Sets the maximum fractional displacement of a node from its ideal lattice position. Works in conjunction with `lattice_disorder_level`.

```matlab
net.architecture.lattice_disorder_maxfrac = 0.5;
```


---

### `lattice_max_del_per_node`

| Type | Args | Default |
|------|------|---------|
| `int` | [0, ∞) | `1` |

Maximum number of bonds that can be deleted per node when topological disorder is active. Bounded above to ensure each node retains at least 2 connections.

```matlab
net.architecture.lattice_max_del_per_node = 2;
```

---

### `lattice_min_degree_keep`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `5` |

Minimum node degree (number of bonds) to maintain during topological disorder. Nodes will not have bonds removed below this threshold.

```matlab
net.architecture.lattice_min_degree_keep = 3;
```

---

### `rho_atom`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `0.0078` |

Target atom number density in the 2D domain (atoms per unit area). Used to determine the total number of nodes to place in the domain.

```matlab
net.architecture.rho_atom = 0.5;
```
