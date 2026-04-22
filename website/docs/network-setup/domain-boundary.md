---
custom_edit_url: null
sidebar_position: 1
---

# Domain & Boundary

These settings define the physical size and shape of the 2D simulation domain, as well as the random seed and output file naming.

---

### `b`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `1.6` |

The fundamental lengthscale of the network. The in-plane domain lengths (`Lx`, `Ly`) are specified in units of `b`. Changing `b` effectively rescales the entire network.

```matlab
net.domain.b = 1.0;
```

---

### `Lx`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | — |

The x-dimension of the simulation domain in units of `b`.

```matlab
net.domain.Lx = 10;
```

---

### `Ly`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | — |

The y-dimension of the simulation domain in units of `b`.

```matlab
net.domain.Ly = 10;
```

---

### `scale`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `1.6` |

A global rescaling factor applied to the domain after generation. Useful for quickly resizing a network without reconfiguring all domain settings.

```matlab
net.domain.scale = 2.0; % doubles the domain size
```

---

### `boundary`

| Type | Args | Default |
|------|------|---------|
| `string` | `'fixed'` \| `'periodic'` | `'fixed'` |

Defines the boundary conditions of the simulation domain.

- **fixed** — nodes near the boundary are clamped. Suitable for tensile test simulations where boundary nodes act as grips.
- **periodic** — the domain wraps in-plane. Suitable for bulk network simulations where edge effects should be avoided.

```matlab
net.domain.boundary = 'periodic';
```

---

### `seed`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | — |

The random number generator seed. Only used when `imanualseed = true`. Setting a fixed seed ensures reproducible network generation.

:::tip Reproducibility
Always set a manual seed when generating networks for publication or benchmarking. This guarantees the exact same topology can be recreated.
:::

```matlab
net.flags.imanualseed = true;
net.domain.seed = 42;
```

---

### `write_location`

| Type | Args | Default |
|------|------|---------|
| `string` | any valid path | `'./networks'` |

The directory path where all output files are written. The directory must exist before calling `genNetwork()`.

```matlab
net.domain.write_location = './output/';
```

---

### `lammps_data_file`

| Type | Args | Default |
|------|------|---------|
| `string` | any string | `'PolyNetwork'` |

Prefix for the LAMMPS data file name. The full filename will be `<prefix>.data`.

```matlab
net.domain.lammps_data_file = 'my_network';
```

---

### `lammps_viz_file`

| Type | Args | Default |
|------|------|---------|
| `string` | any string | `'PolyVisual'` |

Prefix for the LAMMPS visualization file name.

```matlab
net.domain.lammps_viz_file = 'my_network_viz';
```

---

### `bond_table_file`

| Type | Args | Default |
|------|------|---------|
| `string` | any string | `'bond'` |

Prefix for the bond potential table file. Only written when `ipotential = true`.

```matlab
net.domain.bond_table_file = 'my_bond_table';
```

---

### `smp_number`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `1` |

Sample number appended to output filenames when `savemode = true`. Useful for generating multiple network realizations in a loop.

```matlab
for i = 1:10
    net.domain.smp_number = i;
    net.generateNetwork();
end
```
