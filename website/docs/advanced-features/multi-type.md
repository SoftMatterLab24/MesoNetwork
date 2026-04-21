---
custom_edit_url: null
sidebar_position: 2
---

# Multi-type Networks

NetworkGen supports exported networks with multiple atom and bond types, allowing heterogeneous LAMMPS labels on the final cleaned network.

Multi-type assignment runs after topology generation, defects, and cleanup. It labels the surviving atoms and bonds that are written to the LAMMPS data file; it does not change the earlier bond-generation algorithms.

If the requested atom or bond fractions cannot be matched exactly under the supplied connectivity rules, NetworkGen uses the closest feasible assignment it can find and reports the realized counts in the log.

---

### `types.enabled`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `true` \| `false` | `false` |

Turns post-cleanup multi-type labeling on or off.

```matlab
net.architecture.types.enabled = true;
```

---

### `types.natom_type`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `1` |

Number of distinct atom types in the network.

```matlab
net.architecture.types.natom_type = 2;
```

---

### `types.nbond_type`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `1` |

Number of distinct bond types in the network.

```matlab
net.architecture.types.nbond_type = 2;
```

---

### `types.atype_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'fixed'` \| `'frac'` | `'frac'` |

Controls how atom type counts are specified.

- **fixed** — specify exact counts via `atom_count`
- **frac** — specify fractions via `atom_frac`

```matlab
net.architecture.types.atype_mode = 'frac';
```

---

### `types.btype_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'fixed'` \| `'frac'` | `'frac'` |

Controls how bond type counts are specified.

- **fixed** — specify exact counts via `bond_count`
- **frac** — specify fractions via `bond_frac`

```matlab
net.architecture.types.btype_mode = 'frac';
```

---

### `types.atom_count`

| Type | Args | Default |
|------|------|---------|
| `int array` | size: [1 x `natom_type`] | — |

Array of atom counts per type. Only used when `atype_mode = 'fixed'`.

```matlab
net.architecture.types.atom_count = [100, 50];
```

---

### `types.bond_count`

| Type | Args | Default |
|------|------|---------|
| `int array` | size: [1 x `nbond_type`] | — |

Array of bond counts per type. Only used when `btype_mode = 'fixed'`.

```matlab
net.architecture.types.bond_count = [200, 80];
```

---

### `types.atom_frac`

| Type | Args | Default |
|------|------|---------|
| `double array` | size: [1 x `natom_type`], sum = 1 | — |

Array of atom type fractions. Must sum to 1. Only used when `atype_mode = 'frac'`.

```matlab
net.architecture.types.atom_frac = [0.7, 0.3];
```

---

### `types.bond_frac`

| Type | Args | Default |
|------|------|---------|
| `double array` | size: [1 x `nbond_type`], sum = 1 | — |

Array of bond type fractions. Must sum to 1. Only used when `btype_mode = 'frac'`.

```matlab
net.architecture.types.bond_frac = [0.6, 0.4];
```

---

### `types.connectivity`

| Type | Args | Default |
|------|------|---------|
| `int array` | size: [N x 4] | `[]` (empty) |

Connectivity rule table for exported bond labels. Each row is:

```text
[atom_type_A, atom_type_B, bond_type, allowed]
```

Atom-type order is symmetric, so a rule for `[1, 2, 1, 0]` also applies to `[2, 1, 1, 0]`.

Any `(atom type, atom type, bond type)` triple that is not listed is assumed to be allowed.

:::note Default behavior
The default mode is permissive: anything is allowed unless you explicitly constrain it.

For example, if you specify only `[1, 2, 1, 0]`, then type-1 bonds are forbidden between atom types 1 and 2, but type-2 bonds between those same atom types are still allowed because that triple was not constrained.
:::

```matlab
% Bond type 1 cannot connect atom types 1 and 2
% Bond type 2 is still allowed for atom types 1 and 2
net.architecture.types.connectivity = [ ...
    1 2 1 0; ...
    1 2 2 1  ...
];
```

:::tip Legacy rule tables
Legacy 3-column rule tables of the form `[atom_type_A, atom_type_B, allowed]` are still accepted for backward compatibility. In that case the same rule is applied to every bond type.
:::

---

### `types.atype_sel_method`

| Type | Args | Default |
|------|------|---------|
| `string array` | `'random'`, size: [1 x `natom_type`] | `'random'` |

Method used to select which atoms are assigned each type. Currently only `'random'` is supported.

---

### `types.btype_sel_method`

| Type | Args | Default |
|------|------|---------|
| `string array` | `'random'`, size: [1 x `nbond_type`] | `'random'` |

Method used to select which bonds are assigned each type. Currently only `'random'` is supported.

---

## Example

```matlab
net.architecture.types.enabled = true;
net.architecture.types.natom_type = 2;
net.architecture.types.nbond_type = 2;
net.architecture.types.atype_mode = 'frac';
net.architecture.types.btype_mode = 'frac';
net.architecture.types.atom_frac = [0.7, 0.3];
net.architecture.types.bond_frac = [0.6, 0.4];

% bond type 1 forbidden for atom-type pair (1,2)
% bond type 2 still allowed for the same pair
net.architecture.types.connectivity = [ ...
    1 2 1 0; ...
    1 2 2 1  ...
];
```
