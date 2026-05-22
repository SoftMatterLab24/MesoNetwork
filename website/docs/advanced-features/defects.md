---
custom_edit_url: null
sidebar_position: 1
---

# Defects

NetworkGen can introduce voids and damage into the network to simulate pre-existing defects or heterogeneous microstructures. Defect generation is activated by setting `idefect = true`.

---

### `defect.density_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'count'` \| `'area_frac'` | `'area_frac'` |

Controls how the number of voids is specified.

- **count** — specify an exact number of voids via `N_voids`
- **area_fraction** — specify what fraction of the domain area should be void via `void_area_frac`

```matlab
net.defect.density_mode = 'area_frac';
```

---

### `defect.n_voids`

| Type | Args | Default |
|------|------|---------|
| `int` | [0, ∞) | `0` |

Number of voids to introduce. Only used when `density_mode = 'count'`.

```matlab
net.defect.n_voids = 10;
```

---

### `defect.void_area_frac`

| Type | Args | Default |
|------|------|---------|
| `double` | [0, 1] | `0.75` |

Target fraction of domain area occupied by voids. Only used when `density_mode = 'area_frac'`.

```matlab
net.defect.void_area_frac = 0.1; % 10% void area
```

---

### `defect.size_dist`

| Type | Args | Default |
|------|------|---------|
| `string` | `'fixed'` \| `'gaussian'` \| `'exponential'` | `'gaussian'` |

Distribution from which void radii are drawn.

```matlab
net.defect.size_dist = 'gaussian';
```

---

### `defect.radius_mean`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `12` |

Mean void radius in units of `b`.

```matlab
net.defect.radius_mean = 2.0;
```

---

### `defect.radius_std`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `4` |

Standard deviation of the void radius distribution. Only used when `size_dist = 'gaussian'`.

```matlab
net.defect.radius_std = 0.5;
```

---

### `defect.radius_min`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, `radius_max`) | `2` |

Minimum allowable void radius.

```matlab
net.defect.radius_min = 0.5;
```

---

### `defect.radius_max`

| Type | Args | Default |
|------|------|---------|
| `double` | (`radius_min`, ∞) | `30` |

Maximum allowable void radius.

```matlab
net.defect.radius_max = 5.0;
```

---

### `defect.shape_roughness`

| Type | Args | Default |
|------|------|---------|
| `double` | [0, 1] | `0.3` |

Controls the irregularity of void boundaries. `0` = perfectly circular voids, `1` = maximally rough boundaries.

```matlab
net.defect.shape_roughness = 0.3;
```

---

### `defect.shape_n_modes`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `2` |

Number of Fourier modes used to generate rough void boundaries. Higher values produce finer-scale roughness.

```matlab
net.defect.shape_n_modes = 8;
```

---

### `defect.void_overlap`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'false'` |

When `true`, allows voids to overlap each other. When `false`, voids are placed without overlap.

```matlab
net.defect.void_overlap = false;
```

---

### `defect.center_distribution`

| Type | Args | Default |
|------|------|---------|
| `string` | `'random'` \| `'uniform'` \| `'clustered'` \| `'manual'` | `'clustered'` |

Controls the spatial distribution of void centers.

- **random** — void centers placed randomly
- **uniform** — void centers distributed evenly across the domain
- **clustered** — void centers grouped around parent locations
- **manual** — void centers specified explicitly via `manual_centers` (only works with `density_mode = 'count'`)

```matlab
net.defect.center_distribution = 'clustered';
```

---

### `defect.manual_centers`

| Type | Args | Default |
|------|------|---------|
| `double array` | [N × 2] matrix | `[]` |

Explicit coordinates of void centers when using `center_distribution = 'manual'`. Must be an N×2 array where N equals the number of voids (`n_voids`) and each row contains [x, y] coordinates.

Only used when `center_distribution = 'manual'` AND `density_mode = 'count'`.

Coordinates outside the domain are automatically clamped to the domain boundaries.

```matlab
net.defect.density_mode = 'count';
net.defect.n_voids = 3;
net.defect.center_distribution = 'manual';
net.defect.manual_centers = [25.0, 30.0; 50.0, 50.0; 75.0, 35.0];
```

---

### `defect.n_cluster_parents`

| Type | Args | Default |
|------|------|---------|
| `double` | [1, ∞) | `2` |

Number of cluster parent locations when `center_distribution = 'clustered'`.

```matlab
net.defect.n_cluster_parents = 3;
```

---

### `defect.cluster_spread`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `10` |

Standard deviation of void center positions around each cluster parent, in units of `b`.

```matlab
net.defect.cluster_spread = 2.0;
```

---

### `defect.margin_frac`

| Type | Args | Default |
|------|------|---------|
| `double` | [0, 0.5) | `0.15` |

Fraction of the domain size kept clear of void centers near the boundary. Prevents voids from straddling the domain edge.

```matlab
net.defect.margin_frac = 0.05;
```

---

### `defect.prune_isolated`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, removes isolated nodes (nodes with no remaining bonds) left behind after void creation.

```matlab
net.defect.prune_isolated = true;
```

---

### `defect.sparse_network`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'false'` |

When `true`, applies additional bond removal around void boundaries to produce a sparser network in those regions.

```matlab
net.defect.sparse_network = false;
```

---

### `defect.bridge_width`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `1` |

Minimum width of material bridges between closely spaced voids. Prevents voids from merging into a single large void.

```matlab
net.defect.bridge_width = 1.0;
```

---

### `defect.wall_thickness`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `18` |

Thickness of the solid wall region preserved at the domain boundary, in units of `b`.

```matlab
net.defect.wall_thickness = 2.0;
```

---

### `defect.clamp_thickness`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `0.12` |

Thickness of the clamped (fixed) region at the domain boundary. Used in conjunction with `boundary = 'fixed'` to define grip regions in tensile test simulations.

```matlab
net.defect.clamp_thickness = 1.5;
```

---

## Examples

### Random void distribution

```matlab
net.flags.idefect = true;
net.defect.density_mode = 'count';
net.defect.n_voids = 5;
net.defect.size_dist = 'gaussian';
net.defect.radius_mean = 2.0;
net.defect.radius_std = 0.5;
net.defect.radius_min = 0.5;
net.defect.radius_max = 4.0;
net.defect.center_distribution = 'random';
net.defect.void_overlap = false;
net.defect.prune_isolated = true;
```

### Manual void placement

```matlab
net.flags.idefect = true;
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
net.defect.prune_isolated = true;
```
