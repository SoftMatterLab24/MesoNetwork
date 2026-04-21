---
custom_edit_url: null
sidebar_position: 5
---

# Bimodal

The `bimodal` assignment mode draws values from a mixture of two distributions, producing two distinct populations. Used when `strand_typology = 'bimodal'` or when a perbond property is configured with a bimodal distribution.

---

## Core settings

These settings must always be configured when using the bimodal mode.

### `bimodal.mean_1`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, `mean_2`) | `35` |

Mean of the first (typically shorter/lower) distribution.

```matlab
net.architecture.strand_typology.bimodal.mean_1 = 3;
```

---

### `bimodal.mean_2`

| Type | Args | Default |
|------|------|---------|
| `double` | (`mean_1`, ∞) | `60` |

Mean of the second (typically longer/higher) distribution.

```matlab
net.architecture.strand_typology.bimodal.mean_2 = 10;
```

---

### `bimodal.std_1`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `10` |

Standard deviation of the first distribution.

```matlab
net.architecture.strand_typology.bimodal.std_1 = 1.0;
```

---

### `bimodal.std_2`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `5` |

Standard deviation of the second distribution.

```matlab
net.architecture.strand_typology.bimodal.std_2 = 2.0;
```

---

### `bimodal.method`

| Type | Args | Default |
|------|------|---------|
| `string` | `'single'` \| `'geom'` \| `'gaussian'` | `'gaussian'` |

The method used to generate the bimodal distribution.

- **single** — assigns the mean value directly with no spread. All bonds in each mode receive exactly `mean_1` or `mean_2` — effectively a deterministic two-value assignment with no distribution.
- **geom** — uses geometric distributions for each mode
- **gaussian** — uses Gaussian distributions for each mode, parameterised by `mean` and `std`

```matlab
net.architecture.strand_typology.bimodal.method = 'gaussian';
```

---

### `bimodal.height_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'prob'` \| `'count'` | `'prob'` |

Controls how the relative contribution of each distribution is specified.

- **prob** — specified as a probability (fraction of bonds from each mode) via `height_prob`
- **count** — specified as a fixed count via `height_count`

```matlab
net.architecture.strand_typology.bimodal.height_mode = 'prob';
```

---

### `bimodal.height_prob`

| Type | Args | Default |
|------|------|---------|
| `double` | [0, 1] | `1.0` |

Probability that a value is drawn from the first distribution. Only used when `height_mode = 'prob'`. The remaining fraction `(1 - height_prob)` is drawn from the second distribution.

```matlab
net.architecture.strand_typology.bimodal.height_prob = 0.4;
```

---

### `bimodal.height_count`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `2` |

Fixed number of values drawn from the first distribution. Only used when `height_mode = 'count'`.

```matlab
net.architecture.strand_typology.bimodal.height_count = 50;
```

---

### `bimodal.long_first`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'false'` |

When `true`, assigns values from the longer/higher distribution first during bond assignment.

```matlab
net.architecture.strand_typology.bimodal.long_first = true;
```

---

## Advanced settings

:::note Scope
The following settings are specifically relevant when using `bimodal` mode for `strand_typology` or `perbond.kuhn`, particularly for double network configurations and precise prestretch control. They are not needed for standard bimodal distributions.
:::

### `bimodal.stdR_1`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `3` |

Standard deviation of the bond end-to-end length for the first mode. Controls the spread of physical bond extension around the mean end-to-end length derived from `mean_1`. Used in conjunction with `lam_1` to set the target prestretch distribution.

```matlab
net.architecture.strand_typology.bimodal.stdR_1 = 0.5;
```

---

### `bimodal.stdR_2`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `10` |

Standard deviation of the bond end-to-end length for the second mode.

```matlab
net.architecture.strand_typology.bimodal.stdR_2 = 0.8;
```

---

### `bimodal.lam_1`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, 1) | `0.2` |

Target mean prestretch for bonds in the first mode. Defines the desired mean end-to-end extension relative to the contour length of strands in this mode. Works alongside `auto_1_flag` and `stdR_1` to fully specify the prestretch distribution.

```matlab
net.architecture.strand_typology.bimodal.lam_1 = 0.5;
```

---

### `bimodal.lam_2`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, 1) | `0.5` |

Target mean prestretch for bonds in the second mode.

```matlab
net.architecture.strand_typology.bimodal.lam_2 = 0.3;
```

---

### `bimodal.auto_1_flag`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'false'` |

When `true`, overrides the user-specified `mean_1` and instead automatically adjusts the Kuhn count for mode 1 to ensure the mean bond prestretch matches `lam_1`. This is useful when you want to control the network's mechanical state (prestretch) directly rather than specifying Kuhn counts explicitly.

:::tip Controlling prestretch deviation
Once `auto_1_flag` ensures the mean prestretch is accurate, you can control the *spread* of the prestretch distribution by setting `bin_window_method = 'manual'` and `manual_deviation_type = 'both'`, then adjusting `stdR_1` to set the end-to-end length deviation.
:::

```matlab
net.architecture.strand_typology.bimodal.auto_1_flag = true;
net.architecture.strand_typology.bimodal.lam_1 = 0.5;
```

---

### `bimodal.auto_2_flag`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'false'` |

When `true`, overrides `mean_2` and automatically adjusts the Kuhn count for mode 2 to ensure the mean bond prestretch matches `lam_2`.

```matlab
net.architecture.strand_typology.bimodal.auto_2_flag = true;
net.architecture.strand_typology.bimodal.lam_2 = 0.3;
```

---

### `bimodal.double_net_flag`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'false'` |

When `true`, generates an interpenetrating double network structure where each distribution mode corresponds to a distinct sub-network.

```matlab
net.architecture.strand_typology.bimodal.double_network_flag = true;
```

---

### `bimodal.alpha`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, 1) | `3.0` |

Mesh ratio for the double network. Defines the size of the second sub-network relative to the first — a smaller `alpha` produces a finer second network embedded within the coarser first. Only used when `double_net_flag = true`.

```matlab
net.architecture.strand_typology.bimodal.alpha = 0.5;
```

---

### `bimodal.bin_window_method`

| Type | Args | Default |
|------|------|---------|
| `string` | `'manual'` \| `'adaptive'` | `'manual'` |

Controls whether the bin window used to sort bonds into each distribution mode is computed automatically or defined manually.

- **adaptive** — the code automatically determines the standard deviation and bin boundaries based on the distribution parameters
- **manual** — the user controls the bin boundaries via `manual_deviation_type`, allowing explicit control over the prestretch deviation

```matlab
net.architecture.strand_typology.bimodal.bin_window_method = 'manual';
```

---

### `bimodal.manual_dev_type`

| Type | Args | Default |
|------|------|---------|
| `string` | `'mixed'` \| `'kuhn'` \| `'both'` | `'mixed'` |

Only used when `bin_window_method = 'manual'`. Specifies which deviation parameters are used to define the bin window.

- **kuhn** — the bin window is defined using only the Kuhn length deviation, controlling the spread in strand stiffness
- **both** — the bin window uses both the Kuhn deviation and the end-to-end length deviation (`stdR_1`/`stdR_2`), giving independent control over both the stiffness and the physical extension spread

:::tip Full prestretch control workflow
For complete control over prestretch statistics:
1. Set `auto_1_flag = true` — ensures mean prestretch matches `lam_1`
2. Set `bin_window_method = 'manual'` and `manual_deviation_type = 'both'`
3. Adjust `stdR_1` to set the desired spread in bond end-to-end lengths
:::

```matlab
net.architecture.strand_typology.bimodal.manual_dev_type = 'both';
```

---

## Example

```matlab
net.architecture.strand_typology.mode = 'bimodal';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.bimodal.method = 'gaussian';
net.architecture.strand_typology.bimodal.mean_1 = 3;
net.architecture.strand_typology.bimodal.mean_2 = 10;
net.architecture.strand_typology.bimodal.std_1 = 1.0;
net.architecture.strand_typology.bimodal.std_2 = 2.0;
net.architecture.strand_typology.bimodal.height_mode = 'prob';
net.architecture.strand_typology.bimodal.height_prob = 0.5;
```

**With prestretch control:**
```matlab
net.architecture.strand_typology.bimodal.auto_1_flag = true;
net.architecture.strand_typology.bimodal.auto_2_flag = true;
net.architecture.strand_typology.bimodal.lam_1 = 0.5;
net.architecture.strand_typology.bimodal.lam_2 = 0.3;
net.architecture.strand_typology.bimodal.bin_window_method = 'manual';
net.architecture.strand_typology.bimodal.manual_dev_type = 'both';
net.architecture.strand_typology.bimodal.stdR_1 = 0.5;
net.architecture.strand_typology.bimodal.stdR_2 = 0.8;
```
