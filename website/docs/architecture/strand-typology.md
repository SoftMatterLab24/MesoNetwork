---
sidebar_position: 2
---

# Strand Typology

Strand typology sets the target bond length distribution across the network. It is one of the most important settings in NetworkGen as it directly determines the mechanical heterogeneity of the generated network.

---

### `strand_typology`

| Type | Args | Default |
|------|------|---------|
| `string` | `'mono'` \| `'polydisperse'` \| `'poly'` \| `'bimodal'` | `'mono'` |

Sets the target bond length distribution type for the network.

- **mono** — attempts to make all bonds the same length. Due to the geometric constraints of random node placement, exact uniformity is not guaranteed, but the generator minimizes length variation as much as possible.
- **polydisperse** — bond lengths are drawn from a continuous distribution. This is the preferred website label; the package also accepts the legacy alias `poly`. Configure the `poly` sub-settings in [Assignment Modes](../assignment-modes/overview).
- **bimodal** — bond lengths are drawn from a mixture of two distributions, producing two distinct strand populations. Requires configuring the `bimodal` sub-settings in [Assignment Modes](../assignment-modes/overview).

:::note Linked settings
Setting `strand_typology` to `polydisperse` (or the legacy alias `poly`) or `bimodal` activates the corresponding **Assignment Mode** sub-settings. See [Assignment Modes](../assignment-modes/overview) for full configuration options.
:::

:::tip Coupling with perbond
Any perbond property can be configured to copy the same distribution used here. For example, when `perbond.kuhn` is set to `auto`, it copies this distribution exactly — ensuring bond stiffness stays consistent with bond length. See [Perbond](../advanced-features/potentials-perbond).
:::

```matlab
% Mono — generator attempts uniform bond lengths
net.architecture.strand_typology.mode = 'mono';

% Polydisperse — configure assignment_mode.poly
net.architecture.strand_typology.mode = 'polydisperse';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.poly.method = 'pmf';
net.architecture.strand_typology.poly.pmf_mean = 5;

% Bimodal — configure assignment_mode.bimodal
net.architecture.strand_typology.mode = 'bimodal';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.bimodal.mean_1 = 3;
net.architecture.strand_typology.bimodal.mean_2 = 8;
```
