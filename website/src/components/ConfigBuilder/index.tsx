import React, { useState, useCallback } from 'react';
import styles from './styles.module.css';

type ConnectivityRule = {
  atomTypeA: number;
  atomTypeB: number;
  bondType: number;
  allowed: boolean;
};

const SECTION_COLORS = {
  domain:       { accent: '#378ADD', bg: '#E6F1FB' },
  architecture: { accent: '#1D9E75', bg: '#E1F5EE' },
  typology:     { accent: '#7F77DD', bg: '#EEEDFE' },
  perbond:      { accent: '#D85A30', bg: '#FAECE7' },
  defect:       { accent: '#BA7517', bg: '#FAEEDA' },
  multitype:    { accent: '#D4537E', bg: '#FBEAF0' },
  potential:    { accent: '#3B6D11', bg: '#EAF3DE' },
  flags:        { accent: '#5F5E5A', bg: '#F1EFE8' },
};

const DEFAULT = {
  b: 1.6, Lx: 10, Ly: 10, boundary: 'fixed', seed: 12345,
  write_location: './networks', lammps_data_file: 'PolyNetwork',
  lammps_viz_file: 'PolyVisual', bond_table_file: 'bond', scale: 1,

  geometry: 'random', rho_atom: 0.0078, max_peratom_bond: 6, min_degree_keep: 2,
  lattice_spacing: 6, spacing_multiplier_mode: 'auto',
  spacing_multiplier: 1.2, lattice_disorder_level: 0,
  lattice_disorder_maxfrac: 0.4, lattice_topo_disorder: false,
  lattice_max_del_per_node: 1, lattice_min_degree_keep: 5,

  typology_mode: 'mono',
  mono_value: 20,
  poly_method: 'pmf', poly_min_value: 1, poly_pmf_mean: 40, poly_pmf_min: 5,
  poly_pmf_max: 120, poly_rounding: 'round', poly_align: 'none',
  poly_range_method: 'rank', poly_target_min: 5, poly_target_max: 120,
  bimodal_method: 'gaussian', bimodal_mean1: 10, bimodal_mean2: 40,
  bimodal_std1: 2, bimodal_std2: 5, bimodal_height_mode: 'prob',
  bimodal_height_prob: 0.5, bimodal_height_count: 2, bimodal_long_first: true,
  bimodal_min_value: 1, bimodal_double_network_flag: false, bimodal_alpha: 3.0,
  bimodal_auto_1_flag: false, bimodal_auto_2_flag: false,
  bimodal_bin_window_method: 'manual', bimodal_manual_dev_type: 'mixed',
  bimodal_stdR_1: 3, bimodal_stdR_2: 10, bimodal_lam_1: 0.2, bimodal_lam_2: 0.5,

  kuhn_auto: true, kuhn_mode: 'mono', kuhn_mono_value: 20,
  kuhn_poly_method: 'pmf', kuhn_poly_min_value: 1, kuhn_poly_pmf_mean: 40,
  kuhn_poly_pmf_min: 20, kuhn_poly_pmf_max: 120, kuhn_poly_rounding: 'round',
  kuhn_poly_align: 'ascend', kuhn_poly_range_method: 'rank',
  kuhn_poly_target_min: 5, kuhn_poly_target_max: 120,
  kuhn_bimodal_method: 'gaussian', kuhn_bimodal_mean1: 10, kuhn_bimodal_mean2: 40,
  kuhn_bimodal_std1: 2, kuhn_bimodal_std2: 5, kuhn_bimodal_height_mode: 'prob',
  kuhn_bimodal_height_prob: 0.5, kuhn_bimodal_height_count: 2, kuhn_bimodal_long_first: true,
  kuhn_bimodal_min_value: 1, kuhn_bimodal_double_network_flag: false, kuhn_bimodal_alpha: 3.0,
  kuhn_bimodal_auto_1_flag: false, kuhn_bimodal_auto_2_flag: false,
  kuhn_bimodal_bin_window_method: 'manual', kuhn_bimodal_manual_dev_type: 'mixed',
  kuhn_bimodal_stdR_1: 3, kuhn_bimodal_stdR_2: 10, kuhn_bimodal_lam_1: 0.2, kuhn_bimodal_lam_2: 0.5,

  idefect: false, defect_density_mode: 'count', defect_n_voids: 5,
  defect_void_area_frac: 0.1, defect_size_dist: 'gaussian',
  defect_radius_mean: 12, defect_radius_std: 4,
  defect_radius_min: 2, defect_radius_max: 30,
  defect_shape_roughness: 0.3, defect_shape_n_modes: 2,
  defect_void_overlap: false, defect_center_dist: 'random',
  defect_n_cluster_parents: 2, defect_cluster_spread: 10,
  defect_margin_frac: 0.15, defect_prune_isolated: true,
  defect_sparse_network: false, defect_wall_thickness: 18,
  defect_clamp_thickness: 0.12, defect_bridge_width: 1,
  defect_thinning: false, defect_thinning_radius: 0,
  defect_thinning_target_frac: 0.4, defect_thinning_min_keep: 0.1,
  defect_bridging: false, defect_bridge_max_dist: 0,
  defect_bridge_void_thresh: 0.25, defect_bridge_perp_width: 0,
  defect_bridge_max_degree: 0, defect_bridge_max_bonds: 0,
  defect_bridge_min_spacing: 0,

  use_multitype: false, natom_type: 1, nbond_type: 1,
  atype_mode: 'frac', btype_mode: 'frac',
  atom_frac_values: [1],
  bond_frac_values: [1],
  atom_count_values: [1],
  bond_count_values: [1],
  connectivity_rules: [] as ConnectivityRule[],

  ipotential: false, pot_k_LD: 0.414, pot_N_rho: 100000,
  pot_rho_min: 0.0, pot_rho_max: 500,

  isave: true, iplot: true, ilog: true,
  idumpsettings: false, iversbose_settings: false,
  savemode: true, imanualseed: false,
  Nreplicates: 1,

  // Python-specific
  py_networkgen_path: '/path/to/NetworkGen',
};

function uniformFractions(n) {
  if (n <= 0) {
    return [];
  }

  if (n === 1) {
    return [1];
  }

  const vals = [];
  let running = 0;
  for (let i = 0; i < n - 1; i += 1) {
    const val = Number((1 / n).toFixed(6));
    vals.push(val);
    running += val;
  }
  vals.push(Number((1 - running).toFixed(6)));
  return vals;
}

function unitCounts(n) {
  return Array.from({ length: Math.max(0, n) }, () => 1);
}

function resizeNumericArray(values, count, fillValue) {
  const next = Array.isArray(values) ? values.slice(0, count) : [];

  while (next.length < count) {
    next.push(fillValue);
  }

  return next;
}

function sanitizeWeightArray(values, count) {
  return resizeNumericArray(values, count, 0)
    .map(value => Math.max(0, Number.isFinite(value) ? value : 0));
}

function sanitizeCountArray(values, count) {
  return resizeNumericArray(values, count, 1)
    .map(value => Math.max(0, Math.round(Number.isFinite(value) ? value : 0)));
}

function normalizeWeightArray(values, count) {
  const sanitized = sanitizeWeightArray(values, count);
  const total = sanitized.reduce((sum, value) => sum + value, 0);

  if (total <= 0) {
    return uniformFractions(count);
  }

  return sanitized.map(value => Number((value / total).toFixed(6)));
}

function sanitizeConnectivityRules(rules, natomType, nbondType) {
  return (Array.isArray(rules) ? rules : []).map(rule => ({
    atomTypeA: Math.min(Math.max(1, Math.round(rule.atomTypeA || 1)), natomType),
    atomTypeB: Math.min(Math.max(1, Math.round(rule.atomTypeB || 1)), natomType),
    bondType: Math.min(Math.max(1, Math.round(rule.bondType || 1)), nbondType),
    allowed: !!rule.allowed,
  }));
}

function connectivityRowsToMatlabLines(rules, indent = '') {
  if (!rules.length) {
    return [`${indent}net.architecture.types.connectivity = [];  %% leave empty to allow all combinations`];
  }

  const lines = [`${indent}net.architecture.types.connectivity = [ ...`];

  rules.forEach((rule, index) => {
    const suffix = index === rules.length - 1 ? '  ...' : '; ...';
    lines.push(
      `${indent}    ${rule.atomTypeA} ${rule.atomTypeB} ${rule.bondType} ${rule.allowed ? 1 : 0}${suffix}`
    );
  });

  lines.push(`${indent}];`);
  return lines;
}

function pushPolyAssignmentConfig(lines, basePath, values) {
  lines.push(`${basePath}.poly.method = '${values.method}';`);
  lines.push(`${basePath}.poly.min_value = ${values.minValue};`);
  lines.push(`${basePath}.poly.rounding = '${values.rounding}';`);
  lines.push(`${basePath}.poly.align_to_length = '${values.alignToLength}';`);

  if (values.method === 'range') {
    lines.push(`${basePath}.poly.range_method = '${values.rangeMethod}';`);
    lines.push(`${basePath}.poly.target_min = ${values.targetMin};`);
    lines.push(`${basePath}.poly.target_max = ${values.targetMax};`);
  }

  if (values.method === 'pmf') {
    lines.push(`${basePath}.poly.pmf_mean = ${values.pmfMean};`);
    lines.push(`${basePath}.poly.pmf_min = ${values.pmfMin};`);
    lines.push(`${basePath}.poly.pmf_max = ${values.pmfMax};`);
  }
}

function pushBimodalAssignmentConfig(lines, basePath, values, formatBool = value => `${value}`) {
  lines.push(`${basePath}.bimodal.method = '${values.method}';`);
  lines.push(`${basePath}.bimodal.mean_1 = ${values.mean1};`);
  lines.push(`${basePath}.bimodal.mean_2 = ${values.mean2};`);

  if (values.method !== 'single') {
    lines.push(`${basePath}.bimodal.std_1 = ${values.std1};`);
    lines.push(`${basePath}.bimodal.std_2 = ${values.std2};`);
  }

  lines.push(`${basePath}.bimodal.height_mode = '${values.heightMode}';`);
  if (values.heightMode === 'prob') {
    lines.push(`${basePath}.bimodal.height_prob = ${values.heightProb};`);
  } else {
    lines.push(`${basePath}.bimodal.height_count = ${values.heightCount};`);
  }

  lines.push(`${basePath}.bimodal.long_first = ${formatBool(values.longFirst)};`);
  lines.push(`${basePath}.bimodal.min_value = ${values.minValue};`);
  lines.push(`${basePath}.bimodal.double_network_flag = ${formatBool(values.doubleNetworkFlag)};`);
  if (values.doubleNetworkFlag) {
    lines.push(`${basePath}.bimodal.alpha = ${values.alpha};`);
  }
  lines.push(`${basePath}.bimodal.auto_1_flag = ${formatBool(values.auto1Flag)};`);
  if (values.auto1Flag) {
    lines.push(`${basePath}.bimodal.lam_1 = ${values.lam1};`);
  }
  lines.push(`${basePath}.bimodal.auto_2_flag = ${formatBool(values.auto2Flag)};`);
  if (values.auto2Flag) {
    lines.push(`${basePath}.bimodal.lam_2 = ${values.lam2};`);
  }
  lines.push(`${basePath}.bimodal.stdR_1 = ${values.stdR1};`);
  lines.push(`${basePath}.bimodal.stdR_2 = ${values.stdR2};`);
  lines.push(`${basePath}.bimodal.bin_window_method = '${values.binWindowMethod}';`);
  if (values.binWindowMethod === 'manual') {
    lines.push(`${basePath}.bimodal.manual_dev_type = '${values.manualDevType}';`);
  }
}

function matlabVector(values) {
  return `[${values.join(' ')}]`;
}

function Row({ label, hint, children }) {
  return (
    <div className={styles.row}>
      <label className={styles.rowLabel}>
        {label}
        {hint && <span className={styles.rowHint}>{hint}</span>}
      </label>
      <div className={styles.rowInput}>{children}</div>
    </div>
  );
}

function Section({ id, title, children, defaultOpen = false }) {
  const [open, setOpen] = useState(defaultOpen);
  const col = SECTION_COLORS[id] || SECTION_COLORS.flags;
  return (
    <div className={styles.section} style={{ borderLeftColor: col.accent }}>
      <button
        className={styles.sectionHead}
        style={{ background: col.bg + '66' }}
        onClick={() => setOpen(o => !o)}
      >
        <span className={styles.sectionDot} style={{ background: col.accent }} />
        <span className={styles.sectionTitle}>{title}</span>
        <span className={styles.sectionToggle}>{open ? '▲' : '▼'}</span>
      </button>
      {open && <div className={styles.sectionBody}>{children}</div>}
    </div>
  );
}

function Sub({ title, children }) {
  return (
    <div className={styles.sub}>
      <div className={styles.subLabel}>{title}</div>
      {children}
    </div>
  );
}

// ── Language toggle pill ───────────────────────────────────────────────────────
function LangToggle({ lang, setLang }) {
  const base: React.CSSProperties = {
    padding: '4px 14px',
    fontSize: 12,
    fontWeight: 500,
    border: 'none',
    cursor: 'pointer',
    transition: 'background 0.15s, color 0.15s',
    fontFamily: 'var(--ifm-font-family-base)',
  };
  const active: React.CSSProperties = {
    background: 'var(--ifm-color-primary)',
    color: '#fff',
  };
  const inactive: React.CSSProperties = {
    background: 'var(--ifm-color-emphasis-200)',
    color: 'var(--ifm-color-content-secondary)',
  };
  return (
    <div style={{
      display: 'inline-flex',
      border: '0.5px solid var(--ifm-color-emphasis-300)',
      borderRadius: 6,
      overflow: 'hidden',
    }}>
      <button
        style={{ ...base, ...(lang === 'matlab' ? active : inactive), borderRadius: '6px 0 0 6px' }}
        onClick={() => setLang('matlab')}
      >MATLAB</button>
      <button
        style={{ ...base, ...(lang === 'python' ? active : inactive), borderRadius: '0 6px 6px 0' }}
        onClick={() => setLang('python')}
      >Python</button>
    </div>
  );
}

export default function ConfigBuilder() {
  const [cfg, setCfg] = useState(DEFAULT);
  const [copied, setCopied] = useState(false);
  const [lang, setLang] = useState<'matlab' | 'python'>('matlab');

  const set = useCallback((key, val) => {
    setCfg(c => ({ ...c, [key]: val }));
  }, []);

  const setArrayValue = useCallback((key, index, rawValue, integer = false) => {
    setCfg(c => {
      const next = Array.isArray(c[key]) ? [...c[key]] : [];
      const value = Number.isFinite(rawValue) ? rawValue : 0;
      next[index] = integer ? Math.max(0, Math.round(value)) : Math.max(0, value);
      return { ...c, [key]: next };
    });
  }, []);

  const resetArrayValues = useCallback((key, values) => {
    setCfg(c => ({ ...c, [key]: values }));
  }, []);

  const normalizeArrayValues = useCallback((key, count) => {
    setCfg(c => ({ ...c, [key]: normalizeWeightArray(c[key], count) }));
  }, []);

  const setAtomTypeCount = useCallback((rawValue) => {
    const nextCount = Math.max(1, Math.round(rawValue || 1));

    setCfg(c => ({
      ...c,
      natom_type: nextCount,
      atom_frac_values: resizeNumericArray(c.atom_frac_values, nextCount, 0),
      atom_count_values: resizeNumericArray(c.atom_count_values, nextCount, 1),
      connectivity_rules: sanitizeConnectivityRules(c.connectivity_rules, nextCount, c.nbond_type),
    }));
  }, []);

  const setBondTypeCount = useCallback((rawValue) => {
    const nextCount = Math.max(1, Math.round(rawValue || 1));

    setCfg(c => ({
      ...c,
      nbond_type: nextCount,
      bond_frac_values: resizeNumericArray(c.bond_frac_values, nextCount, 0),
      bond_count_values: resizeNumericArray(c.bond_count_values, nextCount, 1),
      connectivity_rules: sanitizeConnectivityRules(c.connectivity_rules, c.natom_type, nextCount),
    }));
  }, []);

  const addConnectivityRule = useCallback(() => {
    setCfg(c => ({
      ...c,
      connectivity_rules: [
        ...sanitizeConnectivityRules(c.connectivity_rules, c.natom_type, c.nbond_type),
        {
          atomTypeA: 1,
          atomTypeB: Math.min(2, c.natom_type),
          bondType: 1,
          allowed: false,
        },
      ],
    }));
  }, []);

  const updateConnectivityRule = useCallback((index, field, value) => {
    setCfg(c => {
      const next = sanitizeConnectivityRules(c.connectivity_rules, c.natom_type, c.nbond_type);

      if (!next[index]) {
        return c;
      }

      next[index] = {
        ...next[index],
        [field]: field === 'allowed' ? !!value : Math.round(value || 1),
      };

      return { ...c, connectivity_rules: next };
    });
  }, []);

  const removeConnectivityRule = useCallback((index) => {
    setCfg(c => ({
      ...c,
      connectivity_rules: sanitizeConnectivityRules(c.connectivity_rules, c.natom_type, c.nbond_type)
        .filter((_, rowIndex) => rowIndex !== index),
    }));
  }, []);

  const sel = (key, opts) => (
    <select value={cfg[key]} onChange={e => set(key, e.target.value)}>
      {opts.map(([v, l]) => <option key={v} value={v}>{l ?? v}</option>)}
    </select>
  );

  const num = (key, min?, max?, step = 1) => (
    <input type="number" value={cfg[key]}
      min={min} max={max} step={step}
      onChange={e => set(key, parseFloat(e.target.value) || 0)} />
  );

  const txt = (key) => (
    <input type="text" value={cfg[key]}
      onChange={e => set(key, e.target.value)} />
  );

  const chk = (key) => (
    <input type="checkbox" checked={cfg[key]}
      onChange={e => set(key, e.target.checked)} />
  );

  const slide = (key, min, max, step = 0.05) => (
    <div style={{ display: 'flex', gap: 8, alignItems: 'center', flex: 1 }}>
      <input type="range" min={min} max={max} step={step} value={cfg[key]}
        onChange={e => set(key, parseFloat(e.target.value))}
        style={{ flex: 1 }} />
      <span className={styles.slideVal}>{Number(cfg[key]).toFixed(step < 1 ? 2 : 0)}</span>
    </div>
  );

  // ── MATLAB script generator ────────────────────────────────────────────────
  function generateMatlab() {
    const c = cfg;
    const t = c.typology_mode;
    const lines: string[] = [];
    const atomTargets = c.atype_mode === 'fixed'
      ? matlabVector(sanitizeCountArray(c.atom_count_values, c.natom_type))
      : matlabVector(sanitizeWeightArray(c.atom_frac_values, c.natom_type));
    const bondTargets = c.btype_mode === 'fixed'
      ? matlabVector(sanitizeCountArray(c.bond_count_values, c.nbond_type))
      : matlabVector(sanitizeWeightArray(c.bond_frac_values, c.nbond_type));
    const connectivityLines = connectivityRowsToMatlabLines(
      sanitizeConnectivityRules(c.connectivity_rules, c.natom_type, c.nbond_type)
    );

    lines.push(`%% NetworkGen configuration script`);
    lines.push(`%% Generated by the NetworkGen config builder`);
    lines.push(`%% https://soft-matter-lab.github.io/networkgen`);
    lines.push(``);
    lines.push(`net = network();`);
    lines.push(`net.Nreplicates = ${c.Nreplicates};`);
    lines.push(``);
    lines.push(`%% ---- Domain ----`);
    lines.push(`net.domain.b                 = ${c.b};`);
    lines.push(`net.domain.Lx                = ${c.Lx};`);
    lines.push(`net.domain.Ly                = ${c.Ly};`);
    lines.push(`net.domain.scale             = ${c.scale};`);
    lines.push(`net.domain.boundary          = '${c.boundary}';`);
    if (c.imanualseed) lines.push(`net.domain.seed              = ${c.seed};`);
    lines.push(`net.domain.write_location    = '${c.write_location}';`);
    lines.push(`net.domain.lammps_data_file  = '${c.lammps_data_file}';`);
    lines.push(`net.domain.lammps_viz_file   = '${c.lammps_viz_file}';`);
    lines.push(`net.domain.bond_table_file   = '${c.bond_table_file}';`);
    lines.push(``);
    lines.push(`%% ---- Architecture ----`);
    lines.push(`net.architecture.geometry           = '${c.geometry}';`);
    lines.push(`net.architecture.rho_atom           = ${c.rho_atom};`);
    lines.push(`net.peratom.Max_peratom_bond        = ${c.max_peratom_bond};`);
    lines.push(`net.peratom.min_degree_keep         = ${c.min_degree_keep};`);
    if (c.geometry === 'hex_lattice') {
      lines.push(`net.architecture.lattice_spacing            = ${c.lattice_spacing};`);
      lines.push(`net.architecture.spacing_multiplier_mode    = '${c.spacing_multiplier_mode}';`);
      if (c.spacing_multiplier_mode === 'manual')
        lines.push(`net.architecture.spacing_multiplier         = ${c.spacing_multiplier};`);
      lines.push(`net.architecture.lattice_disorder_level     = ${c.lattice_disorder_level};`);
      lines.push(`net.architecture.lattice_disorder_maxfrac   = ${c.lattice_disorder_maxfrac};`);
      lines.push(`net.architecture.lattice_max_del_per_node   = ${c.lattice_max_del_per_node};`);
      lines.push(`net.architecture.lattice_min_degree_keep    = ${c.lattice_min_degree_keep};`);
    }
    lines.push(``);
    lines.push(`%% ---- Strand typology ----`);
    lines.push(`net.architecture.strand_typology.mode = '${t}';`);
    if (t === 'mono') {
      lines.push(`net.architecture.strand_typology.mono.value = ${c.mono_value};`);
    } else if (t === 'polydisperse') {
      pushPolyAssignmentConfig(lines, 'net.architecture.strand_typology', {
        method: c.poly_method,
        minValue: c.poly_min_value,
        rounding: c.poly_rounding,
        alignToLength: c.poly_align,
        rangeMethod: c.poly_range_method,
        targetMin: c.poly_target_min,
        targetMax: c.poly_target_max,
        pmfMean: c.poly_pmf_mean,
        pmfMin: c.poly_pmf_min,
        pmfMax: c.poly_pmf_max,
      });
    } else if (t === 'bimodal') {
      pushBimodalAssignmentConfig(lines, 'net.architecture.strand_typology', {
        method: c.bimodal_method,
        mean1: c.bimodal_mean1,
        mean2: c.bimodal_mean2,
        std1: c.bimodal_std1,
        std2: c.bimodal_std2,
        heightMode: c.bimodal_height_mode,
        heightProb: c.bimodal_height_prob,
        heightCount: c.bimodal_height_count,
        longFirst: c.bimodal_long_first,
        minValue: c.bimodal_min_value,
        doubleNetworkFlag: c.bimodal_double_network_flag,
        alpha: c.bimodal_alpha,
        auto1Flag: c.bimodal_auto_1_flag,
        auto2Flag: c.bimodal_auto_2_flag,
        lam1: c.bimodal_lam_1,
        lam2: c.bimodal_lam_2,
        stdR1: c.bimodal_stdR_1,
        stdR2: c.bimodal_stdR_2,
        binWindowMethod: c.bimodal_bin_window_method,
        manualDevType: c.bimodal_manual_dev_type,
      });
    }
    lines.push(``);
    lines.push(`%% ---- Perbond ----`);
    lines.push(`net.perbond.kuhn.auto = ${c.kuhn_auto};`);
    if (!c.kuhn_auto) {
      lines.push(`net.perbond.kuhn.mode = '${c.kuhn_mode}';`);
      if (c.kuhn_mode === 'mono')
        lines.push(`net.perbond.kuhn.mono.value = ${c.kuhn_mono_value};`);
      else if (c.kuhn_mode === 'polydisperse') {
        pushPolyAssignmentConfig(lines, 'net.perbond.kuhn', {
          method: c.kuhn_poly_method,
          minValue: c.kuhn_poly_min_value,
          rounding: c.kuhn_poly_rounding,
          alignToLength: c.kuhn_poly_align,
          rangeMethod: c.kuhn_poly_range_method,
          targetMin: c.kuhn_poly_target_min,
          targetMax: c.kuhn_poly_target_max,
          pmfMean: c.kuhn_poly_pmf_mean,
          pmfMin: c.kuhn_poly_pmf_min,
          pmfMax: c.kuhn_poly_pmf_max,
        });
      } else if (c.kuhn_mode === 'bimodal') {
        pushBimodalAssignmentConfig(lines, 'net.perbond.kuhn', {
          method: c.kuhn_bimodal_method,
          mean1: c.kuhn_bimodal_mean1,
          mean2: c.kuhn_bimodal_mean2,
          std1: c.kuhn_bimodal_std1,
          std2: c.kuhn_bimodal_std2,
          heightMode: c.kuhn_bimodal_height_mode,
          heightProb: c.kuhn_bimodal_height_prob,
          heightCount: c.kuhn_bimodal_height_count,
          longFirst: c.kuhn_bimodal_long_first,
          minValue: c.kuhn_bimodal_min_value,
          doubleNetworkFlag: c.kuhn_bimodal_double_network_flag,
          alpha: c.kuhn_bimodal_alpha,
          auto1Flag: c.kuhn_bimodal_auto_1_flag,
          auto2Flag: c.kuhn_bimodal_auto_2_flag,
          lam1: c.kuhn_bimodal_lam_1,
          lam2: c.kuhn_bimodal_lam_2,
          stdR1: c.kuhn_bimodal_stdR_1,
          stdR2: c.kuhn_bimodal_stdR_2,
          binWindowMethod: c.kuhn_bimodal_bin_window_method,
          manualDevType: c.kuhn_bimodal_manual_dev_type,
        });
      }
    }
    if (c.idefect) {
      lines.push(``);
      lines.push(`%% ---- Defects ----`);
      lines.push(`net.defect.density_mode       = '${c.defect_density_mode}';`);
      if (c.defect_density_mode === 'count')
        lines.push(`net.defect.n_voids            = ${c.defect_n_voids};`);
      else
        lines.push(`net.defect.void_area_frac     = ${c.defect_void_area_frac};`);
      lines.push(`net.defect.size_dist          = '${c.defect_size_dist}';`);
      lines.push(`net.defect.radius_mean        = ${c.defect_radius_mean};`);
      if (c.defect_size_dist !== 'fixed') {
        lines.push(`net.defect.radius_std         = ${c.defect_radius_std};`);
        lines.push(`net.defect.radius_min         = ${c.defect_radius_min};`);
        lines.push(`net.defect.radius_max         = ${c.defect_radius_max};`);
      }
      lines.push(`net.defect.shape_roughness    = ${c.defect_shape_roughness};`);
      lines.push(`net.defect.shape_n_modes      = ${c.defect_shape_n_modes};`);
      lines.push(`net.defect.void_overlap       = ${c.defect_void_overlap};`);
      lines.push(`net.defect.center_distribution = '${c.defect_center_dist}';`);
      if (c.defect_center_dist === 'clustered') {
        lines.push(`net.defect.n_cluster_parents  = ${c.defect_n_cluster_parents};`);
        lines.push(`net.defect.cluster_spread     = ${c.defect_cluster_spread};`);
      }
      lines.push(`net.defect.margin_frac        = ${c.defect_margin_frac};`);
      lines.push(`net.defect.prune_isolated     = ${c.defect_prune_isolated};`);
      lines.push(`net.defect.sparse_network     = ${c.defect_sparse_network};`);
      lines.push(`net.defect.wall_thickness     = ${c.defect_wall_thickness};`);
      lines.push(`net.defect.clamp_thickness    = ${c.defect_clamp_thickness};`);
      lines.push(`net.defect.bridge_width       = ${c.defect_bridge_width};`);
      if (c.defect_thinning) {
        lines.push(`net.defect.thinning           = true;`);
        lines.push(`net.defect.thinning_radius    = ${c.defect_thinning_radius};`);
        lines.push(`net.defect.thinning_target_frac = ${c.defect_thinning_target_frac};`);
        lines.push(`net.defect.thinning_min_keep  = ${c.defect_thinning_min_keep};`);
      }
      if (c.defect_bridging) {
        lines.push(`net.defect.bridging           = true;`);
        lines.push(`net.defect.bridge_max_dist    = ${c.defect_bridge_max_dist};`);
        lines.push(`net.defect.bridge_void_thresh = ${c.defect_bridge_void_thresh};`);
        lines.push(`net.defect.bridge_perp_width  = ${c.defect_bridge_perp_width};`);
        lines.push(`net.defect.bridge_max_degree  = ${c.defect_bridge_max_degree};`);
        lines.push(`net.defect.bridge_max_bonds   = ${c.defect_bridge_max_bonds};`);
        lines.push(`net.defect.bridge_min_spacing = ${c.defect_bridge_min_spacing};`);
      }
    }
    if (c.use_multitype) {
      lines.push(``);
      lines.push(`%% ---- Multi-type ----`);
      lines.push(`net.architecture.types.enabled     = true;`);
      lines.push(`net.architecture.types.natom_type  = ${c.natom_type};`);
      lines.push(`net.architecture.types.nbond_type  = ${c.nbond_type};`);
      lines.push(`net.architecture.types.atype_mode  = '${c.atype_mode}';`);
      lines.push(`net.architecture.types.btype_mode  = '${c.btype_mode}';`);
      if (c.atype_mode === 'fixed')
        lines.push(`net.architecture.types.atom_count  = ${atomTargets};`);
      else
        lines.push(`net.architecture.types.atom_frac   = ${atomTargets};`);
      if (c.btype_mode === 'fixed')
        lines.push(`net.architecture.types.bond_count  = ${bondTargets};`);
      else
        lines.push(`net.architecture.types.bond_frac   = ${bondTargets};`);
      lines.push(...connectivityLines);
    }
    if (c.ipotential) {
      lines.push(``);
      lines.push(`%% ---- Potential (pair local/density) ----`);
      lines.push(`net.pot.k_LD    = ${c.pot_k_LD};`);
      lines.push(`net.pot.N_rho   = ${c.pot_N_rho};`);
      lines.push(`net.pot.rho_min = ${c.pot_rho_min};`);
      lines.push(`net.pot.rho_max = ${c.pot_rho_max};`);
    }
    lines.push(``);
    lines.push(`%% ---- Flags ----`);
    lines.push(`net.flags.isave      = ${c.isave};`);
    lines.push(`net.flags.iplot      = ${c.iplot};`);
    lines.push(`net.flags.ilog       = ${c.ilog};`);
    lines.push(`net.flags.savemode   = ${c.savemode};`);
    lines.push(`net.flags.imanualseed = ${c.imanualseed};`);
    lines.push(`net.flags.idefect    = ${c.idefect};`);
    lines.push(`net.flags.ipotential = ${c.ipotential};`);
    lines.push(`net.flags.idumpsettings = ${c.idumpsettings};`);
    lines.push(`net.flags.iversbose_settings = ${c.iversbose_settings};`);
    lines.push(``);
    lines.push(`%% ---- Generate ----`);
    lines.push(`net.generateNetwork();`);

    return lines.join('\n');
  }

  // ── Python script generator ────────────────────────────────────────────────
  function generatePython() {
    const c = cfg;
    const t = c.typology_mode;
    const lines: string[] = [];
    const atomTargets = c.atype_mode === 'fixed'
      ? matlabVector(sanitizeCountArray(c.atom_count_values, c.natom_type))
      : matlabVector(sanitizeWeightArray(c.atom_frac_values, c.natom_type));
    const bondTargets = c.btype_mode === 'fixed'
      ? matlabVector(sanitizeCountArray(c.bond_count_values, c.nbond_type))
      : matlabVector(sanitizeWeightArray(c.bond_frac_values, c.nbond_type));
    const connectivityLines = connectivityRowsToMatlabLines(
      sanitizeConnectivityRules(c.connectivity_rules, c.natom_type, c.nbond_type),
      '    '
    );

    const boolStr = (v: boolean) => v ? 'true' : 'false';

    lines.push(`# NetworkGen configuration script (Python / oct2py)`);
    lines.push(`# Generated by the NetworkGen config builder`);
    lines.push(`# https://soft-matter-lab.github.io/networkgen`);
    lines.push(`#`);
    lines.push(`# Requirements: Python 3.10+, Octave 7.0+, oct2py, numpy`);
    lines.push(``);
    lines.push(`from oct2py import octave`);
    lines.push(``);
    lines.push(`octave.addpath(octave.genpath('${c.py_networkgen_path}'))`);
    lines.push(``);
    lines.push(`octave.eval("""`);

    // All the MATLAB config lines indented inside the triple-quoted string
    const ml: string[] = [];
    ml.push(`    net = network();`);
    ml.push(`    net.Nreplicates = ${c.Nreplicates};`);
    ml.push(``);
    ml.push(`    %% ---- Domain ----`);
    ml.push(`    net.domain.b                 = ${c.b};`);
    ml.push(`    net.domain.Lx                = ${c.Lx};`);
    ml.push(`    net.domain.Ly                = ${c.Ly};`);
    ml.push(`    net.domain.scale             = ${c.scale};`);
    ml.push(`    net.domain.boundary          = '${c.boundary}';`);
    if (c.imanualseed) ml.push(`    net.domain.seed              = ${c.seed};`);
    ml.push(`    net.domain.write_location    = '${c.write_location}';`);
    ml.push(`    net.domain.lammps_data_file  = '${c.lammps_data_file}';`);
    ml.push(`    net.domain.lammps_viz_file   = '${c.lammps_viz_file}';`);
    ml.push(`    net.domain.bond_table_file   = '${c.bond_table_file}';`);
    ml.push(``);
    ml.push(`    %% ---- Architecture ----`);
    ml.push(`    net.architecture.geometry           = '${c.geometry}';`);
    ml.push(`    net.architecture.rho_atom           = ${c.rho_atom};`);
    ml.push(`    net.peratom.Max_peratom_bond        = ${c.max_peratom_bond};`);
    ml.push(`    net.peratom.min_degree_keep         = ${c.min_degree_keep};`);
    if (c.geometry === 'hex_lattice') {
      ml.push(`    net.architecture.lattice_spacing            = ${c.lattice_spacing};`);
      ml.push(`    net.architecture.spacing_multiplier_mode    = '${c.spacing_multiplier_mode}';`);
      if (c.spacing_multiplier_mode === 'manual')
        ml.push(`    net.architecture.spacing_multiplier         = ${c.spacing_multiplier};`);
      ml.push(`    net.architecture.lattice_disorder_level     = ${c.lattice_disorder_level};`);
      ml.push(`    net.architecture.lattice_disorder_maxfrac   = ${c.lattice_disorder_maxfrac};`);
      ml.push(`    net.architecture.lattice_max_del_per_node   = ${c.lattice_max_del_per_node};`);
      ml.push(`    net.architecture.lattice_min_degree_keep    = ${c.lattice_min_degree_keep};`);
    }
    ml.push(``);
    ml.push(`    %% ---- Strand typology ----`);
    ml.push(`    net.architecture.strand_typology.mode = '${t}';`);
    if (t === 'mono') {
      ml.push(`    net.architecture.strand_typology.mono.value = ${c.mono_value};`);
    } else if (t === 'polydisperse') {
      pushPolyAssignmentConfig(ml, '    net.architecture.strand_typology', {
        method: c.poly_method,
        minValue: c.poly_min_value,
        rounding: c.poly_rounding,
        alignToLength: c.poly_align,
        rangeMethod: c.poly_range_method,
        targetMin: c.poly_target_min,
        targetMax: c.poly_target_max,
        pmfMean: c.poly_pmf_mean,
        pmfMin: c.poly_pmf_min,
        pmfMax: c.poly_pmf_max,
      });
    } else if (t === 'bimodal') {
      pushBimodalAssignmentConfig(ml, '    net.architecture.strand_typology', {
        method: c.bimodal_method,
        mean1: c.bimodal_mean1,
        mean2: c.bimodal_mean2,
        std1: c.bimodal_std1,
        std2: c.bimodal_std2,
        heightMode: c.bimodal_height_mode,
        heightProb: c.bimodal_height_prob,
        heightCount: c.bimodal_height_count,
        longFirst: c.bimodal_long_first,
        minValue: c.bimodal_min_value,
        doubleNetworkFlag: c.bimodal_double_network_flag,
        alpha: c.bimodal_alpha,
        auto1Flag: c.bimodal_auto_1_flag,
        auto2Flag: c.bimodal_auto_2_flag,
        lam1: c.bimodal_lam_1,
        lam2: c.bimodal_lam_2,
        stdR1: c.bimodal_stdR_1,
        stdR2: c.bimodal_stdR_2,
        binWindowMethod: c.bimodal_bin_window_method,
        manualDevType: c.bimodal_manual_dev_type,
      }, boolStr);
    }
    ml.push(``);
    ml.push(`    %% ---- Perbond ----`);
    ml.push(`    net.perbond.kuhn.auto = ${boolStr(c.kuhn_auto)};`);
    if (!c.kuhn_auto) {
      ml.push(`    net.perbond.kuhn.mode = '${c.kuhn_mode}';`);
      if (c.kuhn_mode === 'mono')
        ml.push(`    net.perbond.kuhn.mono.value = ${c.kuhn_mono_value};`);
      else if (c.kuhn_mode === 'polydisperse') {
        pushPolyAssignmentConfig(ml, '    net.perbond.kuhn', {
          method: c.kuhn_poly_method,
          minValue: c.kuhn_poly_min_value,
          rounding: c.kuhn_poly_rounding,
          alignToLength: c.kuhn_poly_align,
          rangeMethod: c.kuhn_poly_range_method,
          targetMin: c.kuhn_poly_target_min,
          targetMax: c.kuhn_poly_target_max,
          pmfMean: c.kuhn_poly_pmf_mean,
          pmfMin: c.kuhn_poly_pmf_min,
          pmfMax: c.kuhn_poly_pmf_max,
        });
      } else if (c.kuhn_mode === 'bimodal') {
        pushBimodalAssignmentConfig(ml, '    net.perbond.kuhn', {
          method: c.kuhn_bimodal_method,
          mean1: c.kuhn_bimodal_mean1,
          mean2: c.kuhn_bimodal_mean2,
          std1: c.kuhn_bimodal_std1,
          std2: c.kuhn_bimodal_std2,
          heightMode: c.kuhn_bimodal_height_mode,
          heightProb: c.kuhn_bimodal_height_prob,
          heightCount: c.kuhn_bimodal_height_count,
          longFirst: c.kuhn_bimodal_long_first,
          minValue: c.kuhn_bimodal_min_value,
          doubleNetworkFlag: c.kuhn_bimodal_double_network_flag,
          alpha: c.kuhn_bimodal_alpha,
          auto1Flag: c.kuhn_bimodal_auto_1_flag,
          auto2Flag: c.kuhn_bimodal_auto_2_flag,
          lam1: c.kuhn_bimodal_lam_1,
          lam2: c.kuhn_bimodal_lam_2,
          stdR1: c.kuhn_bimodal_stdR_1,
          stdR2: c.kuhn_bimodal_stdR_2,
          binWindowMethod: c.kuhn_bimodal_bin_window_method,
          manualDevType: c.kuhn_bimodal_manual_dev_type,
        }, boolStr);
      }
    }
    if (c.idefect) {
      ml.push(``);
      ml.push(`    %% ---- Defects ----`);
      ml.push(`    net.defect.density_mode       = '${c.defect_density_mode}';`);
      if (c.defect_density_mode === 'count')
        ml.push(`    net.defect.n_voids            = ${c.defect_n_voids};`);
      else
        ml.push(`    net.defect.void_area_frac     = ${c.defect_void_area_frac};`);
      ml.push(`    net.defect.size_dist          = '${c.defect_size_dist}';`);
      ml.push(`    net.defect.radius_mean        = ${c.defect_radius_mean};`);
      if (c.defect_size_dist !== 'fixed') {
        ml.push(`    net.defect.radius_std         = ${c.defect_radius_std};`);
        ml.push(`    net.defect.radius_min         = ${c.defect_radius_min};`);
        ml.push(`    net.defect.radius_max         = ${c.defect_radius_max};`);
      }
      ml.push(`    net.defect.shape_roughness    = ${c.defect_shape_roughness};`);
      ml.push(`    net.defect.shape_n_modes      = ${c.defect_shape_n_modes};`);
      ml.push(`    net.defect.void_overlap       = ${boolStr(c.defect_void_overlap)};`);
      ml.push(`    net.defect.center_distribution = '${c.defect_center_dist}';`);
      if (c.defect_center_dist === 'clustered') {
        ml.push(`    net.defect.n_cluster_parents  = ${c.defect_n_cluster_parents};`);
        ml.push(`    net.defect.cluster_spread     = ${c.defect_cluster_spread};`);
      }
      ml.push(`    net.defect.margin_frac        = ${c.defect_margin_frac};`);
      ml.push(`    net.defect.prune_isolated     = ${boolStr(c.defect_prune_isolated)};`);
      ml.push(`    net.defect.sparse_network     = ${boolStr(c.defect_sparse_network)};`);
      ml.push(`    net.defect.wall_thickness     = ${c.defect_wall_thickness};`);
      ml.push(`    net.defect.clamp_thickness    = ${c.defect_clamp_thickness};`);
      ml.push(`    net.defect.bridge_width       = ${c.defect_bridge_width};`);
      if (c.defect_thinning) {
        ml.push(`    net.defect.thinning           = true;`);
        ml.push(`    net.defect.thinning_radius    = ${c.defect_thinning_radius};`);
        ml.push(`    net.defect.thinning_target_frac = ${c.defect_thinning_target_frac};`);
        ml.push(`    net.defect.thinning_min_keep  = ${c.defect_thinning_min_keep};`);
      }
      if (c.defect_bridging) {
        ml.push(`    net.defect.bridging           = true;`);
        ml.push(`    net.defect.bridge_max_dist    = ${c.defect_bridge_max_dist};`);
        ml.push(`    net.defect.bridge_void_thresh = ${c.defect_bridge_void_thresh};`);
        ml.push(`    net.defect.bridge_perp_width  = ${c.defect_bridge_perp_width};`);
        ml.push(`    net.defect.bridge_max_degree  = ${c.defect_bridge_max_degree};`);
        ml.push(`    net.defect.bridge_max_bonds   = ${c.defect_bridge_max_bonds};`);
        ml.push(`    net.defect.bridge_min_spacing = ${c.defect_bridge_min_spacing};`);
      }
    }
    if (c.use_multitype) {
      ml.push(``);
      ml.push(`    %% ---- Multi-type ----`);
      ml.push(`    net.architecture.types.enabled     = true;`);
      ml.push(`    net.architecture.types.natom_type  = ${c.natom_type};`);
      ml.push(`    net.architecture.types.nbond_type  = ${c.nbond_type};`);
      ml.push(`    net.architecture.types.atype_mode  = '${c.atype_mode}';`);
      ml.push(`    net.architecture.types.btype_mode  = '${c.btype_mode}';`);
      if (c.atype_mode === 'fixed')
        ml.push(`    net.architecture.types.atom_count  = ${atomTargets};`);
      else
        ml.push(`    net.architecture.types.atom_frac   = ${atomTargets};`);
      if (c.btype_mode === 'fixed')
        ml.push(`    net.architecture.types.bond_count  = ${bondTargets};`);
      else
        ml.push(`    net.architecture.types.bond_frac   = ${bondTargets};`);
      ml.push(...connectivityLines);
    }
    if (c.ipotential) {
      ml.push(``);
      ml.push(`    %% ---- Potential ----`);
      ml.push(`    net.pot.k_LD    = ${c.pot_k_LD};`);
      ml.push(`    net.pot.N_rho   = ${c.pot_N_rho};`);
      ml.push(`    net.pot.rho_min = ${c.pot_rho_min};`);
      ml.push(`    net.pot.rho_max = ${c.pot_rho_max};`);
    }
    ml.push(``);
    ml.push(`    %% ---- Flags ----`);
    ml.push(`    net.flags.isave       = ${boolStr(c.isave)};`);
    ml.push(`    net.flags.iplot       = false;  % plotting disabled in Python/Octave`);
    ml.push(`    net.flags.ilog        = ${boolStr(c.ilog)};`);
    ml.push(`    net.flags.savemode    = ${boolStr(c.savemode)};`);
    ml.push(`    net.flags.imanualseed = ${boolStr(c.imanualseed)};`);
    ml.push(`    net.flags.idefect     = ${boolStr(c.idefect)};`);
    ml.push(`    net.flags.ipotential  = ${boolStr(c.ipotential)};`);
    ml.push(`    net.flags.idumpsettings = ${boolStr(c.idumpsettings)};`);
    ml.push(`    net.flags.iversbose_settings = ${boolStr(c.iversbose_settings)};`);
    ml.push(``);
    ml.push(`    %% ---- Generate ----`);
    ml.push(`    net.generateNetwork();`);

    lines.push(ml.join('\n'));
    lines.push(`""")`);

    return lines.join('\n');
  }

  function generateScript() {
    return lang === 'matlab' ? generateMatlab() : generatePython();
  }

  function handleCopy() {
    navigator.clipboard.writeText(generateScript());
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  }

  function handleDownload() {
    const ext = lang === 'matlab' ? '.m' : '.py';
    const blob = new Blob([generateScript()], { type: 'text/plain' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `networkgen_config${ext}`;
    a.click();
  }

  const t = cfg.typology_mode;
  const atomFractionValues = sanitizeWeightArray(cfg.atom_frac_values, cfg.natom_type);
  const bondFractionValues = sanitizeWeightArray(cfg.bond_frac_values, cfg.nbond_type);
  const atomCountValues = sanitizeCountArray(cfg.atom_count_values, cfg.natom_type);
  const bondCountValues = sanitizeCountArray(cfg.bond_count_values, cfg.nbond_type);
  const connectivityRules = sanitizeConnectivityRules(cfg.connectivity_rules, cfg.natom_type, cfg.nbond_type);
  const atomFractionSum = atomFractionValues.reduce((sum, value) => sum + value, 0);
  const bondFractionSum = bondFractionValues.reduce((sum, value) => sum + value, 0);
  const atomTypeOptions = Array.from({ length: cfg.natom_type }, (_, index) => index + 1);
  const bondTypeOptions = Array.from({ length: cfg.nbond_type }, (_, index) => index + 1);
  const topologyPolyFields = {
    method: 'poly_method',
    minValue: 'poly_min_value',
    rounding: 'poly_rounding',
    align: 'poly_align',
    rangeMethod: 'poly_range_method',
    targetMin: 'poly_target_min',
    targetMax: 'poly_target_max',
    pmfMean: 'poly_pmf_mean',
    pmfMin: 'poly_pmf_min',
    pmfMax: 'poly_pmf_max',
  };
  const kuhnPolyFields = {
    method: 'kuhn_poly_method',
    minValue: 'kuhn_poly_min_value',
    rounding: 'kuhn_poly_rounding',
    align: 'kuhn_poly_align',
    rangeMethod: 'kuhn_poly_range_method',
    targetMin: 'kuhn_poly_target_min',
    targetMax: 'kuhn_poly_target_max',
    pmfMean: 'kuhn_poly_pmf_mean',
    pmfMin: 'kuhn_poly_pmf_min',
    pmfMax: 'kuhn_poly_pmf_max',
  };
  const topologyBimodalFields = {
    method: 'bimodal_method',
    mean1: 'bimodal_mean1',
    mean2: 'bimodal_mean2',
    std1: 'bimodal_std1',
    std2: 'bimodal_std2',
    heightMode: 'bimodal_height_mode',
    heightProb: 'bimodal_height_prob',
    heightCount: 'bimodal_height_count',
    longFirst: 'bimodal_long_first',
    minValue: 'bimodal_min_value',
    doubleNetworkFlag: 'bimodal_double_network_flag',
    alpha: 'bimodal_alpha',
    auto1: 'bimodal_auto_1_flag',
    auto2: 'bimodal_auto_2_flag',
    lam1: 'bimodal_lam_1',
    lam2: 'bimodal_lam_2',
    stdR1: 'bimodal_stdR_1',
    stdR2: 'bimodal_stdR_2',
    binWindowMethod: 'bimodal_bin_window_method',
    manualDevType: 'bimodal_manual_dev_type',
  };
  const kuhnBimodalFields = {
    method: 'kuhn_bimodal_method',
    mean1: 'kuhn_bimodal_mean1',
    mean2: 'kuhn_bimodal_mean2',
    std1: 'kuhn_bimodal_std1',
    std2: 'kuhn_bimodal_std2',
    heightMode: 'kuhn_bimodal_height_mode',
    heightProb: 'kuhn_bimodal_height_prob',
    heightCount: 'kuhn_bimodal_height_count',
    longFirst: 'kuhn_bimodal_long_first',
    minValue: 'kuhn_bimodal_min_value',
    doubleNetworkFlag: 'kuhn_bimodal_double_network_flag',
    alpha: 'kuhn_bimodal_alpha',
    auto1: 'kuhn_bimodal_auto_1_flag',
    auto2: 'kuhn_bimodal_auto_2_flag',
    lam1: 'kuhn_bimodal_lam_1',
    lam2: 'kuhn_bimodal_lam_2',
    stdR1: 'kuhn_bimodal_stdR_1',
    stdR2: 'kuhn_bimodal_stdR_2',
    binWindowMethod: 'kuhn_bimodal_bin_window_method',
    manualDevType: 'kuhn_bimodal_manual_dev_type',
  };

  function renderTargetEditor(title, key, values, mode, summary, resetValues, count) {
    const isFraction = mode === 'frac';

    return (
      <div className={styles.targetCard}>
        <div className={styles.targetCardHead}>
          <div>
            <div className={styles.targetCardTitle}>{title}</div>
            <div className={styles.targetCardSubtitle}>
              {isFraction ? 'Fractions / weights per exported type' : 'Exact target count per exported type'}
            </div>
          </div>
          <div className={styles.targetCardActions}>
            {isFraction && (
              <button
                type="button"
                className={styles.actionButton}
                onClick={() => normalizeArrayValues(key, count)}
              >
                Normalize
              </button>
            )}
            <button
              type="button"
              className={styles.actionButton}
              onClick={() => resetArrayValues(key, resetValues)}
            >
              {isFraction ? 'Equalize' : 'Fill 1'}
            </button>
          </div>
        </div>

        <div className={styles.targetGrid}>
          {values.map((value, index) => (
            <label key={`${key}-${index}`} className={styles.targetCell}>
              <span className={styles.targetCellLabel}>Type {index + 1}</span>
              <input
                type="number"
                min={0}
                step={isFraction ? 0.01 : 1}
                value={value}
                onChange={e => setArrayValue(key, index, parseFloat(e.target.value) || 0, !isFraction)}
              />
            </label>
          ))}
        </div>

        <div className={styles.targetSummary}>
          <span className={isFraction && Math.abs(summary - 1) > 0.001 ? styles.targetSummaryWarn : undefined}>
            {isFraction ? `Current sum: ${summary.toFixed(3)}` : `Current total: ${summary.toFixed(0)}`}
          </span>
          {isFraction && <span>Any positive weights are normalized internally by NetworkGen.</span>}
        </div>
      </div>
    );
  }

  function renderPolyEditor(fields) {
    const method = cfg[fields.method];

    return (
      <>
        <Row label="Method">{sel(fields.method, [['pmf'], ['range'], ['geom']])}</Row>
        <Row label="Minimum value">{num(fields.minValue, 1)}</Row>
        {method === 'pmf' && <>
          <Row label="PMF mean">{num(fields.pmfMean, 1)}</Row>
          <Row label="PMF min">{num(fields.pmfMin, 1)}</Row>
          <Row label="PMF max">{num(fields.pmfMax, 1)}</Row>
        </>}
        {method === 'range' && <>
          <Row label="Range method">{sel(fields.rangeMethod, [['rank'], ['linear']])}</Row>
          <Row label="Target min">{num(fields.targetMin, 1)}</Row>
          <Row label="Target max">{num(fields.targetMax, 1)}</Row>
        </>}
        <Row label="Rounding">{sel(fields.rounding, [['round'], ['ceil'], ['floor']])}</Row>
        <Row label="Align to length">{sel(fields.align, [['none'], ['ascend']])}</Row>
      </>
    );
  }

  function renderBimodalEditor(fields) {
    const method = cfg[fields.method];
    const heightMode = cfg[fields.heightMode];
    const binWindowMethod = cfg[fields.binWindowMethod];

    return (
      <>
        <Row label="Method">{sel(fields.method, [['gaussian'], ['geom'], ['single', 'single (fixed mean)']])}</Row>
        <Row label="Mean 1">{num(fields.mean1, 1)}</Row>
        <Row label="Mean 2">{num(fields.mean2, 1)}</Row>
        {method !== 'single' && <>
          <Row label="Std 1">{num(fields.std1, 0, null, 0.5)}</Row>
          <Row label="Std 2">{num(fields.std2, 0, null, 0.5)}</Row>
        </>}
        <Row label="Minimum value">{num(fields.minValue, 1)}</Row>
        <Row label="Height mode">{sel(fields.heightMode, [['prob'], ['count']])}</Row>
        {heightMode === 'prob'
          ? <Row label="Fraction in mode 2">{slide(fields.heightProb, 0.05, 0.95, 0.05)}</Row>
          : <Row label="Count in mode 2">{num(fields.heightCount, 0, null, 1)}</Row>
        }
        <Row label="Long first">{chk(fields.longFirst)}</Row>

        <Sub title="Advanced bimodal">
          <Row label="Double network">{chk(fields.doubleNetworkFlag)}</Row>
          {cfg[fields.doubleNetworkFlag] && <Row label="Alpha">{num(fields.alpha, 0.01, null, 0.1)}</Row>}
          <Row label="Auto mode 1">{chk(fields.auto1)}</Row>
          {cfg[fields.auto1] && <Row label="lam_1">{num(fields.lam1, 0, 1, 0.01)}</Row>}
          <Row label="Auto mode 2">{chk(fields.auto2)}</Row>
          {cfg[fields.auto2] && <Row label="lam_2">{num(fields.lam2, 0, 1, 0.01)}</Row>}
          <Row label="stdR_1">{num(fields.stdR1, 0, null, 0.1)}</Row>
          <Row label="stdR_2">{num(fields.stdR2, 0, null, 0.1)}</Row>
          <Row label="Bin window">{sel(fields.binWindowMethod, [['manual'], ['adaptive']])}</Row>
          {binWindowMethod === 'manual' && (
            <Row label="Manual deviation">{sel(fields.manualDevType, [['mixed'], ['kuhn'], ['both']])}</Row>
          )}
        </Sub>
      </>
    );
  }

  return (
    <div className={styles.outer}>
      <div className={styles.formCol}>

        {/* Language toggle — sits above the first section */}
        <div style={{ display: 'flex', alignItems: 'center', gap: 12, paddingBottom: 4 }}>
          <span style={{ fontSize: 12, color: 'var(--ifm-color-content-secondary)', fontWeight: 500 }}>
            Output language
          </span>
          <LangToggle lang={lang} setLang={setLang} />
        </div>

        {/* Python-only path settings */}
        {lang === 'python' && (
          <Section id="flags" title="Python / Octave settings" defaultOpen>
            <div className={styles.note}>
              Set the path to your Octave executable and your NetworkGen folder.
              These are written into the generated Python script.
            </div>
            <Row label="NetworkGen path" hint="folder containing .m files">
              {txt('py_networkgen_path')}
            </Row>
          </Section>
        )}

        <Section id="domain" title="Domain" defaultOpen>
          <Row label="b" hint="lengthscale">{num('b', 0.1, null, 0.1)}</Row>
          <Row label="Lx" hint="units of b">{num('Lx', 1)}</Row>
          <Row label="Ly" hint="units of b">{num('Ly', 1)}</Row>
          <Row label="Scale">{num('scale', 0.1, null, 0.1)}</Row>
          <Row label="Boundary">{sel('boundary', [['fixed'], ['periodic']])}</Row>
          <Row label="Networks to generate" hint="sets net.Nreplicates">{num('Nreplicates', 1, null, 1)}</Row>
          <div className={styles.note}>
            NetworkGen currently generates 2D networks, so the domain builder only exposes the in-plane size.
          </div>
          <Row label="Manual seed">
            {chk('imanualseed')}
          </Row>
          {cfg.imanualseed && <Row label="Seed">{num('seed', 1)}</Row>}
          <Row label="Output folder">{txt('write_location')}</Row>
          <Row label="Data file prefix">{txt('lammps_data_file')}</Row>
          <Row label="Viz file prefix">{txt('lammps_viz_file')}</Row>
          <Row label="Bond table prefix">{txt('bond_table_file')}</Row>
          <div className={styles.note}>
            Batch outputs already get unique replicate suffixes automatically. Leave sample numbering at the package default unless you are managing your own outer loop by hand.
          </div>
        </Section>

        <Section id="architecture" title="Architecture">
          <Row label="Geometry">{sel('geometry', [['random'], ['hex_lattice', 'hex lattice']])}</Row>
          <Row label="rho_atom" hint="atoms/unit area">{num('rho_atom', 0.0001, null, 0.0001)}</Row>
          <Row label="Max bonds/atom">{num('max_peratom_bond', 3, null, 1)}</Row>
          <Row label="Min degree keep">{num('min_degree_keep', 1, null, 1)}</Row>
          {cfg.geometry === 'hex_lattice' && (
            <Sub title="Lattice settings">
              <Row label="Lattice spacing">{num('lattice_spacing', 1, null, 0.5)}</Row>
              <Row label="Spacing mode">{sel('spacing_multiplier_mode', [['auto'], ['manual']])}</Row>
              {cfg.spacing_multiplier_mode === 'manual' && (
                <Row label="Spacing multiplier">{num('spacing_multiplier', 0, null, 0.1)}</Row>
              )}
              <Row label="Disorder level">{slide('lattice_disorder_level', 0, 1, 0.05)}</Row>
              <Row label="Disorder maxfrac">{slide('lattice_disorder_maxfrac', 0, 1, 0.05)}</Row>
              <Row label="Max del/node">{num('lattice_max_del_per_node', 0, null, 1)}</Row>
              <Row label="Min degree keep">{num('lattice_min_degree_keep', 3, null, 1)}</Row>
            </Sub>
          )}
        </Section>

        <Section id="typology" title="Strand typology">
          <Row label="Mode">{sel('typology_mode', [['mono'], ['polydisperse'], ['bimodal']])}</Row>
          {t === 'mono' && (
            <Sub title="Mono">
              <Row label="Kuhn value">{num('mono_value', 1)}</Row>
            </Sub>
          )}
          {t === 'polydisperse' && (
            <Sub title="Polydisperse">
              {renderPolyEditor(topologyPolyFields)}
            </Sub>
          )}
          {t === 'bimodal' && (
            <Sub title="Bimodal">
              {renderBimodalEditor(topologyBimodalFields)}
            </Sub>
          )}
        </Section>

        <Section id="perbond" title="Perbond">
          <Row label="Kuhn auto" hint="copies typology">{chk('kuhn_auto')}</Row>
          {!cfg.kuhn_auto && (
            <Sub title="Kuhn manual distribution">
              <Row label="Mode">{sel('kuhn_mode', [['mono'], ['polydisperse'], ['bimodal']])}</Row>
              {cfg.kuhn_mode === 'mono' && (
                <Row label="Kuhn value">{num('kuhn_mono_value', 1)}</Row>
              )}
              {cfg.kuhn_mode === 'polydisperse' && renderPolyEditor(kuhnPolyFields)}
              {cfg.kuhn_mode === 'bimodal' && renderBimodalEditor(kuhnBimodalFields)}
            </Sub>
          )}
        </Section>

        <Section id="defect" title="Defects">
          <Row label="Enable defects">{chk('idefect')}</Row>
          {cfg.idefect && (<>
            <Row label="Density mode">{sel('defect_density_mode', [['count'], ['area_frac', 'area fraction']])}</Row>
            {cfg.defect_density_mode === 'count'
              ? <Row label="N voids">{num('defect_n_voids', 0)}</Row>
              : <Row label="Void area frac">{num('defect_void_area_frac', 0, 1, 0.01)}</Row>
            }
            <Row label="Size distribution">{sel('defect_size_dist', [['gaussian'], ['fixed'], ['exponential']])}</Row>
            <Row label="Radius mean">{num('defect_radius_mean', 0, null, 0.5)}</Row>
            {cfg.defect_size_dist !== 'fixed' && <>
              <Row label="Radius std">{num('defect_radius_std', 0, null, 0.5)}</Row>
              <Row label="Radius min">{num('defect_radius_min', 0, null, 0.5)}</Row>
              <Row label="Radius max">{num('defect_radius_max', 0, null, 0.5)}</Row>
            </>}
            <Sub title="Shape">
              <Row label="Roughness">{slide('defect_shape_roughness', 0, 1, 0.05)}</Row>
              <Row label="N modes">{num('defect_shape_n_modes', 1)}</Row>
              <Row label="Void overlap">{chk('defect_void_overlap')}</Row>
            </Sub>
            <Sub title="Placement">
              <Row label="Center distribution">{sel('defect_center_dist', [['random'], ['uniform'], ['clustered']])}</Row>
              {cfg.defect_center_dist === 'clustered' && <>
                <Row label="N cluster parents">{num('defect_n_cluster_parents', 1)}</Row>
                <Row label="Cluster spread">{num('defect_cluster_spread', 0, null, 0.5)}</Row>
              </>}
              <Row label="Margin frac">{num('defect_margin_frac', 0, 0.5, 0.01)}</Row>
              <Row label="Bridge width">{num('defect_bridge_width', 0, null, 0.5)}</Row>
            </Sub>
            <Sub title="Cleanup">
              <Row label="Prune isolated">{chk('defect_prune_isolated')}</Row>
              <Row label="Sparse network">{chk('defect_sparse_network')}</Row>
              <Row label="Wall thickness">{num('defect_wall_thickness', 0, null, 0.5)}</Row>
              <Row label="Clamp thickness">{num('defect_clamp_thickness', 0, null, 0.01)}</Row>
            </Sub>
            <Sub title="Advanced passes">
              <Row label="Density thinning">{chk('defect_thinning')}</Row>
              {cfg.defect_thinning && <>
                <Row label="Thinning radius">{num('defect_thinning_radius', 0, null, 0.5)}</Row>
                <Row label="Target keep frac">{num('defect_thinning_target_frac', 0, 1, 0.01)}</Row>
                <Row label="Min keep frac">{num('defect_thinning_min_keep', 0, 1, 0.01)}</Row>
              </>}
              <Row label="Constriction bridging">{chk('defect_bridging')}</Row>
              {cfg.defect_bridging && <>
                <Row label="Bridge max dist">{num('defect_bridge_max_dist', 0, null, 0.5)}</Row>
                <Row label="Void threshold">{num('defect_bridge_void_thresh', 0, 1, 0.01)}</Row>
                <Row label="Perp width">{num('defect_bridge_perp_width', 0, null, 0.5)}</Row>
                <Row label="Max degree">{num('defect_bridge_max_degree', 0, null, 1)}</Row>
                <Row label="Max bonds">{num('defect_bridge_max_bonds', 0, null, 1)}</Row>
                <Row label="Min spacing">{num('defect_bridge_min_spacing', 0, null, 0.5)}</Row>
              </>}
            </Sub>
          </>)}
        </Section>

        <Section id="multitype" title="Multi-type">
          <Row label="Enable multi-type">{chk('use_multitype')}</Row>
          {cfg.use_multitype && (<>
            <Row label="N atom types">
              <input
                type="number"
                value={cfg.natom_type}
                min={1}
                step={1}
                onChange={e => setAtomTypeCount(parseFloat(e.target.value) || 1)}
              />
            </Row>
            <Row label="N bond types">
              <input
                type="number"
                value={cfg.nbond_type}
                min={1}
                step={1}
                onChange={e => setBondTypeCount(parseFloat(e.target.value) || 1)}
              />
            </Row>
            <Row label="Atom type mode">{sel('atype_mode', [['frac', 'fraction'], ['fixed', 'fixed count']])}</Row>
            <Row label="Bond type mode">{sel('btype_mode', [['frac', 'fraction'], ['fixed', 'fixed count']])}</Row>

            <Sub title="Type targets">
              <div className={styles.targetEditorGrid}>
                {renderTargetEditor(
                  cfg.atype_mode === 'frac' ? 'Atom fractions' : 'Atom counts',
                  cfg.atype_mode === 'frac' ? 'atom_frac_values' : 'atom_count_values',
                  cfg.atype_mode === 'frac' ? atomFractionValues : atomCountValues,
                  cfg.atype_mode,
                  cfg.atype_mode === 'frac'
                    ? atomFractionSum
                    : atomCountValues.reduce((sum, value) => sum + value, 0),
                  cfg.atype_mode === 'frac' ? uniformFractions(cfg.natom_type) : unitCounts(cfg.natom_type),
                  cfg.natom_type
                )}
                {renderTargetEditor(
                  cfg.btype_mode === 'frac' ? 'Bond fractions' : 'Bond counts',
                  cfg.btype_mode === 'frac' ? 'bond_frac_values' : 'bond_count_values',
                  cfg.btype_mode === 'frac' ? bondFractionValues : bondCountValues,
                  cfg.btype_mode,
                  cfg.btype_mode === 'frac'
                    ? bondFractionSum
                    : bondCountValues.reduce((sum, value) => sum + value, 0),
                  cfg.btype_mode === 'frac' ? uniformFractions(cfg.nbond_type) : unitCounts(cfg.nbond_type),
                  cfg.nbond_type
                )}
              </div>
            </Sub>

            <Sub title="Connectivity rules">
              <div className={styles.note}>
                Add rows of the form `[atomTypeA atomTypeB bondType allowed]`. Leave the table empty to allow all combinations. Atom-type order is symmetric, so `(1,2)` is treated the same as `(2,1)`.
              </div>

              <div className={styles.ruleToolbar}>
                <button type="button" className={styles.actionButton} onClick={addConnectivityRule}>
                  Add rule
                </button>
                {connectivityRules.length > 0 && (
                  <button
                    type="button"
                    className={`${styles.actionButton} ${styles.actionButtonGhost}`}
                    onClick={() => resetArrayValues('connectivity_rules', [])}
                  >
                    Clear rules
                  </button>
                )}
              </div>

              {connectivityRules.length === 0 ? (
                <div className={styles.ruleEmptyState}>
                  No explicit rules yet. The generated script will leave `types.connectivity = []`, which means all type combinations are allowed.
                </div>
              ) : (
                <div className={styles.ruleTableWrap}>
                  <table className={styles.ruleTable}>
                    <thead>
                      <tr>
                        <th>Atom A</th>
                        <th>Atom B</th>
                        <th>Bond type</th>
                        <th>Action</th>
                        <th />
                      </tr>
                    </thead>
                    <tbody>
                      {connectivityRules.map((rule, index) => (
                        <tr key={`rule-${index}`}>
                          <td>
                            <select
                              value={rule.atomTypeA}
                              onChange={e => updateConnectivityRule(index, 'atomTypeA', parseFloat(e.target.value) || 1)}
                            >
                              {atomTypeOptions.map(typeId => (
                                <option key={`rule-a-${typeId}`} value={typeId}>Type {typeId}</option>
                              ))}
                            </select>
                          </td>
                          <td>
                            <select
                              value={rule.atomTypeB}
                              onChange={e => updateConnectivityRule(index, 'atomTypeB', parseFloat(e.target.value) || 1)}
                            >
                              {atomTypeOptions.map(typeId => (
                                <option key={`rule-b-${typeId}`} value={typeId}>Type {typeId}</option>
                              ))}
                            </select>
                          </td>
                          <td>
                            <select
                              value={rule.bondType}
                              onChange={e => updateConnectivityRule(index, 'bondType', parseFloat(e.target.value) || 1)}
                            >
                              {bondTypeOptions.map(typeId => (
                                <option key={`rule-bond-${typeId}`} value={typeId}>Type {typeId}</option>
                              ))}
                            </select>
                          </td>
                          <td>
                            <select
                              value={rule.allowed ? 'allow' : 'forbid'}
                              onChange={e => updateConnectivityRule(index, 'allowed', e.target.value === 'allow')}
                            >
                              <option value="forbid">Forbid</option>
                              <option value="allow">Allow</option>
                            </select>
                          </td>
                          <td className={styles.ruleDeleteCell}>
                            <button
                              type="button"
                              className={`${styles.actionButton} ${styles.actionButtonGhost}`}
                              onClick={() => removeConnectivityRule(index)}
                            >
                              Remove
                            </button>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </Sub>
          </>)}
        </Section>

        <Section id="potential" title="Potential">
          <Row label="Enable potential">{chk('ipotential')}</Row>
          {cfg.ipotential && (<>
            <div className={styles.note}>Generates lookup table for LAMMPS <code>pair_style local/density</code> (cohesion). Equilibrium separation for <code>pair_style bpm/spring</code> is written to the log file.</div>
            <Row label="k_LD">{num('pot_k_LD', 0, null, 0.001)}</Row>
            <Row label="N_rho" hint="table resolution">{num('pot_N_rho', 100, null, 1000)}</Row>
            <Row label="rho_min">{num('pot_rho_min', 0, null, 0.01)}</Row>
            <Row label="rho_max">{num('pot_rho_max', 0, null, 1)}</Row>
          </>)}
        </Section>

        <Section id="flags" title="Flags & output">
          <Row label="Save files">{chk('isave')}</Row>
          {lang === 'matlab' && <Row label="Plot">{chk('iplot')}</Row>}
          {lang === 'python' && (
            <Row label="Plot" hint="disabled in Python/Octave">
              <span style={{ fontSize: 12, color: 'var(--ifm-color-secondary-darkest)' }}>
                always off (no display)
              </span>
            </Row>
          )}
          <Row label="Write log">{chk('ilog')}</Row>
          <Row label="Auto-name files">{chk('savemode')}</Row>
          <Row label="Dump settings">{chk('idumpsettings')}</Row>
          <Row label="Verbose settings dump">{chk('iversbose_settings')}</Row>
        </Section>

      </div>

      <div className={styles.codeCol}>
        <div className={styles.codeWrap}>
          <div className={styles.codeHead}>
            <span className={styles.codeTitle}>
              {lang === 'matlab' ? 'MATLAB script (.m)' : 'Python script (.py)'}
            </span>
            <div className={styles.codeActions}>
              <button onClick={handleCopy}>{copied ? 'Copied!' : 'Copy'}</button>
              <button onClick={handleDownload}>
                Download {lang === 'matlab' ? '.m' : '.py'}
              </button>
            </div>
          </div>
          <pre className={styles.code}>{generateScript()}</pre>
        </div>
      </div>

    </div>
  );
}
