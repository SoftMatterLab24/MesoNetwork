import React, { useState } from 'react';
import styles from './styles.module.css';

const CLASSES = [
  {
    id: 'network',
    name: 'network',
    kind: 'handle class',
    color: '#378ADD',
    description: 'Root class. Instantiate this to configure and generate a network.',
    fields: [
      { name: 'Nreplicates', type: 'int', default: '1', desc: 'Number of network replicates to generate in a single call.' },
      { name: 'flags', type: 'struct', default: '—', desc: 'Boolean control flags: isave, iplot, ilog, savemode, imanualseed, idefect, ipotential.' },
      { name: 'domain', type: 'struct', default: '—', desc: 'Domain settings: b, Lx, Ly, scale, boundary, seed, write_location, file prefixes.' },
      { name: 'peratom', type: 'struct', default: '—', desc: 'Per-atom constraints: Max_peratom_bond, min_degree_keep.' },
      { name: 'perbond', type: 'bondstyle', default: 'bondstyle()', desc: 'Per-bond property assignments. Subclass — see bondstyle.', isRef: true, refId: 'bondstyle' },
      { name: 'defect', type: 'struct', default: '—', desc: 'Defect settings: density_mode, n_voids, size_dist, radius_mean, center_distribution, etc.' },
      { name: 'pot', type: 'struct', default: '—', desc: 'Potential table settings: k_LD, N_rho, rho_min, rho_max.' },
      { name: 'architecture', type: 'architecture', default: 'architecture()', desc: 'Network geometry and typology. Subclass — see architecture.', isRef: true, refId: 'architecture' },
      { name: 'log', type: 'networklog', default: 'networklog()', desc: 'Replicate-aware logging helper used for printed output and settings dumps.' },
    ],
    methods: [
      { name: 'generateNetwork()', desc: 'Main generation loop. Runs all build steps for each replicate.' },
    ],
  },
  {
    id: 'architecture',
    name: 'architecture',
    kind: 'class',
    color: '#1D9E75',
    description: 'Controls node placement geometry, lattice settings, and strand typology.',
    fields: [
      { name: 'geometry', type: 'string', default: "'random'", desc: "Node placement pattern. Options: 'random', 'hex_lattice'." },
      { name: 'strand_typology', type: 'assignmentmode', default: 'assignmentmode()', desc: 'Target bond length distribution. Uses assignmentmode framework.', isRef: true, refId: 'assignmentmode' },
      { name: 'types', type: 'struct', default: '—', desc: 'Multi-type settings: enabled flag, atom/bond target arrays, and connectivity rows [atomTypeA atomTypeB bondType allowed].' },
      { name: 'lattice_spacing', type: 'double', default: '6', desc: 'Nominal node spacing in units of b. Used when geometry = hex_lattice.' },
      { name: 'spacing_multiplier_mode', type: 'string', default: "'auto'", desc: "How spacing multiplier is set. Options: 'auto', 'manual'." },
      { name: 'spacing_multiplier', type: 'double', default: '1.2', desc: 'Manual scale factor applied to lattice_spacing.' },
      { name: 'lattice_disorder_level', type: 'double', default: '1', desc: 'Magnitude of random positional perturbations. Range [0, 1].' },
      { name: 'lattice_disorder_maxfrac', type: 'double', default: '0.4', desc: 'Max fractional displacement from ideal lattice position.' },
      { name: 'lattice_max_del_per_node', type: 'int', default: '1', desc: 'Max bonds deletable per node during topological disorder.' },
      { name: 'lattice_min_degree_keep', type: 'int', default: '5', desc: 'Minimum node degree to preserve during topological disorder.' },
      { name: 'rho_atom', type: 'double', default: '0.0078', desc: 'Target atom number density (atoms per unit area in 2D).' },
    ],
  },
  {
    id: 'bondstyle',
    name: 'bondstyle',
    kind: 'class',
    color: '#1D9E75',
    description: 'Defines per-bond physical properties and their distributions.',
    fields: [
      { name: 'type', type: 'int', default: '1', desc: 'Bond style type index for LAMMPS output.' },
      { name: 'kuhn', type: 'assignmentmode', default: 'assignmentmode()', desc: 'Kuhn segment count per bond. When auto=true, copies strand_typology distribution exactly.', isRef: true, refId: 'assignmentmode' },
    ],
  },
  {
    id: 'assignmentmode',
    name: 'assignmentmode',
    kind: 'class (shared)',
    color: '#7F77DD',
    description: 'General statistical distribution framework. Used by strand_typology and all perbond properties.',
    fields: [
      { name: 'auto', type: 'boolean', default: 'true', desc: 'When true, copies the parent property distribution (e.g. kuhn copies strand_typology).' },
      { name: 'mode', type: 'string', default: "'mono'", desc: "Active distribution mode. Options: 'mono', 'uniform', 'poly', 'polydisperse' (alias), 'bimodal'." },
      { name: 'mono', type: 'struct', default: '{value: 20}', desc: 'Single fixed value assignment. Field: value.' },
      { name: 'uniform', type: 'struct', default: '{min: 5, max: 40}', desc: 'Flat distribution. Fields: min_value, max_value.' },
      { name: 'poly', type: 'struct', default: '{method: pmf, ...}', desc: "Continuous distribution. Fields: method ('geom','range','pmf'), pmf_mean, pmf_min, pmf_max, rounding, align_to_length, etc." },
      { name: 'bimodal', type: 'struct', default: '{mean_1: 35, ...}', desc: "Two-population mixture. Core: mean_1/2, std_1/2, method ('single','geom','gaussian'), height_mode, height_prob, height_count. Advanced: lam_1/2, stdR_1/2, double_network_flag, alpha, auto_1/2_flag, bin_window_method." },
    ],
  },
];

const CONNECTIONS = [
  { from: 'network', field: 'architecture', to: 'architecture', color: '#1D9E75' },
  { from: 'network', field: 'perbond', to: 'bondstyle', color: '#1D9E75' },
  { from: 'architecture', field: 'strand_typology', to: 'assignmentmode', color: '#7F77DD', dashed: true },
  { from: 'bondstyle', field: 'kuhn', to: 'assignmentmode', color: '#7F77DD', dashed: true },
];

function ClassCard({ cls, expanded, onToggle, highlight, onFieldHover }) {
  return (
    <div
      className={styles.card}
      style={{
        borderLeftColor: cls.color,
        boxShadow: highlight ? `0 0 0 1.5px ${cls.color}` : undefined,
        transition: 'box-shadow 0.2s',
      }}
    >
      <div className={styles.cardHeader} onClick={onToggle} style={{ cursor: 'pointer' }}>
        <div className={styles.cardTitle}>
          <span className={styles.className}>{cls.name}</span>
          <span className={styles.classKind}>{cls.kind}</span>
        </div>
        <div className={styles.cardMeta}>
          <span className={styles.desc}>{cls.description}</span>
          <span className={styles.toggle}>{expanded ? '▲ collapse' : '▼ expand'}</span>
        </div>
      </div>

      {expanded && (
        <div className={styles.cardBody}>
          <table className={styles.fieldTable}>
            <thead>
              <tr>
                <th>field</th>
                <th>type</th>
                <th>default</th>
                <th>description</th>
              </tr>
            </thead>
            <tbody>
              {cls.fields.map(f => (
                <tr
                  key={f.name}
                  onMouseEnter={() => f.isRef && onFieldHover(f.refId)}
                  onMouseLeave={() => onFieldHover(null)}
                  style={{ background: f.isRef ? `${cls.color}10` : undefined }}
                >
                  <td>
                    <code className={styles.fieldName} style={f.isRef ? { color: f.type === 'assignmentmode' ? '#7F77DD' : '#1D9E75', fontWeight: 500 } : {}}>
                      {f.name}
                    </code>
                  </td>
                  <td><code className={styles.fieldType}>{f.type}</code></td>
                  <td><code className={styles.fieldDefault}>{f.default}</code></td>
                  <td className={styles.fieldDesc}>{f.desc}</td>
                </tr>
              ))}
            </tbody>
          </table>
          {cls.methods && (
            <div className={styles.methods}>
              <div className={styles.methodsLabel}>methods</div>
              {cls.methods.map(m => (
                <div key={m.name} className={styles.method}>
                  <code className={styles.methodName}>{m.name}</code>
                  <span className={styles.methodDesc}>{m.desc}</span>
                </div>
              ))}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export default function ClassDiagram() {
  const [expanded, setExpanded] = useState({ network: true });
  const [highlight, setHighlight] = useState(null);

  const toggle = (id) => setExpanded(e => ({ ...e, [id]: !e[id] }));

  return (
    <div className={styles.wrapper}>
      <div className={styles.legend}>
        {[
          { color: '#378ADD', label: 'main class' },
          { color: '#1D9E75', label: 'subclass (composition)' },
          { color: '#7F77DD', label: 'shared utility class' },
        ].map(l => (
          <span key={l.label} className={styles.legendItem}>
            <span className={styles.legendDot} style={{ background: l.color }} />
            {l.label}
          </span>
        ))}
      </div>

      <div className={styles.connectionNote}>
        <strong>Hover</strong> over a reference field to highlight the linked class.
        Fields shown in <span style={{ color: '#1D9E75', fontWeight: 500 }}>green</span> or{' '}
        <span style={{ color: '#7F77DD', fontWeight: 500 }}>purple</span> are subclass references.
      </div>

      {CLASSES.map(cls => (
        <ClassCard
          key={cls.id}
          cls={cls}
          expanded={!!expanded[cls.id]}
          onToggle={() => toggle(cls.id)}
          highlight={highlight === cls.id}
          onFieldHover={setHighlight}
        />
      ))}

      <div className={styles.compositionNote}>
        <strong>Composition</strong>
        {CONNECTIONS.map(c => (
          <div key={`${c.from}-${c.field}`} className={styles.connRow}>
            <code>{c.from}.{c.field}</code>
            <span className={styles.arrow} style={{ color: c.color, borderBottom: c.dashed ? '1.5px dashed' : '1.5px solid', borderColor: c.color }}>→</span>
            <code style={{ color: c.color }}>{c.to}</code>
          </div>
        ))}
      </div>
    </div>
  );
}
