import React, { useEffect, useRef, useState } from 'react';
import styles from './ClassMap.module.css';

const CLASSES = [
  {
    id: 'network',
    name: 'network',
    kind: 'handle',
    color: '#378ADD',
    style: { left: 220, top: 20, width: 250 },
    fields: [
      { name: 'Nreplicates', type: 'int' },
      { name: 'flags', type: 'struct' },
      { name: 'domain', type: 'struct' },
      { name: 'peratom', type: 'struct' },
      { name: 'defect', type: 'struct' },
      { name: 'pot', type: 'struct' },
      { name: 'architecture', type: '→ architecture', isRef: true, refColor: '#1D9E75', fieldId: 'f-arch' },
      { name: 'perbond', type: '→ bondstyle', isRef: true, refColor: '#1D9E75', fieldId: 'f-pb' },
    ],
  },
  {
    id: 'architecture',
    name: 'architecture',
    kind: '',
    color: '#1D9E75',
    style: { left: 20, top: 290, width: 240 },
    fields: [
      { name: 'geometry', type: 'string' },
      { name: 'types', type: 'struct' },
      { name: 'lattice_spacing', type: 'double' },
      { name: 'lattice_disorder_level', type: 'double' },
      { name: 'rho_atom', type: 'double' },
      { name: 'strand_typology', type: '→ assignmentmode', isRef: true, refColor: '#7F77DD', fieldId: 'f-st' },
    ],
  },
  {
    id: 'bondstyle',
    name: 'bondstyle',
    kind: '',
    color: '#1D9E75',
    style: { left: 430, top: 290, width: 200 },
    fields: [
      { name: 'type', type: 'int' },
      { name: 'kuhn', type: '→ assignmentmode', isRef: true, refColor: '#7F77DD', fieldId: 'f-kuhn' },
    ],
  },
  {
    id: 'assignmentmode',
    name: 'assignmentmode',
    kind: 'shared utility',
    color: '#7F77DD',
    style: { left: 170, top: 510, width: 320 },
    fields: [
      { name: 'auto', type: 'boolean' },
      { name: 'mode', type: "string {'mono'|'uniform'|'poly'|'polydisperse'|'bimodal'}" },
      { name: 'mono', type: 'struct {value}' },
      { name: 'uniform', type: 'struct {min_value, max_value}' },
      { name: 'poly', type: 'struct {method, pmf_mean...}' },
      { name: 'bimodal', type: 'struct {mean_1, mean_2, lam_1...}' },
    ],
  },
];

const CONNECTIONS = [
  { from: 'network',      to: 'architecture',   fieldId: 'f-arch',  color: '#1D9E75', dash: false },
  { from: 'network',      to: 'bondstyle',       fieldId: 'f-pb',    color: '#1D9E75', dash: false },
  { from: 'architecture', to: 'assignmentmode',  fieldId: 'f-st',    color: '#7F77DD', dash: true },
  { from: 'bondstyle',    to: 'assignmentmode',  fieldId: 'f-kuhn',  color: '#7F77DD', dash: true },
];

export default function ClassMap() {
  const wrapRef = useRef(null);
  const svgRef = useRef(null);
  const [active, setActive] = useState(null);
  const [highlightedFields, setHighlightedFields] = useState([]);

  function getBox(id) {
    const wrap = wrapRef.current;
    const el = wrap?.querySelector(`[data-cls="${id}"]`);
    if (!el || !wrap) return null;
    const wr = wrap.getBoundingClientRect();
    const er = el.getBoundingClientRect();
    return {
      x: er.left - wr.left,
      y: er.top - wr.top,
      w: er.width,
      h: er.height,
      cx: er.left - wr.left + er.width / 2,
      cy: er.top - wr.top + er.height / 2,
      t: er.top - wr.top,
      b: er.bottom - wr.top,
      l: er.left - wr.left,
      r: er.right - wr.left,
    };
  }

  function drawArrows(activeId) {
    const svg = svgRef.current;
    if (!svg) return;
    while (svg.children.length > 1) svg.removeChild(svg.lastChild);

    CONNECTIONS.forEach(c => {
      const isActive = !activeId || c.from === activeId || c.to === activeId;
      const op = isActive ? 1 : 0.1;

      const f = getBox(c.from);
      const t = getBox(c.to);
      if (!f || !t) return;

      let d = '';
      if (c.from === 'network' && c.to === 'architecture') {
        const px = f.l + 50;
        d = `M${px} ${f.b} L${px} ${f.b + 25} L${t.cx} ${f.b + 25} L${t.cx} ${t.t}`;
      } else if (c.from === 'network' && c.to === 'bondstyle') {
        const px = f.r - 50;
        d = `M${px} ${f.b} L${px} ${f.b + 25} L${t.cx} ${f.b + 25} L${t.cx} ${t.t}`;
      } else if (c.from === 'architecture' && c.to === 'assignmentmode') {
        const px = t.l + 60;
        d = `M${px} ${f.b} L${px} ${t.t}`;
      } else if (c.from === 'bondstyle' && c.to === 'assignmentmode') {
        const px = t.r - 60;
        d = `M${px} ${f.b} L${px} ${t.t}`;
      }

      const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      path.setAttribute('d', d);
      path.setAttribute('fill', 'none');
      path.setAttribute('stroke', c.color);
      path.setAttribute('stroke-width', isActive ? '1.5' : '1');
      path.setAttribute('stroke-opacity', String(op));
      path.setAttribute('marker-end', 'url(#arrowhead)');
      if (c.dash) path.setAttribute('stroke-dasharray', '5 3');
      svg.appendChild(path);
    });
  }

  function handleClick(id) {
    const next = active === id ? null : id;
    setActive(next);
    if (next) {
      const fields = CONNECTIONS
        .filter(c => c.from === next || c.to === next)
        .map(c => c.fieldId);
      setHighlightedFields(fields);
    } else {
      setHighlightedFields([]);
    }
  }

  useEffect(() => {
    const t = setTimeout(() => drawArrows(active), 60);
    return () => clearTimeout(t);
  }, [active]);

  useEffect(() => {
    const obs = new ResizeObserver(() => drawArrows(active));
    if (wrapRef.current) obs.observe(wrapRef.current);
    return () => obs.disconnect();
  }, [active]);

  const connectedClasses = active
    ? new Set([active, ...CONNECTIONS.filter(c => c.from === active || c.to === active).flatMap(c => [c.from, c.to])])
    : null;

  return (
    <div className={styles.outer}>
      <div className={styles.legend}>
        <span className={styles.legendItem}><span className={styles.solid} style={{ background: '#1D9E75' }} />composition</span>
        <span className={styles.legendItem}><span className={styles.dashed} style={{ borderColor: '#7F77DD' }} />shared utility</span>
        <span className={styles.legendItem} style={{ marginLeft: 'auto', opacity: 0.6 }}>click a class to highlight connections</span>
      </div>

      <div className={styles.wrap} ref={wrapRef}>
        {CLASSES.map(cls => {
          const isHighlighted = connectedClasses ? connectedClasses.has(cls.id) : false;
          const isDimmed = connectedClasses ? !connectedClasses.has(cls.id) : false;
          return (
            <div
              key={cls.id}
              data-cls={cls.id}
              className={`${styles.cls} ${isHighlighted ? styles.highlighted : ''} ${isDimmed ? styles.dimmed : ''}`}
              style={{
                left: cls.style.left,
                top: cls.style.top,
                width: cls.style.width,
                borderLeftColor: cls.color,
                ['--hl-col' as string]: cls.color,
              }}
              onClick={() => handleClick(cls.id)}
            >
              <div className={styles.clsHeader}>
                <span className={styles.clsName}>{cls.name}</span>
                {cls.kind && <span className={styles.clsKind}>{cls.kind}</span>}
              </div>
              <div className={styles.clsBody}>
                {cls.fields.map(f => (
                  <div
                    key={f.name}
                    className={`${styles.field} ${f.isRef ? styles.fieldRef : ''}`}
                    style={{
                      color: f.isRef ? f.refColor : undefined,
                      background: f.isRef && highlightedFields.includes(f.fieldId)
                        ? `${f.refColor}22`
                        : undefined,
                    }}
                  >
                    {f.name}
                    <span className={styles.ftype}>{f.type}</span>
                  </div>
                ))}
              </div>
            </div>
          );
        })}

        <svg ref={svgRef} className={styles.svg}>
          <defs>
            <marker id="arrowhead" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">
              <path d="M2 1L8 5L2 9" fill="none" stroke="context-stroke" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" />
            </marker>
          </defs>
        </svg>
      </div>
    </div>
  );
}
