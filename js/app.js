/*PHYS_START*/
const PHYS = (() => {
  const R = 8.314, SIG = 5.670e-8, T0 = 298.15;
  const RHO = 7874, QFE = 6.7e6, LF = 2.47e5, TM = 1811.15, CP_OX = 800, OXR = 1.382;
  const ETA = 0.6;
  const A_TUBE = Math.PI * 0.025 * 0.025, L_BASE = 0.03, R_CONTACT = 0.08;
  const EPS = 0.8, H_CONV = 15;
  const T_FLAME = 3773.15, T_CAP = 2673.15;
  const A1 = 2.07e-4, E1 = 40000, S0 = 1e-9;          // low-temperature thin-film growth
  const EA = 170000, KP0 = 2.84e-4, SMIN = 5e-9;       // parabolic (diffusion through scale)
  const T_W = 843.15, T_FEO = 1650.15;                 // 570 °C wüstite, 1377 °C FeO melts
  const KSUP = 0.08, ACAP = 0.04;
  const S_NATIVE = 3e-9;                               // native oxide already on the metal                      // oxygen supply limit (mass transfer)

  function create(form, massG) {
    const m0 = massG / 1000;
    let r, Atot, Arad;
    if (form.d) { r = form.d / 2; Atot = 3 * m0 / (RHO * r); Arad = A_TUBE * 1.5; }
    else { r = Math.cbrt(3 * m0 / (4 * Math.PI * RHO)); Atot = 4 * Math.PI * r * r; Arad = 0.85 * Atot; }
    return { m0, r, Atot, Arad, T: T0, s: Math.min(S_NATIVE, r * 0.05), X: 0, fL: 0, maxFL: 0, time: 0,
      Tmax: T0, ignited: false, ignT: null, ignTime: null, melted: false, meltTime: null,
      dTdt: 0, dXdt: 0, dfLdt: 0, ignAcc: 0, ignCand: null };
  }
  const cpFe = T => T >= TM ? 820 : Math.min(450 + 0.3 * (T - T0), 700);
  function heatCap(st) { return st.m0 * (1 - st.X) * cpFe(st.T) + st.m0 * st.X * OXR * CP_OX + 1e-6; }
  function baseG(k) { return k * A_TUBE / L_BASE; }
  function conductance(k) { return 1 / (1 / baseG(k) + R_CONTACT); }
  function baseTopT(k, qCond) { return T0 + qCond / baseG(k); }
  function losses(st, p, T) {
    const cond = conductance(p.k) * (T - T0);
    const rad = EPS * SIG * st.Arad * (T ** 4 - T0 ** 4);
    const conv = H_CONV * st.Arad * (T - T0);
    return { cond, rad, conv, total: cond + rad + conv };
  }
  function torchIn(p, T) { return (p.torchOn ? p.P : 0) * Math.max(0, (T_FLAME - T) / (T_FLAME - T0)); }
  const ZERO = { v: 0, mdot: 0, Q: 0, limited: false };
  function reaction(st, p, T, s) {
    const x = p.xO2;
    if (x <= 0 || s >= st.r * 0.999999) return ZERO;
    const fO = Math.sqrt(x / 0.21);
    const v1 = A1 * Math.exp(-E1 / (R * T)) * Math.exp(-s / S0) * fO;
    const fw = 0.02 + 0.98 / (1 + Math.exp(-(T - T_W) / 12));
    const prot = 1 / (1 + Math.exp((T - T_FEO) / 25));
    const Kp = KP0 * Math.exp(-EA / (R * T)) * fw * fO;
    const v2 = Kp / (s * prot + SMIN);
    const core = 1 - s / st.r, core2 = core * core;
    const mkin = RHO * st.Atot * core2 * (v1 + v2);
    const fsup = Math.log(1 / (1 - Math.min(x, 0.99)));  // Stefan flow: no inert gas to diffuse through in pure O2
    const msup = KSUP * fsup * Math.min(st.Atot, ACAP) * core2;
    const mdot = mkin * msup / (mkin + msup);
    return { v: mdot / (RHO * st.Atot * core2), mdot, Q: mdot * QFE, limited: mkin > msup };
  }
  // Oxide depth the sample would reach after lingering ~tau seconds at T (used for the heat-balance chart).
  function sQuasi(st, p, T, tau) {
    const x = p.xO2;
    if (x <= 0) return st.s;
    const fO = Math.sqrt(x / 0.21);
    const sl = S0 * Math.log(1 + A1 * Math.exp(-E1 / (R * T)) * fO * tau / S0);
    const fw = 0.02 + 0.98 / (1 + Math.exp(-(T - T_W) / 12));
    const sp = Math.sqrt(2 * KP0 * Math.exp(-EA / (R * T)) * fw * fO * tau);
    return Math.max(st.s, Math.min(Math.max(sl, sp), st.r * 0.1));
  }
  function applyHeat(st, dE) {
    const C = heatCap(st);
    const mFe = st.m0 * (1 - st.X);
    const canMelt = mFe > 1e-7;
    if (dE >= 0) {
      if (st.T < TM) {
        const need = (TM - st.T) * C;
        if (dE < need || !canMelt) { st.T += dE / C; return; }
        st.T = TM; dE -= need;
      }
      if (canMelt && st.fL < 1 && st.T <= TM + 1e-6) {
        const L = (1 - st.fL) * mFe * LF;
        if (dE < L) { st.fL += dE / (mFe * LF); return; }
        st.fL = 1; dE -= L;
      }
      st.T += dE / C;
    } else {
      let q = -dE;
      if (st.T > TM + 1e-6) {
        const avail = (st.T - TM) * C;
        if (q < avail) { st.T -= q / C; return; }
        st.T = TM; q -= avail;
      }
      if (canMelt && st.fL > 0) {
        const L = st.fL * mFe * LF;
        if (q < L) { st.fL -= q / (mFe * LF); return; }
        st.fL = 0; q -= L;
      }
      st.T -= q / C;
    }
  }
  function deriv(st, p) {
    const rx = reaction(st, p, st.T, st.s);
    const lo = losses(st, p, st.T);
    const Pin = torchIn(p, st.T);
    const gen = ETA * rx.Q;
    return { rx, lo, Pin, gen, net: Pin + gen - lo.total };
  }
  // Advance up to dtTotal seconds; stops early if the per-call change limits are hit (auto slow motion).
  function step(st, p, dtTotal, lim) {
    lim = lim || { dT: Infinity, dX: Infinity, dfL: Infinity };
    const T_a = st.T, X_a = st.X, fL_a = st.fL;
    let t = 0, n = 0, d = deriv(st, p);
    while (t < dtTotal && n < 8000) {
      n++;
      const C = heatCap(st);
      const dTdt = d.net / C;
      let dt = Math.min(dtTotal - t, 0.5);
      if (Math.abs(dTdt) > 1e-12) dt = Math.min(dt, 1.5 / Math.abs(dTdt));
      if (d.rx.v > 0) {
        dt = Math.min(dt, (0.05 * st.s + 2e-10) / d.rx.v);
        const core = 1 - st.s / st.r;
        const dXdt = 3 * core * core * d.rx.v / st.r;
        if (dXdt > 0) dt = Math.min(dt, 0.002 / dXdt);
      }
      dt = Math.max(dt, 1e-6);
      st.s = Math.min(st.r, st.s + d.rx.v * dt);
      st.X = 1 - Math.pow(1 - st.s / st.r, 3);
      applyHeat(st, d.net * dt);
      if (st.T > T_CAP) st.T = T_CAP;
      st.time += dt; t += dt;
      d = deriv(st, p);
      if (st.T > st.Tmax) st.Tmax = st.T;
      if (!st.ignited) {
        // self-sustaining runaway: reaction heat alone beats all losses and keeps heating fast, for >1 s
        if (d.gen > 150 && d.gen > d.lo.total && d.net / heatCap(st) > 2) {
          if (!st.ignCand) st.ignCand = { T: st.T, t: st.time };
          st.ignAcc += dt;
          if (st.ignAcc > 1) { st.ignited = true; st.ignT = st.ignCand.T; st.ignTime = st.ignCand.t; }
        } else { st.ignAcc = 0; st.ignCand = null; }
      }
      if (st.fL > st.maxFL) st.maxFL = st.fL;
      if (!st.melted && (st.fL > 0.02 || (st.X > 0.9 && st.T > 1870))) { st.melted = true; st.meltTime = st.time; }
      if (Math.abs(st.T - T_a) > lim.dT || st.X - X_a > lim.dX || Math.abs(st.fL - fL_a) > lim.dfL) break;
    }
    if (t > 0) { st.dTdt = (st.T - T_a) / t; st.dXdt = (st.X - X_a) / t; st.dfLdt = Math.abs(st.fL - fL_a) / t; }
    st.last = d;
    return t;
  }
  return { create, reaction, sQuasi, losses, torchIn, heatCap, step, conductance, baseTopT, ETA, T0, TM, T_CAP };
})();
/*PHYS_END*/

(() => {
'use strict';
const ETA = PHYS.ETA, T0 = PHYS.T0;
const $ = id => document.getElementById(id);
const FONT = '"IBM Plex Sans KR", system-ui, -apple-system, "Apple SD Gothic Neo", "Malgun Gothic", sans-serif';

const FORMS = {
  nano:   { label: '나노분말', sub: '0.1 μm', d: 0.1e-6, col: '#2b2d30', rad: [0.9, 1.5] },
  fine:   { label: '미세분말', sub: '5 μm',   d: 5e-6,   col: '#56595e', rad: [1.1, 1.9] },
  powder: { label: '분말',     sub: '50 μm',  d: 50e-6,  col: '#787d83', rad: [1.4, 2.4] },
  coarse: { label: '굵은 분말', sub: '0.5 mm', d: 500e-6, col: '#8a8f95', rad: [2.2, 3.4] },
  bulk:   { label: '덩어리',   sub: '한 조각', d: 0,      col: '#9da3aa' },
};
const METALS = [
  { name: '스테인리스', k: 16,  mp: 1450, c1: '#9aa1a8', c2: '#596068' },
  { name: '순철',       k: 80,  mp: 1538, c1: '#7b8086', c2: '#45494e' },
  { name: '황동',       k: 110, mp: 930,  c1: '#d8b356', c2: '#8a6a1f' },
  { name: '알루미늄',   k: 237, mp: 660,  c1: '#d6dbe0', c2: '#8e959d' },
  { name: '구리',       k: 401, mp: 1085, c1: '#de8656', c2: '#8f431c' },
];
const K_MIN = 5, K_MAX = 420;
const PRESETS = [
  { label: '느린 산화', form: 'powder', mass: 100, P: 500,  o2: 21,  k: 16 },
  { label: '분말 점화', form: 'powder', mass: 100, P: 1000, o2: 21,  k: 16 },
  { label: '같은 토치, 구리 받침', form: 'powder', mass: 100, P: 1000, o2: 21, k: 401 },
  { label: '덩어리 용융', form: 'bulk', mass: 100, P: 5000, o2: 21, k: 16 },
  { label: '순산소 속 덩어리', form: 'bulk', mass: 100, P: 3200, o2: 100, k: 16 },
  { label: '나노분말 자연발화', form: 'nano', mass: 100, P: 0, o2: 21, k: 16 },
];

const uFromP = P => Math.round(1000 * Math.sqrt(P / 10000));
const PFromU = u => Math.round(10000 * (u / 1000) ** 2 / 10) * 10;
const uFromK = k => Math.round(1000 * Math.log(k / K_MIN) / Math.log(K_MAX / K_MIN));
const kFromU = u => K_MIN * Math.pow(K_MAX / K_MIN, u / 1000);
const nearestMetal = k => METALS.reduce((a, b) => Math.abs(Math.log(b.k / k)) < Math.abs(Math.log(a.k / k)) ? b : a);

const fmtW = w => w < 1 ? '<1 W' : w < 1000 ? `${Math.round(w)} W` : `${(w / 1000).toFixed(w < 10000 ? 1 : 0)} kW`;
const fmtLen = nm => nm < 1000 ? `${Math.round(nm)} nm` : nm < 1e6 ? `${(nm / 1000).toFixed(nm < 1e4 ? 1 : 0)} μm` : `${(nm / 1e6).toFixed(1)} mm`;
const fmtTime = s => { s = Math.floor(s); const m = Math.floor(s / 60); return `${m}:${String(s % 60).padStart(2, '0')}`; };
const C = T => T - 273.15;
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));

// ---------- state ----------
let formKey = 'powder', massG = 100;
const p = { P: 500, xO2: 0.21, k: 16, torchOn: true };
let st = PHYS.create(FORMS[formKey], massG);
let status = 'idle';         // idle | running | paused | done
let speed = 10, effSpeed = 0;
let hist = [], lastPushT = -1;

// ---------- controls ----------
const rP = $('rP'), rO = $('rO'), rK = $('rK'), rM = $('rM');
function syncLabels() {
  $('oP').textContent = p.P === 0 ? '꺼짐 (0 W)' : fmtW(p.P);
  $('oO').textContent = `${Math.round(p.xO2 * 100)} %` + (Math.round(p.xO2 * 100) === 21 ? ' (공기)' : '');
  const m = nearestMetal(p.k);
  const exact = Math.abs(Math.log(m.k / p.k)) < 0.04;
  $('oK').textContent = `${p.k < 20 ? p.k.toFixed(1) : Math.round(p.k)} W/m·K`;
  document.querySelectorAll('#metals .chip').forEach(b => b.setAttribute('aria-pressed', String(exact && b.dataset.k == m.k)));
  $('oM').textContent = `${massG} g`;
  document.querySelectorAll('#forms .chip').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.f === formKey)));
  const locked = status !== 'idle';
  document.querySelectorAll('#forms .chip').forEach(b => b.disabled = locked);
  rM.disabled = locked;
  $('lockHint').textContent = locked ? '리셋 후 바꿀 수 있어요' : '';
}
function setControlsFromP() {
  rP.value = uFromP(p.P); rO.value = Math.round(p.xO2 * 100); rK.value = uFromK(p.k); rM.value = massG;
  syncLabels();
}
rP.addEventListener('input', () => { p.P = PFromU(+rP.value); syncLabels(); });
rO.addEventListener('input', () => { p.xO2 = +rO.value / 100; syncLabels(); });
rK.addEventListener('input', () => { p.k = kFromU(+rK.value); syncLabels(); });
rM.addEventListener('input', () => { if (status !== 'idle') return; massG = +rM.value; rebuildSample(); syncLabels(); });

function buildChips() {
  const fBox = $('forms');
  for (const [key, f] of Object.entries(FORMS)) {
    const b = document.createElement('button');
    b.className = 'chip'; b.dataset.f = key; b.type = 'button';
    b.innerHTML = `${f.label}<small>${f.sub}</small>`;
    b.addEventListener('click', () => { if (status !== 'idle') return; formKey = key; rebuildSample(); syncLabels(); });
    fBox.appendChild(b);
  }
  const mBox = $('metals');
  for (const m of METALS) {
    const b = document.createElement('button');
    b.className = 'chip'; b.dataset.k = m.k; b.type = 'button';
    b.textContent = `${m.name} ${m.k}`;
    b.addEventListener('click', () => { p.k = m.k; rK.value = uFromK(m.k); syncLabels(); });
    mBox.appendChild(b);
  }
  const pBox = $('presets');
  for (const pr of PRESETS) {
    const b = document.createElement('button');
    b.className = 'chip'; b.type = 'button'; b.textContent = pr.label;
    b.addEventListener('click', () => {
      reset();
      formKey = pr.form; massG = pr.mass; p.P = pr.P; p.xO2 = pr.o2 / 100; p.k = pr.k; p.torchOn = true;
      rebuildSample(); setControlsFromP();
    });
    pBox.appendChild(b);
  }
}

function rebuildSample() { st = PHYS.create(FORMS[formKey], massG); hist = []; lastPushT = -1; buildBed(); }

const btnRun = $('btnRun'), btnTorch = $('btnTorch');
function syncButtons() {
  btnRun.textContent = status === 'idle' ? '스위치 켜기' : status === 'running' ? '일시정지' : status === 'paused' ? '계속' : '다시 실행';
  btnRun.classList.toggle('on', status === 'running');
  btnTorch.disabled = status === 'idle';
  btnTorch.textContent = p.torchOn ? '토치 끄기' : '토치 켜기';
  syncLabels();
}
btnRun.addEventListener('click', () => {
  if (status === 'idle') { p.torchOn = true; status = 'running'; pushHist(true); }
  else if (status === 'running') status = 'paused';
  else if (status === 'paused') status = 'running';
  else { reset(); p.torchOn = true; status = 'running'; pushHist(true); }
  syncButtons();
});
btnTorch.addEventListener('click', () => { p.torchOn = !p.torchOn; syncButtons(); });
$('btnReset').addEventListener('click', () => { reset(); });
function reset() {
  status = 'idle'; p.torchOn = true; rebuildSample();
  sparks.length = 0; smoke.length = 0;
  syncButtons();
}
document.querySelectorAll('#speeds button').forEach(b => b.addEventListener('click', () => {
  speed = +b.dataset.s;
  document.querySelectorAll('#speeds button').forEach(x => x.setAttribute('aria-pressed', String(x === b)));
}));

function pushHist(force) {
  const minDt = Math.max(0.05, st.time / 900);
  if (!force && st.time - lastPushT < minDt) return;
  hist.push({ t: st.time, T: C(st.T), X: st.X });
  lastPushT = st.time;
}

// ---------- stage drawing ----------
const cv = $('stage'), g = cv.getContext('2d'), stageBox = $('stageBox');
const VW = 360, VH = 540;
const GX0 = 118, GX1 = 242, GT = 70, BY = 440, BASE_B = 488, SINK_B = 532, CX = 180, TIP_Y = 318, PIPE_Y = 350;
let dpr = 1, sScale = 1;
function resizeStage() {
  dpr = Math.min(window.devicePixelRatio || 1, 2.5);
  const cw = stageBox.clientWidth || VW;
  const maxH = window.innerWidth < 860 ? window.innerHeight * 0.6 : 1e9;
  const w = Math.floor(Math.min(cw, maxH * VW / VH, 460));
  sScale = w / VW;
  cv.style.width = w + 'px'; cv.style.height = Math.round(w * VH / VW) + 'px';
  cv.width = Math.round(w * dpr); cv.height = Math.round(w * VH / VW * dpr);
}

const hex = h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
const mix = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
const rgba = (c, a = 1) => `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${a})`;
const shade = (c, f) => f >= 0 ? mix(c, [255, 255, 255], f) : mix(c, [0, 0, 0], -f);
function interp(stops, x) {
  if (x <= stops[0][0]) return stops[0][1];
  for (let i = 1; i < stops.length; i++) if (x <= stops[i][0]) {
    const [x0, c0] = stops[i - 1], [x1, c1] = stops[i];
    return mix(c0, c1, (x - x0) / (x1 - x0));
  }
  return stops[stops.length - 1][1];
}
// interference (temper) colours by oxide thickness in nm
const OX_TAIL = [[30, '#d8be78'], [42, '#cf9a4b'], [52, '#a35e3a'], [62, '#7b4a8e'], [72, '#3f60ae'], [86, '#7095be'], [130, '#6a737c'], [450, '#3b3c3f'], [5000, '#2c2826']].map(([x, c]) => [x, hex(c)]);
function sampleColor() {
  const base = hex(FORMS[formKey].col);
  const dox = 2.1 * st.s * 1e9;
  let c = interp([[0, base], [20, base], ...OX_TAIL], dox);
  if (st.X > 0.25) c = mix(c, hex('#4d2b22'), Math.min(1, (st.X - 0.25) / 0.75) * 0.7);
  return c;
}
const GLOW = [[460, [110, 0, 0], 0], [540, [150, 14, 0], 0.32], [680, [205, 36, 6], 0.55], [820, [240, 86, 12], 0.72], [980, [255, 135, 35], 0.83], [1180, [255, 186, 82], 0.9], [1400, [255, 222, 150], 0.95], [1700, [255, 246, 226], 1]];
function glowAt(Tc) {
  if (Tc <= GLOW[0][0]) return null;
  for (let i = 1; i < GLOW.length; i++) if (Tc <= GLOW[i][0]) {
    const a = GLOW[i - 1], b = GLOW[i], t = (Tc - a[0]) / (b[0] - a[0]);
    return { c: mix(a[1], b[1], t), a: a[2] + (b[2] - a[2]) * t };
  }
  return { c: GLOW[GLOW.length - 1][1], a: 1 };
}

let bed = [], bedH = 40;
function rng(seed) { let s = seed >>> 0; return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }; }
function bedTop(x) { const u = (x - GX0) / (GX1 - GX0); return BY - bedH * (0.86 + 0.14 * Math.sin(Math.PI * u)); }
function buildBed() {
  bedH = Math.min(74, 16 + 0.28 * massG);
  bed = [];
  if (formKey === 'bulk') return;
  const r = rng(11), rad = FORMS[formKey].rad;
  let tries = 0;
  while (bed.length < 360 && tries < 9000) {
    tries++;
    const x = GX0 + 2 + r() * (GX1 - GX0 - 4), y = BY - 2 - r() * bedH;
    const pr = rad[0] + r() * (rad[1] - rad[0]);
    if (y - pr * 0.6 < bedTop(x)) continue;
    bed.push({ x, y, r: pr, sh: (r() - 0.5) * 0.36, seed: r() * 1000 });
  }
  bed.sort((a, b) => a.y - b.y);
}
function bulkDims() { const f = Math.cbrt(massG / 100); return { w: 62 * f, h: 40 * f }; }
function sampleShape() {
  const m = st.maxFL;
  if (formKey !== 'bulk') return { x0: GX0, w: GX1 - GX0, h: bedH * 0.55, m };
  const { w, h } = bulkDims(), pw = GX1 - GX0, ph = w * h / pw * 0.9;
  const W = w + (pw - w) * m, H = h + (ph - h) * m;
  return { x0: CX - W / 2, w: W, h: H, m };
}
function topAt(x) {
  if (formKey !== 'bulk') {
    const pool = BY - bedH * 0.55;
    return bedTop(x) + (pool - bedTop(x)) * st.maxFL;
  }
  const s = sampleShape();
  return BY - s.h;
}
function rrect(x, y, w, h, r) {
  r = Math.max(0, Math.min(r, w / 2, h / 2));
  g.beginPath();
  g.moveTo(x + r, y); g.lineTo(x + w - r, y); g.quadraticCurveTo(x + w, y, x + w, y + r);
  g.lineTo(x + w, y + h - r); g.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
  g.lineTo(x + r, y + h); g.quadraticCurveTo(x, y + h, x, y + h - r);
  g.lineTo(x, y + r); g.quadraticCurveTo(x, y, x + r, y); g.closePath();
}
function blob(x0, w, h, col, liquid, t, alpha, rTop) {
  g.save(); g.globalAlpha = alpha;
  const y0 = BY - h;
  g.beginPath();
  g.moveTo(x0, BY);
  g.lineTo(x0, y0 + rTop);
  g.quadraticCurveTo(x0, y0, x0 + rTop, y0);
  const steps = 18;
  for (let i = 0; i <= steps; i++) {
    const x = x0 + rTop + (w - 2 * rTop) * i / steps;
    const wob = liquid > 0 ? Math.sin(x * 0.18 + t * 5) * 0.9 * liquid : 0;
    g.lineTo(x, y0 + wob);
  }
  g.quadraticCurveTo(x0 + w, y0, x0 + w, y0 + rTop);
  g.lineTo(x0 + w, BY); g.closePath();
  const gr = g.createLinearGradient(x0, y0, x0 + w, BY);
  gr.addColorStop(0, rgba(shade(col, 0.22))); gr.addColorStop(0.55, rgba(col)); gr.addColorStop(1, rgba(shade(col, -0.38)));
  g.fillStyle = gr; g.fill();
  if (liquid > 0.02) {
    g.strokeStyle = `rgba(255,255,255,${0.5 * liquid})`; g.lineWidth = 1.3;
    g.beginPath();
    for (let i = 0; i <= steps; i++) {
      const x = x0 + rTop + 4 + (w - 2 * rTop - 8) * i / steps;
      const y = y0 + 3 + Math.sin(x * 0.18 + t * 5) * 0.9 * liquid;
      i ? g.lineTo(x, y) : g.moveTo(x, y);
    }
    g.stroke();
  } else {
    g.fillStyle = 'rgba(255,255,255,0.14)';
    g.fillRect(x0 + 6, y0 + 4, Math.max(0, w * 0.35), 2);
  }
  g.restore();
}

function drawSample(t, col, gl, burnI) {
  const body = gl ? mix(col, gl.c, gl.a) : col;
  const s = sampleShape();
  if (formKey !== 'bulk') {
    if (s.m < 1) {
      g.save(); g.globalAlpha = 1 - s.m;
      g.beginPath(); g.moveTo(GX0, BY);
      for (let x = GX0; x <= GX1; x += 4) g.lineTo(x, bedTop(x) + 1.5);
      g.lineTo(GX1, BY); g.closePath();
      g.fillStyle = rgba(shade(body, -0.4)); g.fill();
      for (const q of bed) {
        let c = shade(body, q.sh);
        if (burnI > 0) {
          const e = Math.pow(0.5 + 0.5 * Math.sin(q.seed * 12.9898 + t * (5 + (q.seed % 5))), 6) * burnI * 1.3;
          if (e > 0.02) c = mix(c, [255, 228, 165], Math.min(1, e));
        }
        g.fillStyle = rgba(c);
        g.beginPath(); g.arc(q.x, q.y, q.r, 0, 6.2832); g.fill();
      }
      g.restore();
    }
    if (s.m > 0) blob(GX0, GX1 - GX0, s.h, body, st.fL, t, s.m, 6);
  } else {
    const rTop = 7 + (s.h / 2 - 7) * Math.min(1, s.m) * 0.6;
    blob(s.x0, s.w, s.h, body, st.fL, t, 1, rTop);
    if (burnI > 0) {
      g.save(); g.globalCompositeOperation = 'lighter';
      for (let i = 0; i < 10; i++) {
        const x = s.x0 + 6 + ((i * 37.7) % (s.w - 12));
        const e = Math.pow(0.5 + 0.5 * Math.sin(i * 7.1 + t * 6), 5) * burnI;
        g.fillStyle = `rgba(255,225,160,${e})`;
        g.beginPath(); g.arc(x, BY - s.h + 2, 2.2, 0, 6.2832); g.fill();
      }
      g.restore();
    }
  }
}

function drawBaseAndSink(t, metal, cond, overheated) {
  // cooling water
  rrect(70, BASE_B, 220, SINK_B - BASE_B, 8);
  let gr = g.createLinearGradient(0, BASE_B, 0, SINK_B);
  gr.addColorStop(0, '#14456b'); gr.addColorStop(1, '#0b2438');
  g.fillStyle = gr; g.fill();
  g.strokeStyle = 'rgba(150,205,255,0.22)'; g.lineWidth = 1;
  for (let k = 0; k < 2; k++) {
    g.beginPath();
    for (let x = 76; x <= 284; x += 4) { const y = BASE_B + 12 + k * 9 + Math.sin(x * 0.07 + t * 1.6 + k) * 1.6; x === 76 ? g.moveTo(x, y) : g.lineTo(x, y); }
    g.stroke();
  }
  g.fillStyle = 'rgba(170,215,255,0.75)'; g.font = `500 10px ${FONT}`; g.textAlign = 'center'; g.textBaseline = 'alphabetic';
  g.fillText('냉각수 25 °C', CX, SINK_B - 7);
  // metal base
  const c1 = hex(metal.c1), c2 = hex(metal.c2);
  gr = g.createLinearGradient(96, 0, 264, 0);
  gr.addColorStop(0, rgba(c2)); gr.addColorStop(0.45, rgba(c1)); gr.addColorStop(1, rgba(c2));
  g.fillStyle = gr; g.fillRect(96, BY, 168, BASE_B - BY);
  g.fillStyle = 'rgba(255,255,255,0.22)'; g.fillRect(96, BY, 168, 1.5);
  if (overheated) {
    const og = g.createLinearGradient(0, BY, 0, BY + 14);
    og.addColorStop(0, 'rgba(255,90,40,0.85)'); og.addColorStop(1, 'rgba(255,90,40,0)');
    g.fillStyle = og; g.fillRect(96, BY, 168, 14);
  }
  // heat flow chevrons
  const a = clamp(Math.log10(1 + cond) / 3.6, 0, 1);
  if (status !== 'idle' && a > 0.05) {
    g.strokeStyle = `rgba(255,120,60,${0.2 + 0.65 * a})`; g.lineWidth = 1.4 + 1.8 * a; g.lineCap = 'round';
    for (const x of [109, 251]) for (let k = 0; k < 2; k++) {
      const y = BY + 6 + ((t * 14 * (0.5 + a) + k * 16) % 30);
      g.beginPath(); g.moveTo(x - 5, y); g.lineTo(x, y + 5); g.lineTo(x + 5, y); g.stroke();
    }
  }
}
function drawBaseLabel(metal, overheated) {
  const c1 = hex(metal.c1);
  const light = (c1[0] + c1[1] + c1[2]) / 3 > 150;
  g.save();
  g.shadowColor = light ? 'rgba(255,255,255,0.6)' : 'rgba(0,0,0,0.6)'; g.shadowBlur = 3;
  g.fillStyle = light ? 'rgba(16,20,24,0.9)' : 'rgba(255,255,255,0.92)';
  g.font = `700 11px ${FONT}`; g.textAlign = 'center';
  g.fillText(`${metal.name} 받침`, CX, BY + 21);
  g.font = `500 10px ${FONT}`;
  if (overheated) g.fillStyle = light ? '#8a1c00' : '#ffd2c2';
  g.fillText(overheated ? `표면이 녹는점 ${metal.mp} °C 초과` : `k = ${p.k < 20 ? p.k.toFixed(1) : Math.round(p.k)} W/m·K`, CX, BY + 36);
  g.restore();
}

function drawTorch(on) {
  g.strokeStyle = '#2f2623'; g.lineWidth = 7; g.lineCap = 'round';
  g.beginPath(); g.moveTo(CX + 6, 10); g.bezierCurveTo(CX + 20, -4, CX + 50, 2, CX + 64, -8); g.stroke();
  rrect(CX - 20, 6, 40, 20, 5);
  let gr = g.createLinearGradient(0, 6, 0, 26);
  gr.addColorStop(0, '#c9a35a'); gr.addColorStop(1, '#6f5524');
  g.fillStyle = gr; g.fill();
  g.fillStyle = on ? '#ff7a32' : '#3a3f45';
  g.beginPath(); g.arc(CX - 10, 16, 3, 0, 6.2832); g.fill();
  gr = g.createLinearGradient(CX - 4, 0, CX + 4, 0);
  gr.addColorStop(0, '#5b6067'); gr.addColorStop(0.5, '#cfd4da'); gr.addColorStop(1, '#5f646b');
  g.fillStyle = gr; g.fillRect(CX - 3.5, 26, 7, TIP_Y - 34);
  rrect(CX - 5, TIP_Y - 10, 10, 10, 2); g.fillStyle = '#a7803a'; g.fill();
  g.fillStyle = on ? '#ffcf9e' : '#7d8992'; g.font = `600 11px ${FONT}`; g.textAlign = 'left';
  g.fillText(on ? `토치 ${fmtW(p.P)}` : (status === 'idle' ? `토치 ${p.P ? fmtW(p.P) : '꺼짐'}` : '토치 꺼짐'), CX + 26, 21);
}

function drawFlame(t, topY) {
  const pw = Math.sqrt(p.P / 10000);
  const gap = topY - TIP_Y;
  if (gap < 4) return;
  const len = gap * Math.min(1, 0.3 + 0.7 * Math.sqrt(p.P / 1800));
  const fl = 1 + 0.06 * Math.sin(t * 37) + 0.04 * Math.sin(t * 23 + 1);
  g.save(); g.globalCompositeOperation = 'lighter';
  const w = 4 + 15 * pw;
  let gr = g.createLinearGradient(0, TIP_Y, 0, TIP_Y + len);
  gr.addColorStop(0, 'rgba(120,170,255,0.6)'); gr.addColorStop(0.55, 'rgba(90,120,255,0.3)'); gr.addColorStop(1, 'rgba(255,150,70,0.18)');
  g.fillStyle = gr;
  const L = Math.min(gap, len * fl);
  g.beginPath(); g.moveTo(CX - 3, TIP_Y);
  g.quadraticCurveTo(CX - w * fl, TIP_Y + L * 0.6, CX - w * 0.7, TIP_Y + L);
  g.lineTo(CX + w * 0.7, TIP_Y + L);
  g.quadraticCurveTo(CX + w * fl, TIP_Y + L * 0.6, CX + 3, TIP_Y); g.closePath(); g.fill();
  const ic = len * (0.24 + 0.22 * pw);
  gr = g.createLinearGradient(0, TIP_Y, 0, TIP_Y + ic);
  gr.addColorStop(0, 'rgba(215,236,255,0.95)'); gr.addColorStop(1, 'rgba(80,140,255,0.15)');
  g.fillStyle = gr;
  g.beginPath(); g.moveTo(CX - 2.6, TIP_Y); g.quadraticCurveTo(CX, TIP_Y + ic * 1.3, CX + 2.6, TIP_Y); g.fill();
  if (len >= gap * 0.98) {
    const sw = 16 + 48 * pw;
    const rg = g.createRadialGradient(CX, topY - 2, 1, CX, topY - 2, sw);
    rg.addColorStop(0, 'rgba(255,200,140,0.55)'); rg.addColorStop(1, 'rgba(255,120,40,0)');
    g.fillStyle = rg; g.beginPath(); g.ellipse(CX, topY - 2, sw, 6 + 6 * pw, 0, 0, 6.2832); g.fill();
  }
  g.restore();
}

function drawGlass() {
  const wall = x => {
    const gr = g.createLinearGradient(x, 0, x + 6, 0);
    gr.addColorStop(0, 'rgba(190,220,255,0.10)'); gr.addColorStop(0.5, 'rgba(220,240,255,0.36)'); gr.addColorStop(1, 'rgba(190,220,255,0.10)');
    g.fillStyle = gr; g.fillRect(x, GT, 6, BY - GT);
  };
  wall(GX0 - 6); wall(GX1);
  g.fillStyle = 'rgba(255,255,255,0.08)'; g.fillRect(GX0 + 5, GT + 12, 3, BY - GT - 70);
  g.fillStyle = '#2a3139';
  g.fillRect(GX0 - 9, GT - 8, (CX - 12) - (GX0 - 9), 8);
  g.fillRect(CX + 12, GT - 8, (GX1 + 9) - (CX + 12), 8);
  g.fillStyle = '#7d8992'; g.font = `500 10px ${FONT}`; g.textAlign = 'left';
  g.fillText('배기', GX1 + 14, GT - 1);
}

function drawPipe() {
  g.fillStyle = '#36404a'; g.fillRect(0, PIPE_Y - 5, GX0 - 2, 10);
  g.fillStyle = 'rgba(255,255,255,0.12)'; g.fillRect(0, PIPE_Y - 4, GX0 - 2, 2);
  g.fillStyle = '#4b5561'; g.beginPath(); g.arc(46, PIPE_Y, 10, 0, 6.2832); g.fill();
  g.strokeStyle = '#8fd3ff'; g.lineWidth = 2.2; g.lineCap = 'round';
  const a = -Math.PI * 0.95 + p.xO2 * Math.PI * 0.9;
  g.beginPath(); g.moveTo(46, PIPE_Y); g.lineTo(46 + Math.cos(a) * 7.5, PIPE_Y + Math.sin(a) * 7.5); g.stroke();
  g.fillStyle = '#9fdcff'; g.font = `700 12px ${FONT}`; g.textAlign = 'center';
  g.fillText(`O₂ ${Math.round(p.xO2 * 100)}%`, 46, PIPE_Y - 17);
}

function drawThermocouple() {
  g.strokeStyle = '#b89a5a'; g.lineWidth = 1.3; g.lineCap = 'round';
  g.beginPath(); g.moveTo(VW + 2, 290); g.bezierCurveTo(300, 300, 290, 400, GX1 + 6, 418);
  g.lineTo(CX + 24, BY - 8); g.stroke();
  g.fillStyle = '#e2c27b'; g.beginPath(); g.arc(CX + 24, BY - 8, 2, 0, 6.2832); g.fill();
  g.fillStyle = '#8b97a0'; g.font = `500 10px ${FONT}`; g.textAlign = 'right';
  g.fillText('열전대', VW - 8, 282);
}

// particles
const sparks = [], smoke = [], gas = [];
let sparkAcc = 0, smokeAcc = 0, gasAcc = 0;
function spawnSpark(x, y, vx, vy, life, child) {
  sparks.push({ x, y, vx: vx ?? (Math.random() - 0.5) * 150, vy: vy ?? -(90 + Math.random() * 240), age: 0, life: life ?? (0.5 + Math.random() * 0.9), child: !!child, split: false });
}
function updateFX(rdt, sparkRate, burnI) {
  gasAcc += rdt * p.xO2 * 30;
  while (gasAcc >= 1) { gasAcc -= 1; gas.push({ x: 0, y: PIPE_Y + (Math.random() - 0.5) * 5, vx: 80 + Math.random() * 30, vy: 0, life: 0 }); }
  for (const q of gas) {
    q.life += rdt;
    if (q.x > GX0) { q.vx *= Math.pow(0.12, rdt); q.vy += (Math.random() - 0.65) * 50 * rdt; }
    q.x += q.vx * rdt; q.y += q.vy * rdt;
  }
  for (let i = gas.length - 1; i >= 0; i--) { const q = gas[i]; if (q.life > 3.2 || q.x > GX1 - 2 || q.y < GT) gas.splice(i, 1); }

  if (status === 'running') { sparkAcc += rdt * sparkRate; smokeAcc += rdt * burnI * 16; }
  let guard = 0;
  while (sparkAcc >= 1 && guard++ < 40) {
    sparkAcc -= 1;
    if (sparks.length > 320) continue;
    const s = sampleShape();
    const x = s.x0 + 6 + Math.random() * Math.max(4, s.w - 12);
    spawnSpark(x, topAt(x) - 1);
  }
  if (sparkAcc > 5) sparkAcc = 0;
  const n = sparks.length;
  for (let i = 0; i < n; i++) {
    const s = sparks[i];
    s.age += rdt; s.vy += 290 * rdt; s.x += s.vx * rdt; s.y += s.vy * rdt;
    if (s.y > GT) {
      if (s.x < GX0 + 1) { s.x = GX0 + 1; s.vx = -s.vx * 0.4; }
      if (s.x > GX1 - 1) { s.x = GX1 - 1; s.vx = -s.vx * 0.4; }
    }
    if (!s.child && !s.split && s.age > 0.55 * s.life && Math.random() < 0.4) {
      s.split = true;
      for (let k = 0; k < 3; k++) spawnSpark(s.x, s.y, s.vx * 0.3 + (Math.random() - 0.5) * 170, s.vy * 0.3 + (Math.random() - 0.5) * 170, 0.16 + Math.random() * 0.2, true);
    }
  }
  for (let i = sparks.length - 1; i >= 0; i--) if (sparks[i].age > sparks[i].life || sparks[i].y > BY) sparks.splice(i, 1);

  let sg = 0;
  while (smokeAcc >= 1 && sg++ < 20) {
    smokeAcc -= 1;
    const s = sampleShape();
    const x = s.x0 + Math.random() * s.w;
    smoke.push({ x, y: topAt(x) - 4, vx: (Math.random() - 0.5) * 10, vy: -(26 + Math.random() * 26), r: 3 + Math.random() * 3, age: 0, life: 2.6 + Math.random() });
  }
  for (const q of smoke) { q.age += rdt; q.x += q.vx * rdt; q.y += q.vy * rdt; q.r += 5 * rdt; if (q.y < GT + 10) q.vx += (q.x < CX ? -1 : 1) * 30 * rdt; }
  for (let i = smoke.length - 1; i >= 0; i--) if (smoke[i].age > smoke[i].life || smoke[i].y < -30) smoke.splice(i, 1);
}

function drawStage(t, vis) {
  g.setTransform(dpr * sScale, 0, 0, dpr * sScale, 0, 0);
  const bg = g.createLinearGradient(0, 0, 0, VH);
  bg.addColorStop(0, '#0a0f14'); bg.addColorStop(1, '#111921');
  g.fillStyle = bg; g.fillRect(0, 0, VW, VH);
  g.fillStyle = 'rgba(255,255,255,0.035)';
  for (let x = 10; x < VW; x += 20) for (let y = 10; y < VH; y += 20) g.fillRect(x, y, 1, 1);

  const metal = nearestMetal(p.k);
  drawBaseAndSink(t, metal, vis.lo.cond, vis.overheated);

  g.fillStyle = `rgba(90,160,255,${0.025 + 0.1 * p.xO2})`;
  g.fillRect(GX0, GT, GX1 - GX0, BY - GT);

  const topC = topAt(CX);
  drawSample(t, vis.col, vis.gl, vis.burnI);

  if (vis.gl || vis.burnI > 0) {
    const gc = vis.gl ? vis.gl.c : [255, 160, 60], ga = vis.gl ? vis.gl.a : 0.4;
    g.save(); g.globalCompositeOperation = 'lighter';
    const Rr = 30 + 120 * ga + 60 * vis.burnI;
    const rg = g.createRadialGradient(CX, topC + 6, 2, CX, topC + 6, Rr);
    rg.addColorStop(0, rgba(gc, 0.38 * ga + 0.25 * vis.burnI)); rg.addColorStop(1, rgba(gc, 0));
    g.fillStyle = rg; g.fillRect(CX - Rr, topC + 6 - Rr, 2 * Rr, 2 * Rr);
    g.restore();
  }

  for (const q of smoke) {
    g.fillStyle = `rgba(150,105,82,${0.24 * (1 - q.age / q.life)})`;
    g.beginPath(); g.arc(q.x, q.y, q.r, 0, 6.2832); g.fill();
  }

  const flameOn = status === 'running' && p.torchOn && p.P > 0;
  drawTorch(flameOn || (status === 'paused' && p.torchOn && p.P > 0));
  if (flameOn || (status === 'paused' && p.torchOn && p.P > 0)) drawFlame(t, topC);
  drawGlass();
  drawPipe();
  g.fillStyle = 'rgba(160,220,255,0.75)';
  for (const q of gas) { g.globalAlpha = q.x > GX0 ? Math.max(0, 1 - q.life / 3.2) : 0.85; g.beginPath(); g.arc(q.x, q.y, 1.4, 0, 6.2832); g.fill(); }
  g.globalAlpha = 1;
  drawThermocouple();
  drawBaseLabel(metal, vis.overheated);

  g.save(); g.globalCompositeOperation = 'lighter'; g.lineCap = 'round';
  for (const s of sparks) {
    const a = 1 - s.age / s.life;
    g.strokeStyle = `rgba(255,${(200 + 55 * a) | 0},${(110 + 110 * a) | 0},${a})`;
    g.lineWidth = s.child ? 0.8 : 1.4;
    g.beginPath(); g.moveTo(s.x, s.y); g.lineTo(s.x - s.vx * 0.035, s.y - s.vy * 0.035); g.stroke();
  }
  g.restore();
}

// ---------- charts ----------
const cT = $('cT'), cS = $('cS');
function sizeChart(c) {
  const w = c.clientWidth, h = c.clientHeight;
  if (c._w !== w || c._h !== h || c._dpr !== dpr) { c._w = w; c._h = h; c._dpr = dpr; c.width = Math.round(w * dpr); c.height = Math.round(h * dpr); }
}
let theme = {};
function readTheme() {
  const cs = getComputedStyle(document.documentElement), v = n => cs.getPropertyValue(n).trim();
  theme = { ink: v('--ink'), muted: v('--muted'), grid: v('--grid'), heat: v('--heat'), cool: v('--cool'), oxide: v('--oxide'), surface: v('--surface') };
}
const pick = (v, arr) => arr.find(a => a >= v) ?? arr[arr.length - 1];

function drawTChart() {
  sizeChart(cT);
  const ctx = cT.getContext('2d'), W = cT._w, H = cT._h;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  const L = 44, Rr = 40, Tt = 10, B = 24, pw = W - L - Rr, ph = H - Tt - B;
  const tMax = pick(Math.max(20, st.time * 1.1), [20, 30, 60, 90, 120, 180, 300, 600, 900, 1200, 1800, 2700, 3600, 5400]);
  const yMax = pick(Math.max(300, C(st.Tmax) + 60), [300, 400, 600, 800, 1000, 1200, 1500, 1800, 2000, 2500]);
  const xs = t => L + pw * t / tMax, ys = v => Tt + ph * (1 - v / yMax), yx = x => Tt + ph * (1 - x);
  ctx.font = `500 11px ${FONT}`; ctx.textBaseline = 'middle'; ctx.lineWidth = 1;
  const yStep = yMax <= 600 ? 100 : yMax <= 1200 ? 200 : 500;
  for (let v = 0; v <= yMax + 1e-6; v += yStep) {
    const y = Math.round(ys(v)) + 0.5;
    ctx.strokeStyle = theme.grid; ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke();
    ctx.fillStyle = theme.muted; ctx.textAlign = 'right'; ctx.fillText(String(v), L - 6, y);
  }
  ctx.textAlign = 'left';
  for (const f of [0, 0.5, 1]) ctx.fillText(`${f * 100}%`, L + pw + 6, yx(f));
  const xStep = pick(tMax / 5, [5, 10, 15, 30, 60, 120, 300, 600, 900, 1200]);
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (let v = 0; v <= tMax + 1e-6; v += xStep) ctx.fillText(tMax <= 180 ? `${v}s` : `${v / 60}분`, xs(v), Tt + ph + 6);
  ctx.textBaseline = 'middle';
  ctx.setLineDash([3, 4]); ctx.strokeStyle = theme.muted; ctx.globalAlpha = 0.6;
  for (const [v, lab] of [[570, '570 FeO 생성'], [1377, '1377 FeO 융해'], [1538, '1538 철 녹는점']]) {
    if (v >= yMax) continue;
    const y = Math.round(ys(v)) + 0.5;
    ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke();
    ctx.fillStyle = theme.muted; ctx.textAlign = 'left'; ctx.fillText(lab, L + 6, y - 7);
  }
  ctx.setLineDash([]); ctx.globalAlpha = 1;
  if (hist.length > 1) {
    ctx.strokeStyle = theme.oxide; ctx.lineWidth = 1.8; ctx.setLineDash([5, 4]);
    ctx.beginPath(); hist.forEach((h, i) => { const x = xs(h.t), y = yx(h.X); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); ctx.stroke();
    ctx.setLineDash([]);
    ctx.strokeStyle = theme.heat; ctx.lineWidth = 2.2; ctx.lineJoin = 'round';
    ctx.beginPath(); hist.forEach((h, i) => { const x = xs(h.t), y = ys(Math.min(h.T, yMax)); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); ctx.stroke();
  }
  const mark = (t, Tc, lab) => {
    const x = xs(t), y = ys(Math.min(Tc, yMax));
    ctx.fillStyle = theme.ink; ctx.beginPath(); ctx.arc(x, y, 4, 0, 6.2832); ctx.fill();
    ctx.font = `700 11px ${FONT}`; ctx.textAlign = x > L + pw - 50 ? 'right' : 'left';
    ctx.fillText(lab, x + (ctx.textAlign === 'left' ? 7 : -7), y - 9);
  };
  if (st.ignited) mark(st.ignTime, C(st.ignT), '점화');
  if (st.melted) { const h = hist.find(q => q.t >= st.meltTime); if (h) mark(h.t, h.T, '용융'); }
}

function drawSChart() {
  sizeChart(cS);
  const ctx = cS.getContext('2d'), W = cS._w, H = cS._h;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  const L = 46, Rr = 12, Tt = 10, B = 24, pw = W - L - Rr, ph = H - Tt - B;
  const TMIN = 25, TMAX = 2500;
  const pp = status === 'idle' ? { ...p, torchOn: true } : p;
  const xs = Tc => L + pw * (Tc - TMIN) / (TMAX - TMIN);
  const ys = w => Tt + ph * (1 - clamp(Math.log10(Math.max(1, w)), 0, 5) / 5);
  ctx.font = `500 11px ${FONT}`; ctx.lineWidth = 1;
  ['1 W', '10 W', '100 W', '1 kW', '10 kW', '100 kW'].forEach((lab, i) => {
    const y = Math.round(ys(10 ** i)) + 0.5;
    ctx.strokeStyle = theme.grid; ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke();
    ctx.fillStyle = theme.muted; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(lab, L - 6, y);
  });
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (let v = 500; v <= TMAX; v += 500) { ctx.textAlign = v === TMAX ? 'right' : 'center'; ctx.fillText(`${v}°C`, v === TMAX ? L + pw : xs(v), Tt + ph + 6); }
  ctx.setLineDash([3, 4]); ctx.strokeStyle = theme.muted; ctx.globalAlpha = 0.5;
  for (const v of [570, 1377, 1538]) { const x = Math.round(xs(v)) + 0.5; ctx.beginPath(); ctx.moveTo(x, Tt); ctx.lineTo(x, Tt + ph); ctx.stroke(); }
  ctx.globalAlpha = 1; ctx.setLineDash([]);

  const N = 500, Ts = [], gen = [], loss = [], torch = [];
  for (let i = 0; i <= N; i++) {
    const Tc = TMIN + (TMAX - TMIN) * i / N, Tk = Tc + 273.15;
    const tq = PHYS.torchIn(pp, Tk);
    Ts.push(Tc); torch.push(tq);
    gen.push(tq + ETA * PHYS.reaction(st, pp, Tk, PHYS.sQuasi(st, pp, Tk, 30)).Q);
    loss.push(PHYS.losses(st, pp, Tk).total);
  }
  const line = (arr, color, w, dash, alpha) => {
    ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = w; ctx.setLineDash(dash || []); ctx.globalAlpha = alpha ?? 1; ctx.lineJoin = 'round';
    ctx.beginPath(); arr.forEach((v, i) => { const x = xs(Ts[i]), y = ys(v); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); ctx.stroke(); ctx.restore();
  };
  if (pp.torchOn && pp.P > 0) line(torch, theme.heat, 1.5, [4, 4], 0.55);
  line(loss, theme.cool, 2.2);
  line(gen, theme.heat, 2.4);
  // crossings
  for (let i = 0; i < N; i++) {
    const f0 = gen[i] - loss[i], f1 = gen[i + 1] - loss[i + 1];
    if ((f0 > 0) === (f1 > 0) || Math.max(gen[i], loss[i]) < 1) continue;
    const u = f0 / (f0 - f1), Tc = Ts[i] + (Ts[i + 1] - Ts[i]) * u, w = loss[i] + (loss[i + 1] - loss[i]) * u;
    const stable = f0 > 0;
    ctx.beginPath(); ctx.arc(xs(Tc), ys(w), 5, 0, 6.2832);
    ctx.fillStyle = stable ? theme.ink : theme.surface; ctx.fill();
    ctx.strokeStyle = theme.ink; ctx.lineWidth = 1.6; ctx.stroke();
    ctx.fillStyle = theme.ink; ctx.font = `600 11px ${FONT}`; ctx.textAlign = 'center';
    ctx.textBaseline = stable ? 'bottom' : 'top';
    ctx.fillText(`${Math.round(Tc)}°C`, clamp(xs(Tc), L + 18, L + pw - 18), ys(w) + (stable ? -8 : 8));
  }
  // current state
  if (status !== 'idle') {
    const Tc = C(st.T), x = xs(clamp(Tc, TMIN, TMAX));
    ctx.save(); ctx.setLineDash([2, 3]); ctx.strokeStyle = theme.ink; ctx.globalAlpha = 0.5;
    ctx.beginPath(); ctx.moveTo(x, Tt); ctx.lineTo(x, Tt + ph); ctx.stroke(); ctx.restore();
    const gNow = PHYS.torchIn(pp, st.T) + ETA * PHYS.reaction(st, pp, st.T, st.s).Q;
    ctx.fillStyle = theme.heat; ctx.beginPath(); ctx.arc(x, ys(gNow), 3.5, 0, 6.2832); ctx.fill();
    ctx.fillStyle = theme.cool; ctx.beginPath(); ctx.arc(x, ys(PHYS.losses(st, pp, st.T).total), 3.5, 0, 6.2832); ctx.fill();
    ctx.fillStyle = theme.ink; ctx.font = `700 11px ${FONT}`; ctx.textAlign = x > L + pw - 60 ? 'right' : 'left'; ctx.textBaseline = 'top';
    ctx.fillText('현재', x + (ctx.textAlign === 'left' ? 5 : -5), Tt + 2);
  }
}

// ---------- HUD & verdict ----------
const VERD = {
  idle: ['⏻', '대기 중', () => '조건을 고른 뒤 스위치를 켜세요. 아래 열수지 도표의 ●이 지금 조건에서 온도가 멈출 지점입니다.'],
  none: ['○', '아직 변화 없음', () => `최고 ${Math.round(C(st.Tmax))} °C. 산화가 눈에 띌 만큼 진행되지 않았습니다. 토치를 올리거나 받침의 열전도도를 낮춰 보세요.`],
  noox: ['∅', '가열만 됨 (산소 없음)', () => '산소가 없어 반응열이 0입니다. 온도는 토치 열과 방열의 균형으로만 정해지고, 1538 °C를 넘으면 타지 않고 그냥 녹습니다.'],
  oxide: ['🟤', '녹 (느린 산화)', () => `발열이 방열을 넘지 못해 온도가 ${Math.round(C(st.T))} °C 근처에서 안정됐고, 표면에만 약 ${fmtLen(2.1 * st.s * 1e9)} 두께의 산화막이 자랐습니다. 관 안이 건조해서 붉은 수화 녹 대신 두께에 따라 밀짚색, 보라, 파랑, 회흑색으로 변하는 산화피막으로 나타납니다.`],
  forced: ['♨️', '빠른 산화 (토치에 기대는 중)', () => `산화가 빠르게 진행 중이지만 반응열이 방열을 한 번도 넘지 못했습니다. 받침이 열을 계속 빼앗고 있어서, 토치를 끄면 곧 식으며 멈춥니다.`],
  burn: ['🔥', '연소 (열폭주)', () => `${Math.round(C(st.ignT))} °C에서 반응열이 방열을 넘어섰습니다. 그 뒤로는 토치를 꺼도 스스로 온도를 올리며 탑니다. 연료가 줄어들면 다시 약해집니다.`],
  melt: ['🫠', '용융', () => '반응열은 끝내 방열을 넘지 못했고, 토치 열만으로 녹는점 1538 °C에 도달해 녹았습니다. 표면적이 작으면 타기 전에 먼저 녹습니다.'],
  burnmelt: ['🔥', '연소 후 용융', () => `${Math.round(C(st.ignT))} °C에서 열폭주가 시작됐고, 연소열로 온도가 치솟아 녹아내렸습니다.`],
};
function classify() {
  if (status === 'idle') return 'idle';
  if (st.ignited && st.melted) return 'burnmelt';
  if (st.ignited) return 'burn';
  if (st.melted) return 'melt';
  const dox = 2.1 * st.s * 1e9;
  if (st.X > 0.1) return 'forced';
  if (p.xO2 <= 0 && dox < 15) return 'noox';
  if (dox >= 15) return 'oxide';
  return 'none';
}
let lastVerdict = '';
function updateVerdict() {
  const k = classify();
  const [e, title, desc] = VERD[k];
  const d = desc();
  const key = k + d;
  if (key !== lastVerdict) {
    lastVerdict = key;
    $('verdict').dataset.k = k;
    $('vE').textContent = e; $('vT').textContent = title; $('vD').textContent = d;
    $('vK').textContent = status === 'running' ? '현재 판정' : status === 'done' ? '최종 판정' : '판정';
  }
  $('sTmax').textContent = status === 'idle' ? '—' : `${Math.round(C(st.Tmax))} °C`;
  $('sIgn').textContent = st.ignited ? `${Math.round(C(st.ignT))} °C` : '—';
  $('sX').textContent = status === 'idle' ? '—' : st.X < 0.001 ? (st.X > 0 ? '0.1% 미만' : '0%') : `${(st.X * 100).toFixed(1)}%`;
  $('sMelt').textContent = !st.melted ? '—' : st.fL >= 1 ? '완전 용융' : st.fL > 0 ? `${Math.round(st.fL * 100)}%` : (st.X > 0.9 ? '산화물 용융' : '응고됨');
}

function visState() {
  const rx = PHYS.reaction(st, p, st.T, st.s);
  const lo = PHYS.losses(st, p, st.T);
  const gen = ETA * rx.Q;
  const burning = rx.Q > 400 && (gen > lo.total * 0.5 || rx.limited);
  const burnI = burning ? Math.min(1, rx.Q / 12000 + 0.15) : 0;
  const sparkRate = rx.Q > 500 ? Math.min(260, rx.Q / 90) : 0;
  const metal = nearestMetal(p.k);
  const overheated = status !== 'idle' && C(PHYS.baseTopT(p.k, lo.cond)) > metal.mp;
  return { rx, lo, gen, burning, burnI, sparkRate, overheated, col: sampleColor(), gl: glowAt(C(st.T)) };
}

function updateHUD(vis) {
  $('hT').innerHTML = `${Math.round(C(st.T))}<small>°C</small>`;
  const dox = 2.1 * st.s * 1e9;
  $('hX').textContent = st.X < 0.001 ? `산화막 ${fmtLen(dox)}` : `산화율 ${(st.X * 100).toFixed(1)}%`;
  $('hQ').textContent = `반응열 ${fmtW(vis.gen)}, 방열 ${fmtW(vis.lo.total)}`;
  const pill = $('pill');
  let txt = '대기', cls = '';
  if (status === 'paused') txt = '일시정지';
  else if (status === 'done') txt = '종료';
  else if (status === 'running') {
    if (vis.burning && st.ignited) { txt = '연소 중'; cls = 'hot'; }
    else if (st.fL > 0 && st.fL < 1) { txt = '녹는 중'; cls = 'melt'; }
    else if (st.fL >= 1) { txt = '액체 상태'; cls = 'melt'; }
    else if (vis.burning) { txt = '격렬한 산화'; cls = 'hot'; }
    else if (st.dXdt > 2e-6 || (st.dTdt < 0.5 && 2.1 * st.s * 1e9 > 15 && vis.gen > 0.5)) { txt = '산화 진행'; cls = 'ox'; }
    else if (p.torchOn && p.P > 0) txt = '가열 중';
    else txt = '식는 중';
  }
  pill.textContent = txt; pill.className = 'pill ' + cls;
  $('clock').innerHTML = `${fmtTime(st.time)}<br>${status === 'running' ? (effSpeed < speed * 0.8 ? `슬로모션 ×${effSpeed < 1 ? effSpeed.toFixed(1) : Math.round(effSpeed)}` : `×${speed}`) : ''}`;
  $('eff').textContent = status === 'running' && effSpeed < speed * 0.8 ? '급변 구간이라 자동으로 느려졌어요' : '';
}

// ---------- loop ----------
let lastTs = performance.now(), frameN = 0;
function frame(ts) {
  const rdt = Math.min(0.05, Math.max(0, (ts - lastTs) / 1000)); lastTs = ts;
  if (status === 'running' && rdt > 0) {
    const adv = PHYS.step(st, p, rdt * speed, { dT: 250 * rdt, dX: 0.2 * rdt, dfL: 0.5 * rdt });
    effSpeed = effSpeed * 0.8 + (adv / rdt) * 0.2;
    pushHist(false);
    if (st.time >= 3600) { status = 'done'; syncButtons(); }
  }
  const vis = visState();
  updateFX(rdt, vis.sparkRate, vis.burnI);
  drawStage(ts / 1000, vis);
  frameN++;
  if (frameN % 3 === 0) { drawTChart(); drawSChart(); }
  if (frameN % 4 === 0) { updateHUD(vis); updateVerdict(); }
  requestAnimationFrame(frame);
}

function onResize() { resizeStage(); cT._w = 0; cS._w = 0; drawTChart(); drawSChart(); }
window.addEventListener('resize', onResize);
const mq = window.matchMedia('(prefers-color-scheme: dark)');
mq.addEventListener && mq.addEventListener('change', readTheme);
new MutationObserver(readTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

buildChips();
readTheme();
rebuildSample();
setControlsFromP();
syncButtons();
resizeStage();
updateVerdict();
requestAnimationFrame(frame);
})();
