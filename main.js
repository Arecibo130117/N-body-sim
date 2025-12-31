// main.js
// 전체 orchestrator: bodies(pool+typed arrays), integrator, collisions, fracture, thermal, render, UI
import * as THREE from "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js";

import { computeAccelerations, integrateVerlet, estimateSystemInvariants } from "./gravityIntegrator.js";
import { broadphaseSpatialHash, narrowphaseSpheres, computeImpactScalars, collisionImpulseResponse } from "./collision.js";
import { fractureAndSpawnFragments, tryReaccumulate } from "./fracture.js";
import { updateThermalState, applyImpactHeating, spawnFireSmokeParticles } from "./thermal.js";
import { createRenderSystem } from "./render.js";
import { createUI } from "./ui.js";

/**
 * Material presets (요구: 최소 3개 이상 + 룩/강도/열반응/발광 색감 차등)
 */
function makeMaterials() {
  const v = (x,y,z)=> new THREE.Vector3(x,y,z);
  return [
    {
      name: "rock",
      rho: 2600,
      restitution: 0.35,
      friction: 0.55,
      Qstar: 8.0,            // 파괴 임계(특정에너지)
      heatCapacity: 800,     // J/kgK
      heatAbsorb: 0.9,
      meltTemp: 1400,
      vaporTemp: 3200,
      plasmaTemp: 9000,
      glowTempStart: 1100,
      emissiveScale: 1.0,
      cooling: 0.18,
      fragmentAlpha: 1.9,
      baseColor: v(0.55, 0.52, 0.50),
      glowColor: v(1.00, 0.55, 0.18),
      trailTint: v(0.85, 0.80, 0.72),
    },
    {
      name: "ice",
      rho: 920,
      restitution: 0.25,
      friction: 0.18,
      Qstar: 3.0,
      heatCapacity: 2100,
      heatAbsorb: 0.7,
      meltTemp: 273,
      vaporTemp: 800,
      plasmaTemp: 6000,
      glowTempStart: 500,
      emissiveScale: 0.85,
      cooling: 0.35,
      fragmentAlpha: 2.2,
      baseColor: v(0.55, 0.70, 0.85),
      glowColor: v(0.65, 0.90, 1.00),
      trailTint: v(0.60, 0.82, 1.00),
    },
    {
      name: "metal",
      rho: 7800,
      restitution: 0.45,
      friction: 0.45,
      Qstar: 14.0,
      heatCapacity: 500,
      heatAbsorb: 1.0,
      meltTemp: 1800,
      vaporTemp: 5000,
      plasmaTemp: 10000,
      glowTempStart: 1200,
      emissiveScale: 1.25,
      cooling: 0.12,
      fragmentAlpha: 1.7,
      baseColor: v(0.70, 0.70, 0.75),
      glowColor: v(1.00, 0.85, 0.35),
      trailTint: v(0.95, 0.95, 1.00),
    },
    {
      name: "gas",
      rho: 0.9,
      restitution: 0.05,
      friction: 0.02,
      Qstar: 0.8,
      heatCapacity: 1000,
      heatAbsorb: 0.6,
      meltTemp: 120,
      vaporTemp: 250,
      plasmaTemp: 4000,
      glowTempStart: 350,
      emissiveScale: 0.95,
      cooling: 0.55,
      fragmentAlpha: 2.4,
      baseColor: v(0.45, 0.55, 0.70),
      glowColor: v(0.70, 0.85, 1.00),
      trailTint: v(0.55, 0.70, 0.95),
    }
  ];
}

/**
 * Bodies pool with TypedArrays (성능/안정성)
 */
function createBodies(capacity, materials) {
  const pos = new Float32Array(capacity * 3);
  const vel = new Float32Array(capacity * 3);
  const acc = new Float32Array(capacity * 3);
  const accNew = new Float32Array(capacity * 3);
  const omega = new Float32Array(capacity * 3);

  const mass = new Float32Array(capacity);
  const radius = new Float32Array(capacity);
  const damage = new Float32Array(capacity);

  // Thermal
  const heat = new Float32Array(capacity);       // Joule
  const temperature = new Float32Array(capacity);// Kelvin
  const meltFrac = new Float32Array(capacity);
  const vaporFrac = new Float32Array(capacity);
  const emissive = new Float32Array(capacity);

  const materialId = new Uint8Array(capacity);
  const active = new Uint8Array(capacity);

  const free = [];
  for (let i = capacity - 1; i >= 0; i--) free.push(i);

  function alloc() {
    if (free.length === 0) return -1;
    const id = free.pop();
    active[id] = 1;
    damage[id] = 0;
    omega[id * 3] = omega[id * 3 + 1] = omega[id * 3 + 2] = 0;
    heat[id] = 0;
    temperature[id] = 300;
    meltFrac[id] = 0;
    vaporFrac[id] = 0;
    emissive[id] = 0;
    return id;
  }

  function release(id) {
    active[id] = 0;
    free.push(id);
  }

  function clearAll() {
    for (let i = 0; i < capacity; i++) active[i] = 0;
    free.length = 0;
    for (let i = capacity - 1; i >= 0; i--) free.push(i);
  }

  function addBody({ p, v, m, r, matName, tempK = 300 }) {
    const id = alloc();
    if (id < 0) return -1;

    pos[id * 3] = p.x; pos[id * 3 + 1] = p.y; pos[id * 3 + 2] = p.z;
    vel[id * 3] = v.x; vel[id * 3 + 1] = v.y; vel[id * 3 + 2] = v.z;

    mass[id] = m;
    radius[id] = r;

    const idx = materials.findIndex(x => x.name === matName);
    materialId[id] = idx >= 0 ? idx : 0;

    temperature[id] = tempK;
    heat[id] = (tempK - 300) * m * materials[materialId[id]].heatCapacity;

    return id;
  }

  return {
    capacity,
    pos, vel, acc, accNew, omega,
    mass, radius, damage,
    heat, temperature, meltFrac, vaporFrac, emissive,
    materialId, active,
    alloc, release, clearAll, addBody,
    freeCount: () => free.length
  };
}

// --- Parameters / Defaults
const params = {
  // Gravity & integrator
  G: 1.0,
  eps: 0.6,
  dt: 0.012,
  substeps: 2,

  enableCollisions: true,
  enableFracture: true,
  enableThermal: true,

  collisionPositionCorrection: 0.35,

  // Visual
  exposure: 1.05,
  bloomStrength: 1.05,
  bloomThreshold: 0.18,
  bloomRadius: 0.35,

  trailLength: 128,
  maxTrails: 128,

  starCount: 7000,
  starRadius: 2600,

  // Budget & degrade
  particleBudgetHard: 14000,
  particleBudgetSoft: 7000,
  maxFragmentsPerEvent: 56,
  maxSpallFragments: 36,
  autoDegrade: true,
  degradeTargetFPS: 35,

  // Damage/fracture tuning
  damageGain: 0.65,
  minFractureMass: 0.7,
  minFragmentMass: 0.06,
  coreKeepChance: 0.7,
  coreKeepFrac: 0.25,
  coreKeepMass: 60,

  // Ejecta parameters
  ejectaEnergyFracSpall: 0.22,
  ejectaEnergyFracFrag: 0.42,
  spallLayerFrac: 0.04,
  spallConeAngleRad: 0.65,
  spallSpeedBase: 6.0,
  spallSpeedJitter: 7.5,
  spallTempBoost: 600,

  fragmentConeAngleRad: 0.75,
  fragmentSpeedBase: 12.0,
  fragmentTempBoost: 1300,

  // Thermal
  ambientTemp: 300,
  heatEtaBase: 0.28,
  fireTempThreshold: 1200,

  // Particles spawn rates
  flashParticlesBase: 120,
  fireParticlesBase: 140,
  smokeParticlesBase: 110,

  // Re-accretion
  reaccumulationEnabled: true,
  reaccumulationDistanceFactor: 0.92,
  reaccumulationMaxRelSpeed: 1.2,
  reaccumulationMaxTemp: 650,

  // Spawn UI
  addBodyMode: false,
  spawnMass: 40,
  spawnRadius: 6,
  spawnTemp: 340,
  spawnSpeedScale: 12,
  spawnMaterialPreset: "rock",

  // Wind for particles
  wind: new THREE.Vector3(0.4, 0.0, 0.2),
};

params.materials = makeMaterials();

// --- Init
const container = document.getElementById("app");
const bodies = createBodies(2200, params.materials);
const renderSys = createRenderSystem(container, params, bodies);

const ui = createUI(params, {
  reset: () => {
    bodies.clearAll();
    // trails cleanup
    for (const [idx] of renderSys.trailSys.trails.entries()) renderSys.trailSys.remove(renderSys.scene, idx);
  },
  spawnAtOrigin: () => spawnBodyAt(new THREE.Vector3(0,0,0), new THREE.Vector3(0,0,0)),
  loadPreset: (name) => loadPreset(name)
});

// --- Spawn (click+drag)
const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();
const spawnPlane = new THREE.Plane();
let dragStart = null; // { worldPos: Vector3 }
let dragNow = null;

function screenToWorldOnPlane(clientX, clientY, plane) {
  pointer.x = (clientX / window.innerWidth) * 2 - 1;
  pointer.y = -(clientY / window.innerHeight) * 2 + 1;
  raycaster.setFromCamera(pointer, renderSys.camera);
  const out = new THREE.Vector3();
  raycaster.ray.intersectPlane(plane, out);
  return out;
}

function updateSpawnPlane() {
  // 카메라 시선에 수직인 평면을 원점(또는 controls target) 부근에 둠
  const normal = new THREE.Vector3();
  renderSys.camera.getWorldDirection(normal);
  const point = renderSys.controls.target.clone();
  spawnPlane.setFromNormalAndCoplanarPoint(normal, point);
}

function spawnBodyAt(p, v) {
  const m = params.spawnMass;
  const r = params.spawnRadius;
  const matName = params.spawnMaterialPreset;
  const temp = params.spawnTemp;

  bodies.addBody({ p, v, m, r, matName, tempK: temp });

  // trail budget: 큰 바디만 자동으로 trail 부여(또는 사용자가 많이 넣으면 제한됨)
  const id = bodies.capacity - bodies.freeCount() - 1; // 근사(정확 id는 addBody return 받는게 정석)
  // 정확히는 addBody 반환값 사용
}

renderSys.renderer.domElement.addEventListener("pointerdown", (e) => {
  if (!params.addBodyMode) return;
  updateSpawnPlane();
  const w = screenToWorldOnPlane(e.clientX, e.clientY, spawnPlane);
  dragStart = { worldPos: w.clone() };
  dragNow = w.clone();

  renderSys.arrow.visible = true;
  renderSys.arrow.position.copy(dragStart.worldPos);
});

renderSys.renderer.domElement.addEventListener("pointermove", (e) => {
  if (!params.addBodyMode || !dragStart) return;
  updateSpawnPlane();
  dragNow = screenToWorldOnPlane(e.clientX, e.clientY, spawnPlane);

  const dir = dragNow.clone().sub(dragStart.worldPos);
  const len = dir.length();
  if (len > 1e-6) {
    dir.normalize();
    renderSys.arrow.setDirection(dir);
    renderSys.arrow.setLength(Math.min(120, len), 8, 4);
  }
});

window.addEventListener("pointerup", (e) => {
  if (!params.addBodyMode || !dragStart || !dragNow) return;

  const dir = dragNow.clone().sub(dragStart.worldPos);
  const len = dir.length();
  const v = len > 1e-6 ? dir.multiplyScalar(params.spawnSpeedScale) : new THREE.Vector3(0,0,0);

  const id = bodies.addBody({
    p: dragStart.worldPos.clone(),
    v,
    m: params.spawnMass,
    r: params.spawnRadius,
    matName: params.spawnMaterialPreset,
    tempK: params.spawnTemp
  });

  // trail 자동 부여(예산 내에서)
  if (id >= 0) renderSys.trailSys.ensure(renderSys.scene, id);

  dragStart = null; dragNow = null;
  renderSys.arrow.visible = false;
});

// --- LOD selection (간이): 큰 바디/활성 바디 일부만 trail 유지
function refreshTrailSelection() {
  // 기존 trail이 너무 많으면 큰 질량 우선 유지
  const activeIds = [];
  for (let i = 0; i < bodies.capacity; i++) if (bodies.active[i]) activeIds.push(i);

  activeIds.sort((a, b) => bodies.mass[b] - bodies.mass[a]);
  const keep = new Set(activeIds.slice(0, params.maxTrails));

  for (const [idx] of renderSys.trailSys.trails.entries()) {
    if (!keep.has(idx)) renderSys.trailSys.remove(renderSys.scene, idx);
  }
  for (let k = 0; k < Math.min(params.maxTrails, activeIds.length); k++) {
    renderSys.trailSys.ensure(renderSys.scene, activeIds[k]);
  }
}

// --- Presets
function loadPreset(name) {
  bodies.clearAll();
  for (const [idx] of renderSys.trailSys.trails.entries()) renderSys.trailSys.remove(renderSys.scene, idx);

  if (name === "galaxy") presetGalaxyDisk();
  if (name === "impact") presetImpact();
  if (name === "reaccrete") presetReaccretion();

  refreshTrailSelection();
}

function presetGalaxyDisk() {
  // disk around origin
  const N = 520;
  const centerMass = 9000;
  const centerR = 18;

  const c = bodies.addBody({
    p: new THREE.Vector3(0,0,0),
    v: new THREE.Vector3(0,0,0),
    m: centerMass,
    r: centerR,
    matName: "metal",
    tempK: 450
  });
  if (c >= 0) renderSys.trailSys.ensure(renderSys.scene, c);

  for (let i = 0; i < N; i++) {
    const a = 120 + 520 * Math.pow(Math.random(), 1.25);
    const ang = Math.random() * Math.PI * 2;
    const h = (Math.random() - 0.5) * 22;

    const x = Math.cos(ang) * a;
    const z = Math.sin(ang) * a;
    const y = h;

    // circular velocity v = sqrt(GM/r)
    const vmag = Math.sqrt(params.G * centerMass / (a + 1e-6)) * (0.85 + 0.3 * Math.random());
    const vx = -Math.sin(ang) * vmag;
    const vz =  Math.cos(ang) * vmag;
    const vy = (Math.random() - 0.5) * 0.05 * vmag;

    const m = 0.6 + 6.0 * Math.pow(Math.random(), 2.0);
    const r = 0.7 + 2.2 * Math.pow(m, 0.33);

    const mat = (Math.random() < 0.12) ? "ice" : "rock";
    const id = bodies.addBody({ p: new THREE.Vector3(x,y,z), v: new THREE.Vector3(vx,vy,vz), m, r, matName: mat, tempK: 280 + 60 * Math.random() });
    if (id >= 0 && Math.random() < 0.15) renderSys.trailSys.ensure(renderSys.scene, id);
  }
}

function presetImpact() {
  // two big bodies high speed + some test debris
  const A = bodies.addBody({
    p: new THREE.Vector3(-90, 0, 0),
    v: new THREE.Vector3(26, 0.2, 0),
    m: 900,
    r: 14,
    matName: "rock",
    tempK: 420
  });
  const B = bodies.addBody({
    p: new THREE.Vector3(90, 0, 0),
    v: new THREE.Vector3(-20, -0.3, 0),
    m: 700,
    r: 12,
    matName: "metal",
    tempK: 380
  });
  if (A >= 0) renderSys.trailSys.ensure(renderSys.scene, A);
  if (B >= 0) renderSys.trailSys.ensure(renderSys.scene, B);

  // small target cloud
  for (let i = 0; i < 120; i++) {
    const p = new THREE.Vector3(
      (Math.random() - 0.5) * 40,
      (Math.random() - 0.5) * 20,
      (Math.random() - 0.5) * 40
    );
    const v = new THREE.Vector3((Math.random() - 0.5) * 2, (Math.random() - 0.5) * 2, (Math.random() - 0.5) * 2);
    const m = 0.3 + 2.0 * Math.random();
    const r = 0.6 + 1.6 * Math.pow(m, 0.33);
    bodies.addBody({ p, v, m, r, matName: (Math.random() < 0.5) ? "ice" : "rock", tempK: 260 });
  }
}

function presetReaccretion() {
  // cluster of fragments with low relative speeds, should gradually merge
  const centerMass = 2200;
  const core = bodies.addBody({
    p: new THREE.Vector3(0,0,0),
    v: new THREE.Vector3(0,0,0),
    m: centerMass,
    r: 16,
    matName: "rock",
    tempK: 380
  });
  if (core >= 0) renderSys.trailSys.ensure(renderSys.scene, core);

  for (let i = 0; i < 680; i++) {
    const r = 60 * Math.pow(Math.random(), 1.1);
    const theta = Math.random() * 2 * Math.PI;
    const phi = Math.acos(2 * Math.random() - 1);

    const p = new THREE.Vector3(
      r * Math.sin(phi) * Math.cos(theta),
      r * Math.cos(phi) * 0.55,
      r * Math.sin(phi) * Math.sin(theta)
    );
    const v = new THREE.Vector3(
      (Math.random() - 0.5) * 1.2,
      (Math.random() - 0.5) * 0.8,
      (Math.random() - 0.5) * 1.2
    );
    const m = 0.25 + 3.0 * Math.pow(Math.random(), 1.8);
    const rr = 0.55 + 1.9 * Math.pow(m, 0.33);
    const mat = (Math.random() < 0.25) ? "ice" : "rock";
    const id = bodies.addBody({ p, v, m, r: rr, matName: mat, tempK: 290 + 40 * Math.random() });
    if (id >= 0 && Math.random() < 0.1) renderSys.trailSys.ensure(renderSys.scene, id);
  }
}

// Default preset
loadPreset("galaxy");

// --- Fragment emitter used by fracture.js
function emitFragment({ pos, vel, mass, radius, materialId, temperatureBoost }) {
  const matName = params.materials[materialId]?.name ?? "rock";
  const id = bodies.addBody({ p: pos, v: vel, m: mass, r: radius, matName, tempK: params.ambientTemp + temperatureBoost });
  if (id >= 0 && Math.random() < 0.35) renderSys.trailSys.ensure(renderSys.scene, id);
}

// --- Merge bookkeeping (re-accretion)
const mergePairs = [];
function markMerge(i, j) { mergePairs.push([i, j]); }

function mergeBodies(i, j) {
  const { pos, vel, mass, radius, active, omega, materialId, heat, meltFrac, vaporFrac, damage } = bodies;
  if (!active[i] || !active[j]) return;

  // Combine mass and momentum
  const mi = mass[i], mj = mass[j];
  const m = mi + mj;

  const pi = new THREE.Vector3(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);
  const pj = new THREE.Vector3(pos[j * 3], pos[j * 3 + 1], pos[j * 3 + 2]);

  const vi = new THREE.Vector3(vel[i * 3], vel[i * 3 + 1], vel[i * 3 + 2]);
  const vj = new THREE.Vector3(vel[j * 3], vel[j * 3 + 1], vel[j * 3 + 2]);

  const P = vi.clone().multiplyScalar(mi).add(vj.clone().multiplyScalar(mj));
  const v = P.multiplyScalar(1 / m);

  // COM position
  const p = pi.clone().multiplyScalar(mi).add(pj.clone().multiplyScalar(mj)).multiplyScalar(1 / m);

  // angular momentum: L = Σ r × (m v) + Iω (here: just keep approximate spin energy)
  // We keep a rough omega by weighted average (simple)
  const oi = new THREE.Vector3(omega[i * 3], omega[i * 3 + 1], omega[i * 3 + 2]);
  const oj = new THREE.Vector3(omega[j * 3], omega[j * 3 + 1], omega[j * 3 + 2]);
  const o = oi.multiplyScalar(mi).add(oj.multiplyScalar(mj)).multiplyScalar(1 / m);

  // new radius: conserve density-ish (pick heavier body's material)
  const mid = (mi >= mj) ? materialId[i] : materialId[j];
  const rho = params.materials[mid].rho;
  const rNew = Math.cbrt(m / (rho * (4 / 3) * Math.PI));

  // write into i, release j
  pos[i * 3] = p.x; pos[i * 3 + 1] = p.y; pos[i * 3 + 2] = p.z;
  vel[i * 3] = v.x; vel[i * 3 + 1] = v.y; vel[i * 3 + 2] = v.z;
  omega[i * 3] = o.x; omega[i * 3 + 1] = o.y; omega[i * 3 + 2] = o.z;

  mass[i] = m;
  radius[i] = rNew;
  materialId[i] = mid;

  // thermal combine
  heat[i] += heat[j];
  meltFrac[i] = Math.min(1, (meltFrac[i] + meltFrac[j]) * 0.5);
  vaporFrac[i] = Math.min(1, (vaporFrac[i] + vaporFrac[j]) * 0.35);
  damage[i] = Math.min(1, 0.35 + 0.25 * Math.random()); // rubble pile: 찌그러짐을 남김

  active[j] = 0;
}

// --- Simulation loop
let lastT = performance.now();
let fpsEMA = 60;
let invMon = { E0: null, L0: null, E: 0, K: 0, U: 0, L: new THREE.Vector3() };
let invTimer = 0;

function autoDegradeIfNeeded(dt) {
  if (!params.autoDegrade) return;
  const fps = 1 / Math.max(1e-6, dt);
  fpsEMA = fpsEMA * 0.9 + fps * 0.1;

  if (fpsEMA < params.degradeTargetFPS) {
    // 자연스러운 디그레이드: trail↓, fragments↓, particles↓, substeps↓
    params.trailLength = Math.max(32, Math.floor(params.trailLength * 0.9));
    params.maxFragmentsPerEvent = Math.max(10, Math.floor(params.maxFragmentsPerEvent * 0.9));
    params.maxSpallFragments = Math.max(10, Math.floor(params.maxSpallFragments * 0.9));
    params.particleBudgetSoft = Math.max(800, Math.floor(params.particleBudgetSoft * 0.9));
    params.substeps = Math.max(1, params.substeps - 1);
  }
}

function stepPhysics(dt) {
  const n = bodies.capacity;

  // integrator substeps
  const h = dt / params.substeps;
  for (let s = 0; s < params.substeps; s++) {
    // compute accelerations at x(t)
    computeAccelerations(bodies.pos, bodies.mass, bodies.active, n, params.G, params.eps, bodies.acc);

    // integrate positions with velocity verlet needs accNew at x(t+dt)
    // 1) half kick+drift in integrateVerlet() after we compute accNew
    // We'll do: predict new positions by half-step method:
    // Here we do standard: v+=a*dt/2, x+=v*dt, then recompute aNew, then v+=aNew*dt/2
    // but we need accNew computed from new positions
    // We'll reuse bodies.accNew buffer:
    // - do half kick + drift first inside integrateVerlet using accNew computed after we set it.
    // Approach: inline: v += a*h/2, x += v*h, compute aNew, then call integrateVerlet to finish swap+half kick.
    // To keep code simple, we do:
    //   - v += a*h/2; x += v*h
    //   - compute accNew
    //   - set bodies.accNew and call integrateVerlet with dt=h but note integrateVerlet expects accNew already computed after drift,
    //     so we temporarily do a trick: integrateVerlet uses acc and accNew; we'll pre-drift manually and then do final half kick only.
    // Instead: do real velocity-verlet with integrateVerlet provided in gravityIntegrator.js:
    //   integrateVerlet does: v+=a*h/2; x+=v*h; (caller computes accNew); then it copies accNew->acc and does v+=a*h/2.
    // So we must compute accNew between its stages. We'll emulate by calling integrateVerlet in two halves is messy.
    // => Simpler: perform the same algorithm here:

    // v += a*h/2
    const half = 0.5 * h;
    for (let i = 0; i < n; i++) {
      if (!bodies.active[i]) continue;
      const k = i * 3;
      bodies.vel[k]     += bodies.acc[k] * half;
      bodies.vel[k + 1] += bodies.acc[k + 1] * half;
      bodies.vel[k + 2] += bodies.acc[k + 2] * half;
    }

    // x += v*h
    for (let i = 0; i < n; i++) {
      if (!bodies.active[i]) continue;
      const k = i * 3;
      bodies.pos[k]     += bodies.vel[k] * h;
      bodies.pos[k + 1] += bodies.vel[k + 1] * h;
      bodies.pos[k + 2] += bodies.vel[k + 2] * h;
    }

    // recompute aNew
    computeAccelerations(bodies.pos, bodies.mass, bodies.active, n, params.G, params.eps, bodies.accNew);

    // v += aNew*h/2, and acc <- accNew
    for (let i = 0; i < n; i++) {
      if (!bodies.active[i]) continue;
      const k = i * 3;
      bodies.vel[k]     += bodies.accNew[k] * half;
      bodies.vel[k + 1] += bodies.accNew[k + 1] * half;
      bodies.vel[k + 2] += bodies.accNew[k + 2] * half;

      bodies.acc[k] = bodies.accNew[k];
      bodies.acc[k + 1] = bodies.accNew[k + 1];
      bodies.acc[k + 2] = bodies.accNew[k + 2];
    }

    // collisions & fracture & thermal in each substep for stability
    if (params.enableCollisions) {
      // broadphase
      // cellSize: 평균 반지름 기반
      const cellSize = 2.5 * Math.max(1.0, estimateMaxRadius(120));
      const pairs = broadphaseSpatialHash(bodies.pos, bodies.radius, bodies.active, n, cellSize);
      const hits = narrowphaseSpheres(bodies.pos, bodies.radius, bodies.active, pairs);

      // impact processing
      for (let hIdx = 0; hIdx < hits.length; hIdx++) {
        const impact = computeImpactScalars(hits[hIdx], bodies.pos, bodies.vel, bodies.mass);

        // collision impulse response (momentum conservation baseline)
        const impulseInfo = collisionImpulseResponse(impact, bodies, params);

        if (params.enableThermal) {
          applyImpactHeating(impact, bodies, params);
          spawnFireSmokeParticles(impact, bodies, params, (req) => {
            renderSys.particleSys.spawnBurst(req, params.particleBudgetSoft);
          });
        }

        if (params.enableFracture) {
          fractureAndSpawnFragments(impact, bodies, params, emitFragment);
        }
      }

      // Re-accretion pass (cheap): reuse broadphase pairs
      mergePairs.length = 0;
      tryReaccumulate(bodies, params, pairs, markMerge);

      // merge (avoid chain merges by single pass)
      for (let mIdx = 0; mIdx < mergePairs.length; mIdx++) {
        const [a, b] = mergePairs[mIdx];
        if (bodies.active[a] && bodies.active[b]) mergeBodies(a, b);
      }
    }

    // thermal decay
    if (params.enableThermal) updateThermalState(bodies, params, h);
  }
}

function estimateMaxRadius(sample = 80) {
  // 샘플링으로 max radius 대략 추정
  let maxR = 1.0;
  const n = bodies.capacity;
  for (let s = 0; s < sample; s++) {
    const i = (Math.random() * n) | 0;
    if (!bodies.active[i]) continue;
    maxR = Math.max(maxR, bodies.radius[i]);
  }
  return maxR;
}

function countActive() {
  let c = 0;
  for (let i = 0; i < bodies.capacity; i++) if (bodies.active[i]) c++;
  return c;
}

function animate() {
  const now = performance.now();
  const dtReal = Math.min(0.05, (now - lastT) / 1000);
  lastT = now;

  autoDegradeIfNeeded(dtReal);

  // physics step with user dt (fixed-ish)
  const dtSim = params.dt;
  stepPhysics(dtSim);

  // refresh trail selection occasionally
  if ((now | 0) % 700 < 16) refreshTrailSelection();

  // invariants monitoring every ~0.5s
  invTimer += dtReal;
  if (invTimer > 0.5) {
    invTimer = 0;
    invMon = estimateSystemInvariants(bodies.pos, bodies.vel, bodies.mass, bodies.active, bodies.capacity, params.G, params.eps, 900);
    if (!invMon.E0) {
      invMon.E0 = invMon.E;
      invMon.L0 = invMon.L.clone();
    }
  }

  // overlay
  const N = countActive();
  const dE = invMon.E0 ? (invMon.E - invMon.E0) : 0;
  const dL = invMon.L0 ? invMon.L.clone().sub(invMon.L0).length() : 0;

  ui.setOverlay(
`FPS(EMA): ${fpsEMA.toFixed(1)}
Bodies: ${N} / cap ${bodies.capacity}
Particles: ${renderSys.particleSys.aliveCount} / soft ${params.particleBudgetSoft} / hard ${params.particleBudgetHard}

E ≈ ${(invMon.E).toExponential(3)}   (ΔE ≈ ${dE.toExponential(2)})
K ≈ ${(invMon.K).toExponential(3)}
U ≈ ${(invMon.U).toExponential(3)}
|L| ≈ ${invMon.L.length().toExponential(3)}   (Δ|L| ≈ ${dL.toExponential(2)})

dt=${params.dt.toFixed(4)}  substeps=${params.substeps}  G=${params.G.toFixed(3)}  ε=${params.eps.toFixed(2)}
Bloom: S=${params.bloomStrength.toFixed(2)} T=${params.bloomThreshold.toFixed(2)} R=${params.bloomRadius.toFixed(2)}
Trail: len=${params.trailLength}  max=${params.maxTrails}
`);
  // render
  renderSys.render(dtReal);

  requestAnimationFrame(animate);
}

requestAnimationFrame(animate);
