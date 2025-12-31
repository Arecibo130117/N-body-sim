// fracture.js
// damage 누적 → 임계치 초과 시 fragmentation
// 파편 질량분포: 파워법칙 근사 + 보존량(질량/운동량) 우선 배분
// spallation(표면 박리) + re-accretion(저속/근접/저온이면 병합)

import * as THREE from "three";

function randUnitVecBiased(coneAxis, coneAngleRad) {
  const u = Math.random();
  const v = Math.random();
  const cosA = Math.cos(coneAngleRad);
  const z = cosA + (1 - cosA) * u;
  const phi = 2 * Math.PI * v;
  const r = Math.sqrt(Math.max(0, 1 - z * z));
  const local = new THREE.Vector3(r * Math.cos(phi), r * Math.sin(phi), z);

  const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), coneAxis.clone().normalize());
  return local.applyQuaternion(q).normalize();
}

function powerLawMasses(mTotal, count, alpha = 1.8, minFrac = 0.002) {
  const masses = new Array(count);
  let sum = 0;
  for (let i = 0; i < count; i++) {
    const x = Math.random();
    const w = Math.pow(Math.max(1e-6, x), -1 / (alpha));
    masses[i] = w;
    sum += w;
  }
  const minM = mTotal * minFrac / count;
  let sum2 = 0;
  for (let i = 0; i < count; i++) {
    masses[i] = Math.max(minM, (masses[i] / sum) * mTotal);
    sum2 += masses[i];
  }
  const scale = mTotal / sum2;
  for (let i = 0; i < count; i++) masses[i] *= scale;
  return masses;
}

function approxRadiusFromMass(m, rho) {
  return Math.cbrt(m / (rho * (4 / 3) * Math.PI));
}

export function fractureAndSpawnFragments(impact, bodies, params, emit) {
  const { i, j, n, Eimp, Q, contactPoint, relVel } = impact;
  const { mass, radius, damage, active, materialId } = bodies;

  const matI = params.materials[materialId[i]];
  const matJ = params.materials[materialId[j]];

  const QsI = matI.Qstar;
  const QsJ = matJ.Qstar;

  damage[i] += (Q / (QsI + 1e-9)) * params.damageGain;
  damage[j] += (Q / (QsJ + 1e-9)) * params.damageGain;

  const doSpall = Q > 0.35 * Math.min(QsI, QsJ);

  const shouldFracI = damage[i] > 1 && mass[i] > params.minFractureMass && active[i];
  const shouldFracJ = damage[j] > 1 && mass[j] > params.minFractureMass && active[j];

  const impactNormal = n.clone();
  const impactSpeed = relVel.length();

  if (doSpall) {
    spallateBody(i, impactNormal.clone().negate(), contactPoint, impactSpeed, bodies, matI, params, emit);
    spallateBody(j, impactNormal.clone(), contactPoint, impactSpeed, bodies, matJ, params, emit);
  }

  if (shouldFracI) fragmentBody(i, impactNormal.clone().negate(), Eimp, bodies, matI, params, emit);
  if (shouldFracJ) fragmentBody(j, impactNormal.clone(), Eimp, bodies, matJ, params, emit);
}

function spallateBody(idx, outwardNormal, contactPoint, impactSpeed, bodies, mat, params, emit) {
  const { pos, vel, mass, radius, active, materialId, damage } = bodies;
  if (!active[idx]) return;

  const m0 = mass[idx];
  const r0 = radius[idx];

  const layerFrac = Math.min(0.08, params.spallLayerFrac * (0.25 + Math.min(1, impactSpeed / 40)));
  const mSpall = m0 * layerFrac;
  if (mSpall < params.minFragmentMass * 4) return;

  const count = Math.min(params.maxSpallFragments, Math.max(6, (mSpall / params.minFragmentMass) | 0));
  const masses = powerLawMasses(mSpall, count, 2.2, 0.001);

  const basePos = new THREE.Vector3(pos[idx * 3], pos[idx * 3 + 1], pos[idx * 3 + 2]);
  const baseVel = new THREE.Vector3(vel[idx * 3], vel[idx * 3 + 1], vel[idx * 3 + 2]);

  const Eej = params.ejectaEnergyFracSpall * 0.5 * mSpall * impactSpeed * impactSpeed;

  const vrels = [];
  let pSum = new THREE.Vector3();

  const cone = params.spallConeAngleRad;
  for (let k = 0; k < count; k++) {
    const dir = randUnitVecBiased(outwardNormal, cone);
    const speed = (params.spallSpeedBase + Math.random() * params.spallSpeedJitter) * (0.6 + 0.4 * Math.random());
    const vrel = dir.multiplyScalar(speed);
    vrels.push(vrel);
    pSum.add(vrel.clone().multiplyScalar(masses[k]));
  }

  const vCorr = pSum.multiplyScalar(1 / mSpall);
  for (let k = 0; k < count; k++) vrels[k].sub(vCorr);

  let curE = 0;
  for (let k = 0; k < count; k++) curE += 0.5 * masses[k] * vrels[k].lengthSq();
  const s = Math.sqrt(Math.max(1e-9, Eej / (curE + 1e-9)));
  for (let k = 0; k < count; k++) vrels[k].multiplyScalar(s);

  const Pspall = new THREE.Vector3();
  for (let k = 0; k < count; k++) Pspall.add(vrels[k].clone().multiplyScalar(masses[k]));

  const newMass = m0 - mSpall;
  if (newMass <= params.minFractureMass * 0.25) return;

  const P0 = baseVel.clone().multiplyScalar(m0);
  const Pnew = P0.clone().sub(Pspall);
  const vNew = Pnew.multiplyScalar(1 / newMass);

  mass[idx] = newMass;
  vel[idx * 3] = vNew.x; vel[idx * 3 + 1] = vNew.y; vel[idx * 3 + 2] = vNew.z;
  radius[idx] = approxRadiusFromMass(newMass, mat.rho);

  damage[idx] += 0.15;

  for (let k = 0; k < count; k++) {
    const m = masses[k];
    const rr = approxRadiusFromMass(m, mat.rho);

    const dir = vrels[k].clone().normalize();
    const spawnPos = basePos.clone()
      .add(outwardNormal.clone().multiplyScalar(r0 * (0.85 + 0.12 * Math.random())))
      .add(dir.clone().multiplyScalar(r0 * 0.15 * Math.random()));

    const spawnVel = baseVel.clone().add(vrels[k]);

    emit({
      pos: spawnPos,
      vel: spawnVel,
      mass: m,
      radius: rr,
      materialId: materialId[idx],
      temperatureBoost: params.spallTempBoost
    });
  }
}

function fragmentBody(idx, outwardNormal, Eimp, bodies, mat, params, emit) {
  const { pos, vel, mass, radius, active, materialId, damage } = bodies;
  if (!active[idx]) return;

  const m0 = mass[idx];
  const r0 = radius[idx];
  const basePos = new THREE.Vector3(pos[idx * 3], pos[idx * 3 + 1], pos[idx * 3 + 2]);
  const baseVel = new THREE.Vector3(vel[idx * 3], vel[idx * 3 + 1], vel[idx * 3 + 2]);

  const budget = params.maxFragmentsPerEvent;
  const count = Math.max(8, Math.min(budget, Math.floor(10 + 22 * Math.min(1, damage[idx]))));

  const keepCore = (m0 > params.coreKeepMass && Math.random() < params.coreKeepChance);
  const coreFrac = keepCore ? params.coreKeepFrac : 0.0;

  const mCore = m0 * coreFrac;
  const mFragTotal = m0 - mCore;
  if (mFragTotal < params.minFragmentMass * 6) return;

  const masses = powerLawMasses(mFragTotal, count, mat.fragmentAlpha, 0.002);

  const Eej = params.ejectaEnergyFracFrag * Eimp * (0.35 + 0.65 * Math.random());

  const vrels = [];
  let pSum = new THREE.Vector3();
  const cone = params.fragmentConeAngleRad;

  for (let k = 0; k < count; k++) {
    const dir = randUnitVecBiased(outwardNormal, cone);
    const speed = params.fragmentSpeedBase
      * (0.6 + 1.2 * Math.random())
      * Math.pow(masses[k] / (mFragTotal / count), -0.18);
    const vrel = dir.multiplyScalar(speed);
    vrels.push(vrel);
    pSum.add(vrel.clone().multiplyScalar(masses[k]));
  }

  const vCorr = pSum.multiplyScalar(1 / mFragTotal);
  for (let k = 0; k < count; k++) vrels[k].sub(vCorr);

  let curE = 0;
  for (let k = 0; k < count; k++) curE += 0.5 * masses[k] * vrels[k].lengthSq();
  const s = Math.sqrt(Math.max(1e-9, Eej / (curE + 1e-9)));
  for (let k = 0; k < count; k++) vrels[k].multiplyScalar(s);

  const Pfr = new THREE.Vector3();
  for (let k = 0; k < count; k++) Pfr.add(vrels[k].clone().multiplyScalar(masses[k]));

  const P0 = baseVel.clone().multiplyScalar(m0);
  const Pcore = P0.clone().sub(Pfr);

  if (keepCore && mCore > params.minFractureMass) {
    mass[idx] = mCore;
    const vCore = Pcore.multiplyScalar(1 / mCore);
    vel[idx * 3] = vCore.x; vel[idx * 3 + 1] = vCore.y; vel[idx * 3 + 2] = vCore.z;
    radius[idx] = approxRadiusFromMass(mCore, mat.rho);
    damage[idx] = 0.35;
  } else {
    active[idx] = 0;
  }

  for (let k = 0; k < count; k++) {
    const m = masses[k];
    if (m < params.minFragmentMass) continue;
    const rr = approxRadiusFromMass(m, mat.rho);

    const dir = vrels[k].clone().normalize();
    const spawnPos = basePos.clone()
      .add(dir.clone().multiplyScalar(r0 * (0.2 + 0.8 * Math.random())))
      .add(outwardNormal.clone().multiplyScalar(r0 * 0.25 * Math.random()));

    const spawnVel = baseVel.clone().add(vrels[k]);

    emit({
      pos: spawnPos,
      vel: spawnVel,
      mass: m,
      radius: rr,
      materialId: materialId[idx],
      temperatureBoost: params.fragmentTempBoost
    });
  }
}

export function tryReaccumulate(bodies, params, pairs, markMerge) {
  if (!params.reaccumulationEnabled) return;

  const { pos, vel, mass, radius, active, temperature, vaporFrac } = bodies;

  for (let p = 0; p < pairs.length; p++) {
    const [i, j] = pairs[p];
    if (!active[i] || !active[j]) continue;

    if (vaporFrac[i] > 0.2 || vaporFrac[j] > 0.2) continue;

    const ix = pos[i * 3], iy = pos[i * 3 + 1], iz = pos[i * 3 + 2];
    const jx = pos[j * 3], jy = pos[j * 3 + 1], jz = pos[j * 3 + 2];
    const dx = jx - ix, dy = jy - iy, dz = jz - iz;

    const d2 = dx * dx + dy * dy + dz * dz;
    const rSum = radius[i] + radius[j];

    const stickDist = rSum * params.reaccumulationDistanceFactor;
    if (d2 > stickDist * stickDist) continue;

    const ivx = vel[i * 3], ivy = vel[i * 3 + 1], ivz = vel[i * 3 + 2];
    const jvx = vel[j * 3], jvy = vel[j * 3 + 1], jvz = vel[j * 3 + 2];
    const rvx = jvx - ivx, rvy = jvy - ivy, rvz = jvz - ivz;
    const vRel = Math.sqrt(rvx * rvx + rvy * rvy + rvz * rvz);

    if (vRel > params.reaccumulationMaxRelSpeed) continue;

    const heat = (temperature[i] + temperature[j]) * 0.5;
    if (heat > params.reaccumulationMaxTemp) continue;

    markMerge(i, j);
  }
}
