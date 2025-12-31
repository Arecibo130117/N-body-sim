// fracture.js
// damage 누적 → 임계치 초과 시 fragmentation
// 파편 질량분포: 파워법칙 근사 + (질량/운동량/각운동량) 보존되게 배분
// 박리(spallation): 표면층 파편 다수 + 크레이터/노멀 변형(LOD-1 시각)
// 재응집: 파편 군집이 안정적으로 붙으면 rubble pile → 병합(옵션)

import * as THREE from "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js";

function randUnitVecBiased(coneAxis, coneAngleRad) {
  // 샘플링: 축 주변 콘 샘플
  const u = Math.random();
  const v = Math.random();
  const cosA = Math.cos(coneAngleRad);
  const z = cosA + (1 - cosA) * u; // [cosA,1]
  const phi = 2 * Math.PI * v;
  const r = Math.sqrt(Math.max(0, 1 - z * z));
  const local = new THREE.Vector3(r * Math.cos(phi), r * Math.sin(phi), z);

  // 로컬 z축을 coneAxis로 회전
  const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), coneAxis.clone().normalize());
  return local.applyQuaternion(q).normalize();
}

function powerLawMasses(mTotal, count, alpha = 1.8, minFrac = 0.002) {
  // m ~ x^(-alpha) 근사. 간단히: 랜덤->가중
  const masses = new Array(count);
  let sum = 0;
  for (let i = 0; i < count; i++) {
    const x = Math.random();
    const w = Math.pow(Math.max(1e-6, x), -1 / (alpha)); // heavy-tail
    masses[i] = w;
    sum += w;
  }
  // normalize and enforce minimum
  const minM = mTotal * minFrac / count;
  let sum2 = 0;
  for (let i = 0; i < count; i++) {
    masses[i] = Math.max(minM, (masses[i] / sum) * mTotal);
    sum2 += masses[i];
  }
  // renormalize to exact mTotal
  const scale = mTotal / sum2;
  for (let i = 0; i < count; i++) masses[i] *= scale;
  return masses;
}

function approxRadiusFromMass(m, rho) {
  // sphere volume: m = rho * 4/3 π r^3
  return Math.cbrt(m / (rho * (4 / 3) * Math.PI));
}

/**
 * fractureAndSpawnFragments()
 * - Q (=Eimp/(m1+m2)) vs Q* 비교로 damage 증가
 * - damage>1이면 fragmentation 수행
 * - 파편 운동량/에너지 배분은 "보존량 우선" + "근사 모델"
 */
export function fractureAndSpawnFragments(impact, bodies, params, emit) {
  const { i, j, n, Eimp, Q, contactPoint, relVel } = impact;
  const { mass, radius, damage, active, materialId, pos, vel, omega } = bodies;

  const matI = params.materials[materialId[i]];
  const matJ = params.materials[materialId[j]];

  // 강도 임계(간이): Q*를 재료별로 사용
  const QsI = matI.Qstar;
  const QsJ = matJ.Qstar;

  // damage 증가량: Q/Q*
  const dI = Q / (QsI + 1e-9);
  const dJ = Q / (QsJ + 1e-9);

  damage[i] += dI * params.damageGain;
  damage[j] += dJ * params.damageGain;

  // spallation: 표면 박리(충돌이 충분히 크면 얇은층 다수 파편)
  const doSpall = Q > 0.35 * Math.min(QsI, QsJ);

  // Fragmentation check
  const shouldFracI = damage[i] > 1 && mass[i] > params.minFractureMass;
  const shouldFracJ = damage[j] > 1 && mass[j] > params.minFractureMass;

  // Impact flash / ejecta cone 파티클은 render 쪽에서 처리하게 이벤트 힌트 제공
  const impactNormal = n.clone();
  const impactSpeed = relVel.length();

  if (doSpall) {
    // 표면층 질량 비율
    spallateBody(i, impactNormal.clone().negate(), contactPoint, impactSpeed, bodies, matI, params, emit);
    spallateBody(j, impactNormal.clone(), contactPoint, impactSpeed, bodies, matJ, params, emit);
  }

  if (shouldFracI) fragmentBody(i, impactNormal.clone().negate(), Eimp, bodies, matI, params, emit);
  if (shouldFracJ) fragmentBody(j, impactNormal.clone(), Eimp, bodies, matJ, params, emit);
}

/**
 * 박리(spallation): 표면층에서 작은 파편 다수 + "크레이터 느낌"은 LOD-1 변형 파라미터로 전달
 */
function spallateBody(idx, outwardNormal, contactPoint, impactSpeed, bodies, mat, params, emit) {
  const { pos, vel, mass, radius, active, materialId, omega, damage } = bodies;
  if (!active[idx]) return;

  const m0 = mass[idx];
  const r0 = radius[idx];

  const layerFrac = Math.min(0.08, params.spallLayerFrac * (0.25 + Math.min(1, impactSpeed / 40)));
  const mSpall = m0 * layerFrac;
  if (mSpall < params.minFragmentMass * 4) return;

  // 개수는 충돌 세기에 따라
  const count = Math.min(params.maxSpallFragments, Math.max(6, (mSpall / params.minFragmentMass) | 0));
  const masses = powerLawMasses(mSpall, count, 2.2, 0.001);

  // COM frame: 파편의 총 운동량 0되도록 방향 샘플 후 보정
  const basePos = new THREE.Vector3(pos[idx * 3], pos[idx * 3 + 1], pos[idx * 3 + 2]);
  const baseVel = new THREE.Vector3(vel[idx * 3], vel[idx * 3 + 1], vel[idx * 3 + 2]);

  // ejecta energy budget (일부만)
  const Eej = params.ejectaEnergyFracSpall * 0.5 * mSpall * impactSpeed * impactSpeed;

  const dirs = [];
  const vrels = [];
  let pSum = new THREE.Vector3();

  const cone = params.spallConeAngleRad;
  for (let k = 0; k < count; k++) {
    const dir = randUnitVecBiased(outwardNormal, cone);
    dirs.push(dir);

    // 초기 속도: 작은 파편일수록 더 빠름
    const speed = (params.spallSpeedBase + Math.random() * params.spallSpeedJitter) * (0.6 + 0.4 * Math.random());
    const vrel = dir.clone().multiplyScalar(speed);
    vrels.push(vrel);
    pSum.add(vrel.clone().multiplyScalar(masses[k]));
  }

  // 총 상대 운동량 0으로 보정
  const vCorr = pSum.multiplyScalar(1 / mSpall);
  for (let k = 0; k < count; k++) vrels[k].sub(vCorr);

  // 상대 운동에너지 맞추기(스케일)
  let curE = 0;
  for (let k = 0; k < count; k++) curE += 0.5 * masses[k] * vrels[k].lengthSq();
  const s = Math.sqrt(Math.max(1e-9, Eej / (curE + 1e-9)));
  for (let k = 0; k < count; k++) vrels[k].multiplyScalar(s);

  // 본체 질량/운동량 차감(보존)
  const Pspall = new THREE.Vector3();
  for (let k = 0; k < count; k++) Pspall.add(vrels[k].clone().multiplyScalar(masses[k]));
  // 본체는 반작용으로 속도 변함
  const newMass = m0 - mSpall;
  if (newMass <= params.minFractureMass * 0.25) return;

  // 본체 선형 운동량 보존
  const P0 = baseVel.clone().multiplyScalar(m0);
  const Pnew = P0.clone().sub(Pspall);
  const vNew = Pnew.multiplyScalar(1 / newMass);

  mass[idx] = newMass;
  vel[idx * 3] = vNew.x; vel[idx * 3 + 1] = vNew.y; vel[idx * 3 + 2] = vNew.z;

  // 반지름 업데이트(밀도 유지)
  radius[idx] = approxRadiusFromMass(newMass, mat.rho);

  // 크레이터/찌그러짐 느낌(LOD-1에서 렌더에 반영하도록 damage를 조금 더)
  damage[idx] += 0.15;

  // 파편 생성
  for (let k = 0; k < count; k++) {
    const m = masses[k];
    const rr = approxRadiusFromMass(m, mat.rho);

    // 표면 근처에서 발사
    const spawnPos = basePos.clone()
      .add(outwardNormal.clone().multiplyScalar(r0 * (0.85 + 0.12 * Math.random())))
      .add(dirs[k].clone().multiplyScalar(r0 * 0.15 * Math.random()));

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
  const { pos, vel, mass, radius, active, materialId, omega, damage } = bodies;
  if (!active[idx]) return;

  const m0 = mass[idx];
  const r0 = radius[idx];
  const basePos = new THREE.Vector3(pos[idx * 3], pos[idx * 3 + 1], pos[idx * 3 + 2]);
  const baseVel = new THREE.Vector3(vel[idx * 3], vel[idx * 3 + 1], vel[idx * 3 + 2]);

  // 파편 개수
  const budget = params.maxFragmentsPerEvent;
  const count = Math.max(8, Math.min(budget, Math.floor(10 + 22 * Math.min(1, damage[idx]))));

  // 본체 일부는 큰 조각으로 남길 수도(완전 분해 방지)
  const keepCore = (m0 > params.coreKeepMass && Math.random() < params.coreKeepChance);
  const coreFrac = keepCore ? params.coreKeepFrac : 0.0;

  const mCore = m0 * coreFrac;
  const mFragTotal = m0 - mCore;
  if (mFragTotal < params.minFragmentMass * 6) return;

  const masses = powerLawMasses(mFragTotal, count, mat.fragmentAlpha, 0.002);

  // Ejecta 에너지 예산 (충돌 에너지 일부)
  const Eej = params.ejectaEnergyFracFrag * Eimp * (0.35 + 0.65 * Math.random());

  const vrels = [];
  let pSum = new THREE.Vector3();
  const cone = params.fragmentConeAngleRad;

  for (let k = 0; k < count; k++) {
    const dir = randUnitVecBiased(outwardNormal, cone);
    // 작은 조각일수록 더 빠름(질량 반비례 근사)
    const speed = params.fragmentSpeedBase * (0.6 + 1.2 * Math.random()) * Math.pow(masses[k] / (mFragTotal / count), -0.18);
    const vrel = dir.multiplyScalar(speed);
    vrels.push(vrel);
    pSum.add(vrel.clone().multiplyScalar(masses[k]));
  }

  // 총 상대 운동량 0
  const vCorr = pSum.multiplyScalar(1 / mFragTotal);
  for (let k = 0; k < count; k++) vrels[k].sub(vCorr);

  // 에너지 맞추기
  let curE = 0;
  for (let k = 0; k < count; k++) curE += 0.5 * masses[k] * vrels[k].lengthSq();
  const s = Math.sqrt(Math.max(1e-9, Eej / (curE + 1e-9)));
  for (let k = 0; k < count; k++) vrels[k].multiplyScalar(s);

  // 보존: 파편 총 운동량을 원래 질량에 반영해 core 또는 제거
  const Pfr = new THREE.Vector3();
  for (let k = 0; k < count; k++) Pfr.add(vrels[k].clone().multiplyScalar(masses[k]));

  const P0 = baseVel.clone().multiplyScalar(m0);
  const Pcore = P0.clone().sub(Pfr);

  // 기존 바디는 core로 남기거나 비활성화
  if (keepCore && mCore > params.minFractureMass) {
    mass[idx] = mCore;
    const vCore = Pcore.multiplyScalar(1 / mCore);
    vel[idx * 3] = vCore.x; vel[idx * 3 + 1] = vCore.y; vel[idx * 3 + 2] = vCore.z;
    radius[idx] = approxRadiusFromMass(mCore, mat.rho);
    damage[idx] = 0.35; // core는 손상 잔존
  } else {
    active[idx] = 0;
  }

  // 파편 생성
  for (let k = 0; k < count; k++) {
    const m = masses[k];
    if (m < params.minFragmentMass) continue;
    const rr = approxRadiusFromMass(m, mat.rho);
    const dir = vrels[k].clone().normalize();
    const spawnPos = basePos.clone()
      .add(dir.multiplyScalar(r0 * (0.2 + 0.8 * Math.random())))
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

/**
 * re-accretion: 가까이 있고 상대속도 낮으면 병합(옵션)
 * - 병합 시 질량/운동량/각운동량 보존
 * - 초기 형상은 찌그러진 타원체(렌더에서 scale로 표현)
 */
export function tryReaccumulate(bodies, params, pairs, markMerge) {
  if (!params.reaccumulationEnabled) return;

  const { pos, vel, mass, radius, active, temperature, vaporFrac } = bodies;

  for (let p = 0; p < pairs.length; p++) {
    const [i, j] = pairs[p];
    if (!active[i] || !active[j]) continue;

    // 뜨거운 가스화 상태는 병합 금지
    if (vaporFrac[i] > 0.2 || vaporFrac[j] > 0.2) continue;

    const ix = pos[i * 3], iy = pos[i * 3 + 1], iz = pos[i * 3 + 2];
    const jx = pos[j * 3], jy = pos[j * 3 + 1], jz = pos[j * 3 + 2];
    const dx = jx - ix, dy = jy - iy, dz = jz - iz;

    const d2 = dx * dx + dy * dy + dz * dz;
    const rSum = radius[i] + radius[j];

    // 약간 더 가까워야 "러블파일"처럼 달라붙음
    const stickDist = rSum * params.reaccumulationDistanceFactor;
    if (d2 > stickDist * stickDist) continue;

    const ivx = vel[i * 3], ivy = vel[i * 3 + 1], ivz = vel[i * 3 + 2];
    const jvx = vel[j * 3], jvy = vel[j * 3 + 1], jvz = vel[j * 3 + 2];
    const rvx = jvx - ivx, rvy = jvy - ivy, rvz = jvz - ivz;
    const vRel = Math.sqrt(rvx * rvx + rvy * rvy + rvz * rvz);

    if (vRel > params.reaccumulationMaxRelSpeed) continue;

    // 너무 뜨거우면 점착력 낮음
    const heat = (temperature[i] + temperature[j]) * 0.5;
    if (heat > params.reaccumulationMaxTemp) continue;

    markMerge(i, j);
  }
}
