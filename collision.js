// collision.js
// Broadphase + Narrowphase + Impact Scalars

import * as THREE from "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js";

/**
 * Spatial hash broadphase for sphere bodies.
 * cellSize는 평균/최대 반지름 기반으로 주는 것이 일반적.
 */
export function broadphaseSpatialHash(pos, radius, active, n, cellSize) {
  const inv = 1 / cellSize;
  const map = new Map(); // key -> [indices]
  const keyOf = (cx, cy, cz) => `${cx},${cy},${cz}`;

  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    const cx = Math.floor(pos[k] * inv);
    const cy = Math.floor(pos[k + 1] * inv);
    const cz = Math.floor(pos[k + 2] * inv);
    const key = keyOf(cx, cy, cz);
    let arr = map.get(key);
    if (!arr) { arr = []; map.set(key, arr); }
    arr.push(i);
  }

  const pairs = [];
  const neighbors = [-1, 0, 1];

  for (const [key, arr] of map.entries()) {
    const [sx, sy, sz] = key.split(",").map(Number);

    for (const dx of neighbors) for (const dy of neighbors) for (const dz of neighbors) {
      const key2 = keyOf(sx + dx, sy + dy, sz + dz);
      const arr2 = map.get(key2);
      if (!arr2) continue;

      // pair 생성 (중복 최소화)
      for (let a = 0; a < arr.length; a++) {
        const i = arr[a];
        for (let b = 0; b < arr2.length; b++) {
          const j = arr2[b];
          if (j <= i) continue;
          // quick radius check with loose bound by cell adjacency already
          pairs.push([i, j]);
        }
      }
    }
  }

  return pairs;
}

/**
 * narrowphase sphere-sphere, returns collisions with contact normal and penetration
 */
export function narrowphaseSpheres(pos, radius, active, pairs) {
  const hits = [];
  for (let p = 0; p < pairs.length; p++) {
    const [i, j] = pairs[p];
    if (!active[i] || !active[j]) continue;

    const ix = pos[i * 3], iy = pos[i * 3 + 1], iz = pos[i * 3 + 2];
    const jx = pos[j * 3], jy = pos[j * 3 + 1], jz = pos[j * 3 + 2];

    const dx = jx - ix, dy = jy - iy, dz = jz - iz;
    const dist2 = dx * dx + dy * dy + dz * dz;
    const rSum = radius[i] + radius[j];

    if (dist2 <= rSum * rSum) {
      const dist = Math.sqrt(dist2) || 1e-8;
      const nx = dx / dist, ny = dy / dist, nz = dz / dist;
      const pen = rSum - dist;
      hits.push({ i, j, n: new THREE.Vector3(nx, ny, nz), dist, pen });
    }
  }
  return hits;
}

/**
 * computeImpactScalars()
 * - 상대속도/충돌각/감소질량/충돌에너지/특정에너지 Q
 */
export function computeImpactScalars(hit, pos, vel, mass) {
  const { i, j, n } = hit;

  const iv = new THREE.Vector3(vel[i * 3], vel[i * 3 + 1], vel[i * 3 + 2]);
  const jv = new THREE.Vector3(vel[j * 3], vel[j * 3 + 1], vel[j * 3 + 2]);

  const rel = jv.clone().sub(iv);
  const vRel = rel.length();
  const vN = rel.dot(n); // normal component (signed)
  const vT = rel.clone().sub(n.clone().multiplyScalar(vN)).length();

  const m1 = mass[i], m2 = mass[j];
  const mu = (m1 * m2) / (m1 + m2);
  const Eimp = 0.5 * mu * vRel * vRel;
  const Q = Eimp / (m1 + m2);

  // contact point (mid along normal)
  const ip = new THREE.Vector3(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);
  const jp = new THREE.Vector3(pos[j * 3], pos[j * 3 + 1], pos[j * 3 + 2]);
  const cp = ip.clone().add(jp).multiplyScalar(0.5);

  return {
    i, j,
    n: n.clone(),
    relVel: rel,
    vRel,
    vN,
    vT,
    mu,
    Eimp,
    Q,
    contactPoint: cp
  };
}

/**
 * collisionImpulseResponse()
 * - 운동량 보존 기반 반발/마찰 근사(열/용융 상태에 따라 탄성↓)
 * - 각운동량 전달(간이): 접선 임펄스가 스핀을 변화시킴 (I 근사)
 */
export function collisionImpulseResponse(impact, bodies, params) {
  const { i, j, n, vN } = impact;
  const { pos, vel, mass, radius, omega, active, meltFrac, materialId } = bodies;
  if (!active[i] || !active[j]) return { impulseMag: 0, frictionImpulse: 0 };

  const m1 = mass[i], m2 = mass[j];
  const inv1 = 1 / m1, inv2 = 1 / m2;

  // Restitution: 재료+용융에 따라 감소
  const eBase = params.materials[materialId[i]].restitution * params.materials[materialId[j]].restitution;
  const soft = 1 - 0.65 * Math.min(1, (meltFrac[i] + meltFrac[j]) * 0.5);
  const e = Math.max(0.02, Math.min(0.85, eBase * soft));

  // Normal impulse magnitude
  const jn = -(1 + e) * vN / (inv1 + inv2);
  if (!Number.isFinite(jn) || jn <= 0) return { impulseMag: 0, frictionImpulse: 0 };

  // Apply linear impulses
  vel[i * 3]     -= (jn * n.x) * inv1;
  vel[i * 3 + 1] -= (jn * n.y) * inv1;
  vel[i * 3 + 2] -= (jn * n.z) * inv1;

  vel[j * 3]     += (jn * n.x) * inv2;
  vel[j * 3 + 1] += (jn * n.y) * inv2;
  vel[j * 3 + 2] += (jn * n.z) * inv2;

  // Friction impulse (Coulomb) on tangential relative velocity
  const iv = new THREE.Vector3(vel[i * 3], vel[i * 3 + 1], vel[i * 3 + 2]);
  const jv = new THREE.Vector3(vel[j * 3], vel[j * 3 + 1], vel[j * 3 + 2]);
  const rel = jv.clone().sub(iv);
  const vN2 = rel.dot(n);
  const vtVec = rel.clone().sub(n.clone().multiplyScalar(vN2));
  const vt = vtVec.length();

  let jt = 0;
  if (vt > 1e-6) {
    const muF = params.materials[materialId[i]].friction * params.materials[materialId[j]].friction;
    const maxJt = muF * jn;
    jt = Math.min(maxJt, vt / (inv1 + inv2)); // simplistic
    const t = vtVec.multiplyScalar(1 / vt);

    vel[i * 3]     += (jt * t.x) * inv1;
    vel[i * 3 + 1] += (jt * t.y) * inv1;
    vel[i * 3 + 2] += (jt * t.z) * inv1;

    vel[j * 3]     -= (jt * t.x) * inv2;
    vel[j * 3 + 1] -= (jt * t.y) * inv2;
    vel[j * 3 + 2] -= (jt * t.z) * inv2;

    // Spin transfer (approx): τ ≈ r × J_t
    // sphere inertia I = 0.4 m r^2
    const ri = radius[i], rj = radius[j];
    const Ii = 0.4 * m1 * ri * ri;
    const Ij = 0.4 * m2 * rj * rj;
    const JtVec = new THREE.Vector3(t.x, t.y, t.z).multiplyScalar(jt);

    // contact lever arms: ±n * r
    const tauI = new THREE.Vector3().copy(n).multiplyScalar(ri).cross(JtVec);
    const tauJ = new THREE.Vector3().copy(n).multiplyScalar(-rj).cross(JtVec.clone().multiplyScalar(-1));

    omega[i * 3]     += tauI.x / Ii;
    omega[i * 3 + 1] += tauI.y / Ii;
    omega[i * 3 + 2] += tauI.z / Ii;

    omega[j * 3]     += tauJ.x / Ij;
    omega[j * 3 + 1] += tauJ.y / Ij;
    omega[j * 3 + 2] += tauJ.z / Ij;
  }

  // Separate interpenetration (position correction, mild)
  const k = params.collisionPositionCorrection;
  const corr = Math.max(0, impact.pen) * k;
  if (corr > 0) {
    pos[i * 3]     -= n.x * corr * inv1 / (inv1 + inv2);
    pos[i * 3 + 1] -= n.y * corr * inv1 / (inv1 + inv2);
    pos[i * 3 + 2] -= n.z * corr * inv1 / (inv1 + inv2);

    pos[j * 3]     += n.x * corr * inv2 / (inv1 + inv2);
    pos[j * 3 + 1] += n.y * corr * inv2 / (inv1 + inv2);
    pos[j * 3 + 2] += n.z * corr * inv2 / (inv1 + inv2);
  }

  return { impulseMag: jn, frictionImpulse: jt };
}
