// gravityIntegrator.js
// 현실성 원칙:
// - 완전 FEM/연소화학 대신 "물리량 기반 근사 모델"로 실시간 설득력 확보
// - 질량/운동량/각운동량 보존을 우선, 열/발광 등은 파생 상태변수로 연동
// - Euler 금지, Velocity Verlet/Leapfrog 계열 사용

import * as THREE from "three";

/**
 * computeAccelerations()
 * a_i = G Σ_{j≠i} m_j (r_j - r_i) / ((r^2 + eps^2)^(3/2))
 */
export function computeAccelerations(pos, mass, active, n, G, eps, accOut) {
  const e2 = eps * eps;
  for (let i = 0; i < n * 3; i++) accOut[i] = 0;

  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;

    const ix = pos[i * 3], iy = pos[i * 3 + 1], iz = pos[i * 3 + 2];
    let ax = 0, ay = 0, az = 0;

    for (let j = 0; j < n; j++) {
      if (i === j || !active[j]) continue;

      const dx = pos[j * 3] - ix;
      const dy = pos[j * 3 + 1] - iy;
      const dz = pos[j * 3 + 2] - iz;
      const r2 = dx * dx + dy * dy + dz * dz + e2;
      const invR = 1 / Math.sqrt(r2);
      const invR3 = invR * invR * invR;

      const s = G * mass[j] * invR3;
      ax += dx * s; ay += dy * s; az += dz * s;
    }

    accOut[i * 3] = ax;
    accOut[i * 3 + 1] = ay;
    accOut[i * 3 + 2] = az;
  }
}

/**
 * integrateVerlet()
 * (참고용) Velocity Verlet:
 * v += a*dt/2; x += v*dt; aNew=acc(x); a=aNew; v += a*dt/2
 */
export function integrateVerlet(pos, vel, acc, accNew, active, n, dt) {
  const h = dt * 0.5;

  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    vel[k]     += acc[k] * h;
    vel[k + 1] += acc[k + 1] * h;
    vel[k + 2] += acc[k + 2] * h;
  }

  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    pos[k]     += vel[k] * dt;
    pos[k + 1] += vel[k + 1] * dt;
    pos[k + 2] += vel[k + 2] * dt;
  }

  for (let i = 0; i < n * 3; i++) acc[i] = accNew[i];

  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    vel[k]     += acc[k] * h;
    vel[k + 1] += acc[k + 1] * h;
    vel[k + 2] += acc[k + 2] * h;
  }
}

/**
 * estimateSystemInvariants()
 * - K는 정확, U는 pair 샘플링 근사 (O(N^2) 회피)
 * - L은 정확 합산
 */
export function estimateSystemInvariants(pos, vel, mass, active, n, G, eps, samplePairs = 900) {
  let K = 0;
  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    const vx = vel[k], vy = vel[k + 1], vz = vel[k + 2];
    K += 0.5 * mass[i] * (vx * vx + vy * vy + vz * vz);
  }

  let U = 0;
  const e2 = eps * eps;
  const alive = [];
  for (let i = 0; i < n; i++) if (active[i]) alive.push(i);

  const m = alive.length;
  if (m > 1) {
    const pairs = Math.min(samplePairs, (m * (m - 1)) / 2);
    for (let s = 0; s < pairs; s++) {
      const a = alive[(Math.random() * m) | 0];
      let b = alive[(Math.random() * m) | 0];
      if (b === a) b = alive[(b + 1) % m];

      const ax = pos[a * 3], ay = pos[a * 3 + 1], az = pos[a * 3 + 2];
      const dx = pos[b * 3] - ax;
      const dy = pos[b * 3 + 1] - ay;
      const dz = pos[b * 3 + 2] - az;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz + e2);

      U += -G * mass[a] * mass[b] / r;
    }
    const totalPairs = (m * (m - 1)) / 2;
    U *= totalPairs / pairs;
  }

  const L = new THREE.Vector3();
  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    const r = new THREE.Vector3(pos[k], pos[k + 1], pos[k + 2]);
    const p = new THREE.Vector3(vel[k], vel[k + 1], vel[k + 2]).multiplyScalar(mass[i]);
    L.add(r.clone().cross(p));
  }

  return { K, U, E: K + U, L };
}
