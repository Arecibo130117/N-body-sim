// gravityIntegrator.js
// "완전한 FEM이나 연소화학을 요구하지 않는다. 대신 실시간에서 가능한 물리량 기반 근사 모델로 구현하되,
// 결과가 시각적으로도 물리적으로도 설득력 있게 보이도록 한다."
// "가능한 한 보존법칙(질량/운동량/각운동량)을 우선 지키고 ..."

import * as THREE from "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js";

/**
 * computeAccelerations()
 * 만유인력 (softening 포함):
 * a_i = G Σ_{j≠i} m_j (r_j - r_i) / ( (r^2 + eps^2)^(3/2) )
 *
 * - 기본은 O(N^2) brute force
 * - N이 커지면 샘플링/LOD에서 자연스럽게 타협(본 프로젝트는 실시간 구현 지향)
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
 * Euler 금지. Velocity Verlet:
 * v(t+dt/2)=v(t)+a(t)*dt/2
 * x(t+dt)=x(t)+v(t+dt/2)*dt
 * a(t+dt)=f(x(t+dt))
 * v(t+dt)=v(t+dt/2)+a(t+dt)*dt/2
 */
export function integrateVerlet(pos, vel, acc, accNew, active, n, dt) {
  const h = dt * 0.5;

  // v += a*dt/2
  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    vel[k]     += acc[k] * h;
    vel[k + 1] += acc[k + 1] * h;
    vel[k + 2] += acc[k + 2] * h;
  }

  // x += v*dt
  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    pos[k]     += vel[k] * dt;
    pos[k + 1] += vel[k + 1] * dt;
    pos[k + 2] += vel[k + 2] * dt;
  }

  // a <- accNew
  for (let i = 0; i < n * 3; i++) acc[i] = accNew[i];

  // v += a*dt/2
  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    vel[k]     += acc[k] * h;
    vel[k + 1] += acc[k + 1] * h;
    vel[k + 2] += acc[k + 2] * h;
  }
}

/**
 * 에너지/각운동량 모니터링용(비용 절약 위해 샘플링 지원)
 */
export function estimateSystemInvariants(pos, vel, mass, active, n, G, eps, samplePairs = 800) {
  let K = 0;
  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    const vx = vel[k], vy = vel[k + 1], vz = vel[k + 2];
    K += 0.5 * mass[i] * (vx * vx + vy * vy + vz * vz);
  }

  // Potential: 랜덤 pair 샘플로 근사 (N^2 풀계산 대신)
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
    // 샘플링 보정(대략)
    const totalPairs = (m * (m - 1)) / 2;
    U *= totalPairs / pairs;
  }

  // Angular momentum
  const L = new THREE.Vector3();
  for (let i = 0; i < n; i++) {
    if (!active[i]) continue;
    const k = i * 3;
    const r = new THREE.Vector3(pos[k], pos[k + 1], pos[k + 2]);
    const p = new THREE.Vector3(vel[k], vel[k + 1], vel[k + 2]).multiplyScalar(mass[i]);
    L.add(r.cross(p));
  }

  return { K, U, E: K + U, L };
}
