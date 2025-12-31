// thermal.js
// Thermal / Glow / Fire / Smoke (상태변수 연동)
// ΔH = η E_imp, ΔT = ΔH/(m c)
// meltFraction, vaporFraction, emissiveLevel -> Bloom 영향

import * as THREE from "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js";

/**
 * updateThermalState()
 * - 방열/냉각(간이) + emissiveLevel 갱신
 */
export function updateThermalState(bodies, params, dt) {
  const { active, mass, temperature, heat, meltFrac, vaporFrac, emissive, materialId } = bodies;
  const mats = params.materials;

  for (let i = 0; i < bodies.capacity; i++) {
    if (!active[i]) continue;

    const mat = mats[materialId[i]];
    // heat -> temperature (c 이용)
    // heat는 "열에너지(J)" 누적; 온도는 절대값 근사(K)로
    const c = mat.heatCapacity;
    temperature[i] = Math.max(params.ambientTemp, params.ambientTemp + heat[i] / (mass[i] * c));

    // melt/vapor fraction (0..1)
    const T = temperature[i];
    if (T > mat.meltTemp) {
      const t = Math.min(1, (T - mat.meltTemp) / (mat.vaporTemp - mat.meltTemp + 1e-9));
      meltFrac[i] = Math.max(meltFrac[i], t);
    }
    if (T > mat.vaporTemp) {
      const t = Math.min(1, (T - mat.vaporTemp) / (mat.plasmaTemp - mat.vaporTemp + 1e-9));
      vaporFrac[i] = Math.max(vaporFrac[i], t);
    }

    // 방열/냉각(간이): 열은 지수감쇠, 뜨거울수록 더 빨리
    const cool = mat.cooling * (0.6 + 1.2 * vaporFrac[i]);
    heat[i] *= Math.exp(-cool * dt);

    // emissiveLevel: 온도 기반, 재료별 emissiveScale
    const tNorm = Math.max(0, (T - mat.glowTempStart) / (mat.plasmaTemp - mat.glowTempStart + 1e-9));
    emissive[i] = Math.min(1, tNorm) * mat.emissiveScale;
  }
}

/**
 * applyImpactHeating()
 * - 충돌에너지 일부를 열로 전환 (η)
 */
export function applyImpactHeating(impact, bodies, params) {
  const { i, j, Eimp } = impact;
  const { active, heat, materialId } = bodies;
  if (!active[i] || !active[j]) return { dHi: 0, dHj: 0 };

  const matI = params.materials[materialId[i]];
  const matJ = params.materials[materialId[j]];

  // η는 재료별 + 충돌강도에 따라 약간 변동
  const eta = params.heatEtaBase * (0.8 + 0.4 * Math.random());

  // 분배는 질량/열용량을 단순 고려: 여기서는 재료 emissiveScale로 가중
  const wI = 0.5 + 0.5 * matI.heatAbsorb;
  const wJ = 0.5 + 0.5 * matJ.heatAbsorb;
  const sum = wI + wJ;

  const dH = eta * Eimp;
  const dHi = dH * (wI / sum);
  const dHj = dH * (wJ / sum);

  heat[i] += dHi;
  heat[j] += dHj;

  return { dHi, dHj };
}

/**
 * spawnFireSmokeParticles()
 * - 파티클 발생률은 T, vaporFraction, 충돌 직후 flash 등에 연동
 * - 여기서는 "스폰 요청"을 render쪽 파티클 시스템으로 넘길 이벤트 생성
 */
export function spawnFireSmokeParticles(impact, bodies, params, spawnParticleBurst) {
  const { i, j, contactPoint, n, Eimp } = impact;
  const { active, temperature, vaporFrac, emissive, materialId } = bodies;
  if (!active[i] || !active[j]) return;

  const Ti = temperature[i], Tj = temperature[j];
  const vi = vaporFrac[i], vj = vaporFrac[j];
  const ei = emissive[i], ej = emissive[j];

  const T = 0.5 * (Ti + Tj);
  const v = 0.5 * (vi + vj);
  const e = 0.5 * (ei + ej);

  // 충돌 순간 플래시: 고온일수록 강함
  const flashN = Math.floor(params.flashParticlesBase * (0.4 + 1.6 * e));
  spawnParticleBurst({
    kind: "flash",
    origin: contactPoint,
    axis: n,
    count: Math.min(params.particleBudgetHard, flashN),
    energy: Eimp
  });

  // 플라즈마/불꽃: vaporFraction & emissive에 비례
  if (T > params.fireTempThreshold || v > 0.08) {
    const fireN = Math.floor(params.fireParticlesBase * (0.3 + 2.2 * e + 1.0 * v));
    spawnParticleBurst({
      kind: "fire",
      origin: contactPoint,
      axis: n,
      count: fireN,
      energy: Eimp
    });
  }

  // 연기: 기화/파편이 많을수록(=vaporFrac) + 시간이 지나며 잔류
  if (v > 0.02) {
    const smokeN = Math.floor(params.smokeParticlesBase * (0.4 + 2.5 * v));
    spawnParticleBurst({
      kind: "smoke",
      origin: contactPoint,
      axis: n,
      count: smokeN,
      energy: Eimp
    });
  }
}
