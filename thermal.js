// thermal.js
// ΔH = η Eimp, ΔT = ΔH/(m c)
// meltFraction, vaporFraction, emissiveLevel -> Bloom 연동
// 파티클(불꽃/연기/플래시)는 상태변수 기반으로 발생

export function updateThermalState(bodies, params, dt) {
  const { active, mass, temperature, heat, meltFrac, vaporFrac, emissive, materialId } = bodies;
  const mats = params.materials;

  for (let i = 0; i < bodies.capacity; i++) {
    if (!active[i]) continue;

    const mat = mats[materialId[i]];
    const c = mat.heatCapacity;

    temperature[i] = Math.max(params.ambientTemp, params.ambientTemp + heat[i] / (mass[i] * c));

    const T = temperature[i];
    if (T > mat.meltTemp) {
      const t = Math.min(1, (T - mat.meltTemp) / (mat.vaporTemp - mat.meltTemp + 1e-9));
      meltFrac[i] = Math.max(meltFrac[i], t);
    }
    if (T > mat.vaporTemp) {
      const t = Math.min(1, (T - mat.vaporTemp) / (mat.plasmaTemp - mat.vaporTemp + 1e-9));
      vaporFrac[i] = Math.max(vaporFrac[i], t);
    }

    const cool = mat.cooling * (0.6 + 1.2 * vaporFrac[i]);
    heat[i] *= Math.exp(-cool * dt);

    const tNorm = Math.max(0, (T - mat.glowTempStart) / (mat.plasmaTemp - mat.glowTempStart + 1e-9));
    emissive[i] = Math.min(1, tNorm) * mat.emissiveScale;
  }
}

export function applyImpactHeating(impact, bodies, params) {
  const { i, j, Eimp } = impact;
  const { active, heat, materialId } = bodies;
  if (!active[i] || !active[j]) return { dHi: 0, dHj: 0 };

  const matI = params.materials[materialId[i]];
  const matJ = params.materials[materialId[j]];

  const eta = params.heatEtaBase * (0.8 + 0.4 * Math.random());

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

export function spawnFireSmokeParticles(impact, bodies, params, spawnParticleBurst) {
  const { i, j, contactPoint, n, Eimp } = impact;
  const { active, temperature, vaporFrac, emissive } = bodies;
  if (!active[i] || !active[j]) return;

  const Ti = temperature[i], Tj = temperature[j];
  const vi = vaporFrac[i], vj = vaporFrac[j];
  const ei = emissive[i], ej = emissive[j];

  const T = 0.5 * (Ti + Tj);
  const v = 0.5 * (vi + vj);
  const e = 0.5 * (ei + ej);

  const flashN = Math.floor(params.flashParticlesBase * (0.4 + 1.6 * e));
  spawnParticleBurst({
    kind: "flash",
    origin: contactPoint,
    axis: n,
    count: Math.min(params.particleBudgetHard, flashN),
    energy: Eimp
  });

  if (T > params.fireTempThreshold || v > 0.08) {
    const fireN = Math.floor(params.fireParticlesBase * (0.3 + 2.2 * e + 1.0 * v));
    spawnParticleBurst({ kind: "fire", origin: contactPoint, axis: n, count: fireN, energy: Eimp });
  }

  if (v > 0.02) {
    const smokeN = Math.floor(params.smokeParticlesBase * (0.4 + 2.5 * v));
    spawnParticleBurst({ kind: "smoke", origin: contactPoint, axis: n, count: smokeN, energy: Eimp });
  }
}
