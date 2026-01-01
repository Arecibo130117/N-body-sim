import * as THREE from "three";
import { Body, MaterialPresets, MaterialType } from "./Body.js";
import { clamp, randRange, seededRand } from "../utils/math.js";

/**
 * Collision outcomes: bounce, merge, fragment.
 * We use a coarse model driven by relative speed vs escape speed and specific impact energy.
 */

function escapeSpeed(m1, m2, r1, r2, G){
  const mu = (m1+m2);
  const r = (r1+r2);
  return Math.sqrt(2 * G * mu / Math.max(1e-9, r));
}

function specificImpactEnergy(vRel, mTarget, mImp){
  const mu = (mTarget*mImp)/(mTarget+mImp);
  const Q = 0.5*mu*vRel*vRel / Math.max(1e-9, (mTarget+mImp));
  return Q;
}

function temperatureKick(body, dE){
  const preset = MaterialPresets[body.material] || MaterialPresets.rock;
  body.heatJ = Math.max(0, body.heatJ + dE);
  body.temperatureK = 3 + body.heatJ / (Math.max(1e-9, preset.cp * body.mass));
}

export function handleCollisions(sim){
  const bodies = sim.bodies;
  const G = sim.params.G;
  const restitution = sim.params.restitution;

  // Broadphase naive O(N^2) (works for <= ~500). For huge N, use a spatial hash.
  for (let i=0;i<bodies.length;i++){
    const a = bodies[i];
    if (!a.alive) continue;
    for (let j=i+1;j<bodies.length;j++){
      const b = bodies[j];
      if (!b.alive) continue;

      const dx = b.position.x - a.position.x;
      const dy = b.position.y - a.position.y;
      const dz = b.position.z - a.position.z;
      const r = a.radius + b.radius;
      const dist2 = dx*dx + dy*dy + dz*dz;
      if (dist2 > r*r) continue;

      const dist = Math.sqrt(dist2) + 1e-9;
      const n = new THREE.Vector3(dx/dist, dy/dist, dz/dist);

      // Relative velocity along normal
      const rv = b.velocity.clone().sub(a.velocity);
      const vRel = rv.length();

      // Resolve penetration (position correction)
      const pen = (r - dist);
      const totalMass = a.mass + b.mass;
      const aMove = (b.mass/totalMass) * pen;
      const bMove = (a.mass/totalMass) * pen;
      a.position.addScaledVector(n, -aMove);
      b.position.addScaledVector(n,  bMove);

      // Decide outcome
      const vEsc = escapeSpeed(a.mass,b.mass,a.radius,b.radius,G);
      const Q = specificImpactEnergy(vRel, a.mass, b.mass);
      // thresholds tuned for dramatic visuals in game units
      const mergeBias = sim.params.mergeBias;
      const fragBias = sim.params.fragBias;
      const mergeCond = (vRel < vEsc * (1.15 + mergeBias));
      const catastrophic = (Q > (2.5e5 * (1 - fragBias))); // J/kg-ish scale
      const grazing = Math.abs(rv.dot(n)) < 0.25*vRel;

      // Add craters & heat regardless
      const impactStrength = clamp(vRel / Math.max(1e-6, vEsc), 0, 3);
      a.addCrater(n.clone().multiplyScalar(-1), impactStrength);
      b.addCrater(n.clone(), impactStrength);

      const dEheat = 0.35 * 0.5 * (a.mass*b.mass/(a.mass+b.mass)) * vRel*vRel; // portion to heat
      // split heat proportional to mass
      temperatureKick(a, dEheat * (b.mass/totalMass));
      temperatureKick(b, dEheat * (a.mass/totalMass));

      if (mergeCond && !catastrophic && !grazing){
        sim.mergeBodies(a, b);
        continue;
      }

      if (catastrophic){
        sim.fragmentBodies(a, b, n, vRel);
        continue;
      }

      // Otherwise bounce with impulse
      const vn = rv.dot(n);
      if (vn > 0) continue; // separating
      const invMassA = 1/Math.max(1e-9, a.mass);
      const invMassB = 1/Math.max(1e-9, b.mass);
      const jImpulse = -(1+restitution)*vn/(invMassA+invMassB);
      const impulse = n.clone().multiplyScalar(jImpulse);
      a.velocity.addScaledVector(impulse, -invMassA);
      b.velocity.addScaledVector(impulse,  invMassB);

      // Tangential friction-ish to reduce jitter
      const tangent = rv.clone().sub(n.clone().multiplyScalar(vn));
      if (tangent.lengthSq() > 1e-9){
        tangent.normalize();
        const jt = -rv.dot(tangent) / (invMassA+invMassB);
        const mu = sim.params.friction;
        const clamped = clamp(jt, -mu*jImpulse, mu*jImpulse);
        const impT = tangent.multiplyScalar(clamped);
        a.velocity.addScaledVector(impT, -invMassA);
        b.velocity.addScaledVector(impT,  invMassB);
      }
    }
  }
}
