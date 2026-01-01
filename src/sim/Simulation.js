import * as THREE from "three";
import { Body, MaterialPresets, MaterialType } from "./Body.js";
import { BarnesHutOctree } from "./barnesHut.js";
import { handleCollisions } from "./collisions.js";
import { clamp, randRange, seededRand } from "../utils/math.js";

export class Simulation {
  constructor(){
    this.bodies = [];
    this.time = 0;

    this.params = {
      // Physics scaling
      G: 6.674e-3,          // game-scaled "gravity constant"
      softening: 2.0,       // prevents singularities
      dt: 1/60,
      timeScale: 600,       // seconds simulated per real second
      integrator: "leapfrog",
      useBarnesHut: true,
      barnesTheta: 0.6,

      // Collisions & materials
      restitution: 0.15,
      friction: 0.10,
      mergeBias: 0.00,
      fragBias: 0.00,
      debrisGravityCutoff: 1e7, // bodies lighter than this don't attract each other (perf)

      // Visual-driven physics toggles
      enableThermal: true,
      enableReaccretion: true,
      reaccreteRadiusFactor: 1.8,
      reaccreteSpeedFactor: 0.55,

      // Limits
      maxBodies: 800,
      maxDebris: 5000,
    };

    this._tree = new BarnesHutOctree({ center: new THREE.Vector3(), halfSize: 2e6, theta: this.params.barnesTheta });
    this._tmpA = new THREE.Vector3();
    this._tmpB = new THREE.Vector3();

    // Debris list (massless visual particles) are handled by renderer; we only emit events.
    this.events = { onDebris: [] };
  }

  onDebris(fn){ this.events.onDebris.push(fn); }

  emitDebris(payload){
    for (const fn of this.events.onDebris) fn(payload);
  }

  reset(){
    this.bodies.length = 0;
    this.time = 0;
  }

  addBody(body){
    if (this.bodies.length >= this.params.maxBodies) return;
    body.prevAccel.set(0,0,0);
    body.accel.set(0,0,0);
    this.bodies.push(body);
  }

  removeBody(body){
    body.alive = false;
  }

  step(realDt){
    const dt = this.params.dt;
    const steps = Math.max(1, Math.min(12, Math.floor((realDt * this.params.timeScale) / dt)));
    const subDt = (realDt * this.params.timeScale) / steps;

    for (let s=0;s<steps;s++){
      this._stepOnce(subDt);
    }
  }

  _stepOnce(dt){
    this.time += dt;

    // Rebuild tree if enabled
    if (this.params.useBarnesHut){
      this._tree.theta = this.params.barnesTheta;
      this._tree.rebuild(this.bodies);
    }

    // Compute accelerations
    for (const b of this.bodies){
      if (!b.alive) continue;
      b.prevAccel.copy(b.accel);
      b.accel.set(0,0,0);
      this._accumulateGravity(b, b.accel);
    }

    // Integrate positions/velocities
    if (this.params.integrator === "leapfrog"){
      // Kick-drift-kick (velocity Verlet-like)
      for (const b of this.bodies){
        if (!b.alive) continue;
        b.velocity.addScaledVector(b.prevAccel, 0.5*dt);
        b.position.addScaledVector(b.velocity, dt);
      }
      // accel at new positions
      if (this.params.useBarnesHut){
        this._tree.rebuild(this.bodies);
      }
      for (const b of this.bodies){
        if (!b.alive) continue;
        b.accel.set(0,0,0);
        this._accumulateGravity(b, b.accel);
      }
      for (const b of this.bodies){
        if (!b.alive) continue;
        b.velocity.addScaledVector(b.accel, 0.5*dt);
      }
    } else {
      // semi-implicit Euler
      for (const b of this.bodies){
        if (!b.alive) continue;
        b.velocity.addScaledVector(b.accel, dt);
        b.position.addScaledVector(b.velocity, dt);
      }
    }

    // Spin (visual)
    for (const b of this.bodies){
      if (!b.alive) continue;
      b.spinAngle = (b.spinAngle + b.spinRate*dt) % (Math.PI*2);
      // Trail history
      b.trail.push(b.position.clone());
      if (b.trail.length > 256) b.trail.shift();
      // crater aging
      for (const c of b.craters) c.age += dt;
    }

    // Collisions
    handleCollisions(this);

    // Thermal
    if (this.params.enableThermal){
      for (const b of this.bodies){
        if (!b.alive) continue;
        b.stepThermal(dt);
      }
    }

    // Reaccretion / clustering for small fragments
    if (this.params.enableReaccretion){
      this._reaccrete(dt);
    }

    // Cleanup dead bodies (but keep stable indices for selection in UI? We'll compact.)
    this.bodies = this.bodies.filter(b => b.alive);
  }

  _accumulateGravity(body, outAccel){
    const G = this.params.G;
    const eps = this.params.softening;
    const p = body.position;

    if (this.params.useBarnesHut){
      this._tree.accumulateAccel(body, G, eps, outAccel);
      return;
    }

    // naive
    for (const other of this.bodies){
      if (!other.alive || other === body) continue;
      // performance: small debris don't attract each other
      if (body.mass < this.params.debrisGravityCutoff && other.mass < this.params.debrisGravityCutoff) continue;

      const dx = other.position.x - p.x;
      const dy = other.position.y - p.y;
      const dz = other.position.z - p.z;
      const r2 = dx*dx + dy*dy + dz*dz + eps*eps;
      const invR = 1/Math.sqrt(r2);
      const invR3 = invR*invR*invR;
      const a = G * other.mass * invR3;
      outAccel.x += dx * a;
      outAccel.y += dy * a;
      outAccel.z += dz * a;
    }
  }

  mergeBodies(a, b){
    if (!a.alive || !b.alive) return;

    const totalMass = a.mass + b.mass;
    const v = a.velocity.clone().multiplyScalar(a.mass).addScaledVector(b.velocity, b.mass).multiplyScalar(1/totalMass);
    const p = a.position.clone().multiplyScalar(a.mass).addScaledVector(b.position, b.mass).multiplyScalar(1/totalMass);

    // Composition heuristic: dominant material wins, or blend rock/ice -> rock
    const mat = (a.mass >= b.mass) ? a.material : b.material;

    // Density & radius from mass-weighted density
    const dens = (a.density*a.mass + b.density*b.mass) / totalMass;

    const merged = new THREE.Vector3();
    const name = a.name + "+" + b.name;

    const m = totalMass;
    const radius = Math.cbrt((3*m)/(4*Math.PI*dens));

    const newBody = new (a.constructor)({
      name,
      mass: m,
      radius,
      density: dens,
      position: p,
      velocity: v,
      material: mat,
      seed: (a.seed ^ (b.seed<<1)) >>> 0,
      spinAxis: a.spinAxis.clone().add(b.spinAxis).normalize(),
      spinRate: (a.spinRate*a.radius + b.spinRate*b.radius)/(a.radius+b.radius+1e-6)
    });

    // Thermal merge
    newBody.heatJ = a.heatJ + b.heatJ;
    newBody.temperatureK = 3 + newBody.heatJ / (MaterialPresets[newBody.material].cp * newBody.mass);

    // Carry some craters
    newBody.craters = [...a.craters.slice(0,12), ...b.craters.slice(0,12)].slice(0, 24);
    newBody.damage = clamp(0.5*(a.damage+b.damage), 0, 1);

    a.alive = false;
    b.alive = false;

    this.addBody(newBody);

    // Emit minor debris for visual
    this.emitDebris({
      kind: "mergeDust",
      origin: p,
      count: 280,
      spread: radius*0.8,
      baseVelocity: v.clone(),
      energy: 0.3
    });
  }

  fragmentBodies(a, b, normal, vRel){
    if (!a.alive || !b.alive) return;

    const totalMass = a.mass + b.mass;
    const p = a.position.clone().multiplyScalar(a.mass).addScaledVector(b.position, b.mass).multiplyScalar(1/totalMass);
    const vcm = a.velocity.clone().multiplyScalar(a.mass).addScaledVector(b.velocity, b.mass).multiplyScalar(1/totalMass);

    // Convert a portion of mass into fragments; keep a few large remnants if possible.
    const catastrophicFactor = clamp((vRel / 40), 0.3, 1.0); // tuned
    const massToFragments = totalMass * (0.55 + 0.35*catastrophicFactor);
    const remnantMass = Math.max(0, totalMass - massToFragments);

    // Determine number of fragments
    const n = Math.floor(10 + 40*catastrophicFactor);
    const rng = seededRand((a.seed*2654435761) ^ (b.seed*1597334677));

    const frags = [];

    // Optional remnant
    if (remnantMass > 0.15*totalMass){
      const dens = (a.density*a.mass + b.density*b.mass) / totalMass;
      const radius = Math.cbrt((3*remnantMass)/(4*Math.PI*dens));
      const rem = new Body({
        name: "Remnant",
        mass: remnantMass,
        radius,
        density: dens,
        position: p.clone().addScaledVector(normal, -radius*0.15),
        velocity: vcm.clone().addScaledVector(normal, -0.15*vRel),
        material: (a.mass>=b.mass? a.material : b.material),
        seed: (a.seed ^ b.seed) >>> 0,
        spinAxis: normal.clone().cross(new THREE.Vector3(0,1,0)).normalize(),
        spinRate: randRange(-0.3, 0.3)
      });
      rem.heatJ = (a.heatJ + b.heatJ) * 0.6;
      frags.push(rem);
    }

    // Generate fragment masses (power-law)
    let remaining = massToFragments;
    const masses = [];
    for (let i=0;i<n;i++){
      const x = rng(); // 0..1
      const m = Math.pow(1-x, 6) * (massToFragments*0.18) + massToFragments*0.001;
      masses.push(m);
    }
    const sumM = masses.reduce((s,m)=>s+m,0);
    for (let i=0;i<masses.length;i++) masses[i] = masses[i] / sumM * massToFragments;

    const baseEject = clamp(vRel, 5, 80);
    const dens = (a.density*a.mass + b.density*b.mass) / totalMass;

    for (let i=0;i<n;i++){
      if (remaining <= 0) break;
      const m = masses[i];
      remaining -= m;

      const radius = Math.cbrt((3*m)/(4*Math.PI*dens));
      const dir = randomUnit(rng).addScaledVector(normal, 0.55).normalize();
      const spd = baseEject * (0.25 + 0.95*rng());
      const pos = p.clone().addScaledVector(dir, (a.radius+b.radius) * (0.55 + 0.8*rng()));
      const vel = vcm.clone().addScaledVector(dir, spd);

      const frag = new Body({
        name: "Fragment",
        mass: m,
        radius,
        density: dens,
        position: pos,
        velocity: vel,
        material: chooseMaterialBlend(a.material, b.material, rng),
        seed: (a.seed ^ (i*2654435761)) >>> 0,
        spinAxis: randomUnit(rng),
        spinRate: randRange(-1.8, 1.8)
      });

      // Thermal: fragments get hot
      frag.heatJ = (a.heatJ + b.heatJ) * (0.3/n) + 0.20 * 0.5 * m * spd*spd;
      frags.push(frag);
    }

    a.alive = false;
    b.alive = false;

    // Emit lots of visual debris (dust / sparks)
    this.emitDebris({
      kind: "catastrophic",
      origin: p,
      count: Math.floor(1200 + 2800*catastrophicFactor),
      spread: (a.radius+b.radius)*1.2,
      baseVelocity: vcm.clone(),
      normal: normal.clone(),
      energy: catastrophicFactor
    });

    for (const f of frags){
      this.addBody(f);
    }
  }

  _reaccrete(dt){
    // Very coarse reaccretion:
    // if two bodies are small-ish, close, and relative speed is low -> merge.
    const bodies = this.bodies;
    const G = this.params.G;

    for (let i=0;i<bodies.length;i++){
      const a = bodies[i];
      if (!a.alive) continue;
      if (a.mass > 1e11) continue; // only small fragments
      for (let j=i+1;j<bodies.length;j++){
        const b = bodies[j];
        if (!b.alive) continue;
        if (b.mass > 1e11) continue;

        const dx = b.position.x - a.position.x;
        const dy = b.position.y - a.position.y;
        const dz = b.position.z - a.position.z;
        const dist = Math.sqrt(dx*dx + dy*dy + dz*dz) + 1e-9;

        const r = (a.radius + b.radius) * this.params.reaccreteRadiusFactor;
        if (dist > r) continue;

        const rv = b.velocity.clone().sub(a.velocity);
        const vRel = rv.length();

        const vEsc = Math.sqrt(2*G*(a.mass+b.mass)/Math.max(1e-9, (a.radius+b.radius)));
        if (vRel < vEsc * this.params.reaccreteSpeedFactor){
          this.mergeBodies(a, b);
          break;
        }
      }
    }
  }
}

function randomUnit(rng=Math.random){
  const u = rng()*2-1;
  const t = rng()*Math.PI*2;
  const r = Math.sqrt(1-u*u);
  return new THREE.Vector3(r*Math.cos(t), u, r*Math.sin(t));
}

function chooseMaterialBlend(a, b, rng){
  if (a === b) return a;
  // rock+ice -> rock or ice depending on rng
  if ((a === MaterialType.ICE && b === MaterialType.ROCK) || (a === MaterialType.ROCK && b === MaterialType.ICE)){
    return rng() < 0.55 ? MaterialType.ROCK : MaterialType.ICE;
  }
  // metal tends to stay metal
  if (a === MaterialType.METAL || b === MaterialType.METAL){
    return rng() < 0.7 ? MaterialType.METAL : (a === MaterialType.METAL ? b : a);
  }
  return rng() < 0.5 ? a : b;
}
