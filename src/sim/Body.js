import * as THREE from "three";

export const MaterialType = Object.freeze({
  ROCK: "rock",
  ICE: "ice",
  METAL: "metal",
  GAS: "gas"
});

// Rough physical-ish parameters for gameplay (not strict astrophysics).
export const MaterialPresets = {
  [MaterialType.ROCK]:  { density: 3500, cp: 900,  emiss: 0.85, igniteK: 2200, ablateK: 1800, meltK: 1400, color: new THREE.Color("#6d6a68") },
  [MaterialType.ICE]:   { density: 1200, cp: 2100, emiss: 0.95, igniteK: 900,  ablateK: 600,  meltK: 273,  color: new THREE.Color("#b6d9ff") },
  [MaterialType.METAL]: { density: 7800, cp: 500,  emiss: 0.35, igniteK: 2500, ablateK: 2200, meltK: 1700, color: new THREE.Color("#8b8f94") },
  [MaterialType.GAS]:   { density: 400,  cp: 1300, emiss: 0.7,  igniteK: 1200, ablateK: 1000, meltK: 1,    color: new THREE.Color("#7fb1ff") }
};

let _nextId = 1;

export class Body {
  constructor({
    name="Body",
    mass=1e12,
    radius=50,
    density=3500,
    position=new THREE.Vector3(),
    velocity=new THREE.Vector3(),
    spinAxis=new THREE.Vector3(0,1,0),
    spinRate=0,
    material=MaterialType.ROCK,
    seed=1
  } = {}){
    this.id = _nextId++;
    this.name = name;
    this.mass = mass;
    this.radius = radius;
    this.density = density;
    this.position = position.clone();
    this.velocity = velocity.clone();
    this.accel = new THREE.Vector3();
    this.prevAccel = new THREE.Vector3();

    // Spin (purely visual + for crater rotation)
    this.spinAxis = spinAxis.clone().normalize();
    this.spinRate = spinRate; // rad/s
    this.spinAngle = 0;

    // Thermal + damage model
    this.material = material;
    this.seed = seed;

    this.temperatureK = 250; // start
    this.heatJ = 0;          // tracked energy store for numerics
    this.burning = false;    // volatile burning (visual)
    this.sublimating = false;

    this.damage = 0;         // 0..1
    this.craters = [];       // capped list of {dir:Vector3, r, depth, age}

    this.alive = true;

    // For trail rendering
    this.trail = [];
  }

  get volume(){ return (4/3) * Math.PI * this.radius**3; }
  get densityFromMass(){ return this.mass / this.volume; }

  setMassFromRadius(){
    this.mass = this.density * this.volume;
  }

  setRadiusFromMass(){
    this.radius = Math.cbrt((3*this.mass)/(4*Math.PI*this.density));
  }

  addCrater(worldNormal, strength){
    // worldNormal points from center to impact point
    const dir = worldNormal.clone().normalize();
    const r = this.radius * (0.06 + 0.14 * Math.min(1, strength));
    const depth = this.radius * (0.01 + 0.06 * Math.min(1, strength));
    const age = 0;
    this.craters.unshift({ dir, r, depth, age });
    if (this.craters.length > 32) this.craters.pop();
    this.damage = Math.min(1, this.damage + 0.12*strength);
  }

  stepThermal(dt){
    const preset = MaterialPresets[this.material] || MaterialPresets.rock;
    // crude radiative cooling to space: dT ~ -k*(T^4 - Tspace^4)
    const T = this.temperatureK;
    const Tspace = 3;
    const sigma = 5.670374419e-8;
    // scale factor to make it visible in "game seconds"
    const k = 1e-7;
    const area = 4*Math.PI*this.radius*this.radius;
    const P = preset.emiss * sigma * area * (Math.pow(T,4) - Math.pow(Tspace,4));
    const dE = -k * P * dt;
    this.heatJ = Math.max(0, this.heatJ + dE);
    this.temperatureK = 3 + this.heatJ / (Math.max(1e-9, preset.cp * this.mass));

    // phase / burning
    this.burning = (this.temperatureK > preset.igniteK && this.material !== MaterialType.METAL);
    this.sublimating = (this.temperatureK > preset.ablateK && this.material === MaterialType.ICE);

    // mass loss (ablation / sublimation) for drama
    if (this.temperatureK > preset.ablateK){
      const over = (this.temperatureK - preset.ablateK);
      const loss = (over / preset.ablateK) * 1e-6 * this.mass * dt; // tuned
      this.mass = Math.max(1e-6, this.mass - loss);
      this.setRadiusFromMass();
      // losing mass consumes heat
      this.heatJ = Math.max(0, this.heatJ - loss * preset.cp * 1500);
    }
  }
}
