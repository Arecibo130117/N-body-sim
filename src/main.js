import "./style.css";
import * as THREE from "three";
import Stats from "stats.js";
import { EffectComposer } from "three/examples/jsm/postprocessing/EffectComposer.js";
import { RenderPass } from "three/examples/jsm/postprocessing/RenderPass.js";
import { UnrealBloomPass } from "three/examples/jsm/postprocessing/UnrealBloomPass.js";
import { SMAAPass } from "three/examples/jsm/postprocessing/SMAAPass.js";

import { Simulation } from "./sim/Simulation.js";
import { Body, MaterialType, MaterialPresets } from "./sim/Body.js";
import { BodyRenderer } from "./render/BodyRenderer.js";
import { Trails } from "./render/Trails.js";
import { ParticleSystem } from "./render/Particles.js";
import { UI } from "./ui/UI.js";
import { randRange, seededRand, clamp } from "./utils/math.js";

const app = document.getElementById("app");

// THREE basics
const renderer = new THREE.WebGLRenderer({ antialias: false, powerPreference: "high-performance" });
renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.outputColorSpace = THREE.SRGBColorSpace;
app.appendChild(renderer.domElement);

const scene = new THREE.Scene();
scene.fog = new THREE.FogExp2(0x070a12, 0.000002);

const camera = new THREE.PerspectiveCamera(55, window.innerWidth/window.innerHeight, 0.1, 2e7);
camera.position.set(0, 900, 1600);

const controls = makeControls(camera, renderer.domElement);

// space backdrop
scene.add(makeStarfield());

// lights
scene.add(new THREE.AmbientLight(0xffffff, 0.18));
const key = new THREE.DirectionalLight(0xffffff, 1.05);
key.position.set(0.6, 1.0, 0.2).normalize();
scene.add(key);

// postprocessing
const composer = new EffectComposer(renderer);
composer.addPass(new RenderPass(scene, camera));
const bloom = new UnrealBloomPass(new THREE.Vector2(window.innerWidth, window.innerHeight), 0.85, 0.6, 0.15);
composer.addPass(bloom);
composer.addPass(new SMAAPass(window.innerWidth*renderer.getPixelRatio(), window.innerHeight*renderer.getPixelRatio()));

// sim
const sim = new Simulation();

// render helpers
const bodyRenderer = new BodyRenderer(scene);
const trails = new Trails(scene);
const particles = new ParticleSystem(scene);

// UI
const ui = new UI(sim, {
  seedDefaultSystem,
  seedRandomSystem,
  spawnSimple
});

// debris events
sim.onDebris((payload)=>{
  particles.emit(payload);
});

// stats
const stats = new Stats();
stats.showPanel(0);
stats.dom.style.position = "absolute";
stats.dom.style.left = "12px";
stats.dom.style.bottom = "12px";
stats.dom.style.zIndex = "5";
app.appendChild(stats.dom);

// selection
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();
let shiftDown = false;

window.addEventListener("keydown", (e)=>{ if (e.key === "Shift") shiftDown = true; });
window.addEventListener("keyup", (e)=>{ if (e.key === "Shift") shiftDown = false; });

renderer.domElement.addEventListener("pointerdown", (e)=>{
  if (e.button !== 0) return;
  const rect = renderer.domElement.getBoundingClientRect();
  mouse.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -(((e.clientY - rect.top) / rect.height) * 2 - 1);

  raycaster.setFromCamera(mouse, camera);
  const picked = bodyRenderer.pick(raycaster, sim.bodies);

  if (shiftDown){
    // Spawn at ray-plane intersection (plane through origin)
    const plane = new THREE.Plane(new THREE.Vector3(0,1,0), 0);
    const hit = new THREE.Vector3();
    raycaster.ray.intersectPlane(plane, hit);

    const dir = new THREE.Vector3().copy(raycaster.ray.direction).normalize();
    // spawn slightly in front
    const pos = hit.clone().addScaledVector(dir, 0.0);

    const b = makeBody({
      material: MaterialType.ROCK,
      name: "Dropped",
      mass: ui.spawn.mass,
      position: pos,
      velocity: dir.clone().multiplyScalar(ui.spawn.speed)
    });
    sim.addBody(b);
    return;
  }

  ui.setSelected(picked);
});

// initial system
seedDefaultSystem();

function tick(){
  requestAnimationFrame(tick);
  stats.begin();

  const realDt = Math.min(0.033, clock.getDelta());

  if (!ui.paused){
    sim.step(realDt);
  }

  // flame spawning
  for (const b of sim.bodies){
    if (b.burning || b.temperatureK > 2200){
      particles.emitFlamesAround(b);
    }
  }

  particles.step(realDt*sim.params.timeScale/60, sim.time);

  bodyRenderer.sync(sim.bodies, sim.time);
  trails.sync(sim.bodies);

  controls.update();

  composer.render();
  stats.end();

  ui.updateStats(1000/Math.max(1, stats._frames ? (stats._startTime - stats._prevTime) : 16), sim);
}

const clock = new THREE.Clock();
tick();

window.addEventListener("resize", ()=>{
  renderer.setSize(window.innerWidth, window.innerHeight);
  camera.aspect = window.innerWidth/window.innerHeight;
  camera.updateProjectionMatrix();
  composer.setSize(window.innerWidth, window.innerHeight);
});

// --- seeding & helpers -------------------------------------------------------

function makeBody({ name="Body", material=MaterialType.ROCK, mass=1e12, position, velocity }){
  const preset = MaterialPresets[material];
  const density = preset.density;
  const radius = Math.cbrt((3*mass)/(4*Math.PI*density));
  const b = new Body({
    name,
    mass,
    radius,
    density,
    position,
    velocity,
    material,
    seed: (Math.random()*1e9)>>>0,
    spinAxis: new THREE.Vector3(Math.random()-0.5, Math.random()-0.5, Math.random()-0.5).normalize(),
    spinRate: randRange(-0.4, 0.4)
  });
  // start with a little heat for nicer emissive
  b.heatJ = preset.cp * b.mass * randRange(80, 350);
  b.temperatureK = 3 + b.heatJ/(preset.cp*b.mass);
  return b;
}

function seedDefaultSystem(){
  sim.reset();

  // A "star" at center (not rendered special, just huge mass & hot)
  const star = makeBody({
    name: "Star",
    material: MaterialType.METAL,
    mass: 2.2e16,
    position: new THREE.Vector3(0,0,0),
    velocity: new THREE.Vector3(0,0,0),
  });
  star.radius = 260;
  star.temperatureK = 4200;
  star.heatJ = MaterialPresets[star.material].cp * star.mass * (star.temperatureK - 3);
  sim.addBody(star);

  // Planets
  addPlanet("Terra", MaterialType.ROCK, 7.5e14, 950, 0.0);
  addPlanet("Aqua", MaterialType.ICE,  2.6e14, 1300, 0.7);
  addPlanet("Gasor", MaterialType.GAS,  1.2e15, 1650, 1.2);
  addPlanet("Vulc", MaterialType.ROCK,  1.8e14, 720,  -0.9);

  // some asteroids
  const rng = seededRand(1337);
  for (let i=0;i<80;i++){
    const r = randRange(500, 2000);
    const ang = rng()*Math.PI*2;
    const y = randRange(-60, 60);
    const pos = new THREE.Vector3(Math.cos(ang)*r, y, Math.sin(ang)*r);
    const v = circularVelocity(pos, star.mass).multiplyScalar(randRange(0.85,1.18));
    const m = randRange(2e9, 2e11);
    sim.addBody(makeBody({ name:"Ast", material: MaterialType.ROCK, mass: m, position: pos, velocity: v }));
  }
}

function seedRandomSystem(){
  sim.reset();
  const rng = seededRand((Math.random()*1e9)>>>0);

  const star = makeBody({
    name: "Star",
    material: MaterialType.METAL,
    mass: randRange(1.4e16, 3.5e16),
    position: new THREE.Vector3(0,0,0),
    velocity: new THREE.Vector3(0,0,0),
  });
  star.radius = randRange(220, 360);
  star.temperatureK = randRange(3500, 5200);
  star.heatJ = MaterialPresets[star.material].cp * star.mass * (star.temperatureK - 3);
  sim.addBody(star);

  const nPlanets = Math.floor(randRange(3,6));
  for (let i=0;i<nPlanets;i++){
    const mat = [MaterialType.ROCK, MaterialType.ICE, MaterialType.GAS][Math.floor(rng()*3)];
    const m = randRange(1e14, 2e15) * (mat === MaterialType.GAS ? 1.5 : 1.0);
    const r = randRange(650, 1800);
    const ang = rng()*Math.PI*2;
    const pos = new THREE.Vector3(Math.cos(ang)*r, randRange(-50,50), Math.sin(ang)*r);
    const v = circularVelocity(pos, star.mass).multiplyScalar(randRange(0.92,1.12));
    sim.addBody(makeBody({ name: "P"+(i+1), material: mat, mass: m, position: pos, velocity: v }));
  }

  const nAst = Math.floor(randRange(60, 140));
  for (let i=0;i<nAst;i++){
    const r = randRange(450, 2400);
    const ang = rng()*Math.PI*2;
    const pos = new THREE.Vector3(Math.cos(ang)*r, randRange(-90,90), Math.sin(ang)*r);
    const v = circularVelocity(pos, star.mass).multiplyScalar(randRange(0.75,1.25));
    const m = randRange(2e9, 1.2e11);
    const mat = rng() < 0.2 ? MaterialType.ICE : MaterialType.ROCK;
    sim.addBody(makeBody({ name:"Ast", material: mat, mass: m, position: pos, velocity: v }));
  }
}

function spawnSimple(kind){
  const center = sim.bodies[0]?.position || new THREE.Vector3();
  const star = sim.bodies[0];
  const r = randRange(800, 2000);
  const ang = Math.random()*Math.PI*2;
  const pos = new THREE.Vector3(Math.cos(ang)*r, randRange(-30,30), Math.sin(ang)*r);
  const v = star ? circularVelocity(pos, star.mass).multiplyScalar(randRange(0.95,1.08)) : new THREE.Vector3();
  const mat = kind === "ice" ? MaterialType.ICE : (kind === "gas" ? MaterialType.GAS : MaterialType.ROCK);
  const m = kind === "gas" ? randRange(8e14, 2e15) : randRange(8e13, 8e14);
  sim.addBody(makeBody({ name: kind.toUpperCase(), material: mat, mass: m, position: pos, velocity: v }));
}

function addPlanet(name, material, mass, orbitalR, phase){
  const pos = new THREE.Vector3(Math.cos(phase)*orbitalR, randRange(-30,30), Math.sin(phase)*orbitalR);
  const v = circularVelocity(pos, sim.bodies[0].mass);
  sim.addBody(makeBody({ name, material, mass, position: pos, velocity: v }));
}

function circularVelocity(pos, centralMass){
  const r = pos.length();
  const speed = Math.sqrt(sim.params.G * centralMass / Math.max(1e-6, r));
  // direction perpendicular to radius in XZ plane
  const dir = new THREE.Vector3(-pos.z, 0, pos.x).normalize();
  return dir.multiplyScalar(speed);
}

function makeStarfield(){
  const geo = new THREE.BufferGeometry();
  const count = 2200;
  const pos = new Float32Array(count*3);
  for (let i=0;i<count;i++){
    const r = 5e6 + Math.random()*7e6;
    const u = Math.random()*2-1;
    const t = Math.random()*Math.PI*2;
    const s = Math.sqrt(1-u*u);
    pos[i*3+0] = r*s*Math.cos(t);
    pos[i*3+1] = r*u;
    pos[i*3+2] = r*s*Math.sin(t);
  }
  geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  const mat = new THREE.PointsMaterial({ size: 1200, sizeAttenuation: true, transparent: true, opacity: 0.6 });
  const pts = new THREE.Points(geo, mat);
  pts.frustumCulled = false;
  return pts;
}

// --- minimal orbit controls (no external deps) -------------------------------

function makeControls(camera, dom){
  let isDown = false;
  let lastX=0, lastY=0;
  let yaw = 0.7, pitch = 0.4, dist = camera.position.length();
  const target = new THREE.Vector3(0,0,0);

  dom.addEventListener("pointerdown", (e)=>{
    if (e.button !== 0) return;
    isDown = true;
    lastX = e.clientX; lastY = e.clientY;
  });
  window.addEventListener("pointerup", ()=>{ isDown = false; });
  window.addEventListener("pointermove", (e)=>{
    if (!isDown) return;
    const dx = e.clientX - lastX;
    const dy = e.clientY - lastY;
    lastX = e.clientX; lastY = e.clientY;
    yaw -= dx * 0.0032;
    pitch -= dy * 0.0032;
    pitch = clamp(pitch, -1.2, 1.2);
  });
  dom.addEventListener("wheel", (e)=>{
    dist *= (e.deltaY > 0) ? 1.08 : 0.92;
    dist = clamp(dist, 150, 1.2e6);
  }, { passive: true });

  // target follows selected body slowly
  function update(){
    const sel = ui.selected;
    if (sel){
      target.lerp(sel.position, 0.08);
    } else {
      target.lerp(new THREE.Vector3(0,0,0), 0.02);
    }
    const x = target.x + dist * Math.cos(pitch) * Math.cos(yaw);
    const y = target.y + dist * Math.sin(pitch);
    const z = target.z + dist * Math.cos(pitch) * Math.sin(yaw);
    camera.position.set(x,y,z);
    camera.lookAt(target);
  }

  return { update };
}
