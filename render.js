// render.js
// Three.js + postprocessing (Bloom 필수, FXAA) + Instanced bodies + Trails + Particles (budget)

import * as THREE from "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js";
import { OrbitControls } from "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/controls/OrbitControls.js";
import { EffectComposer } from "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/postprocessing/EffectComposer.js";
import { RenderPass } from "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/postprocessing/RenderPass.js";
import { UnrealBloomPass } from "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/postprocessing/UnrealBloomPass.js";
import { ShaderPass } from "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/postprocessing/ShaderPass.js";
import { FXAAShader } from "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/shaders/FXAAShader.js";

function makeStarfield(count = 6000, radius = 2200) {
  const geo = new THREE.BufferGeometry();
  const pos = new Float32Array(count * 3);
  const col = new Float32Array(count * 3);

  for (let i = 0; i < count; i++) {
    // 구면 랜덤 + 약간의 밀도 변형
    const u = Math.random(), v = Math.random();
    const theta = 2 * Math.PI * u;
    const phi = Math.acos(2 * v - 1);
    const r = radius * (0.75 + 0.25 * Math.pow(Math.random(), 0.22));

    const x = r * Math.sin(phi) * Math.cos(theta);
    const y = r * Math.sin(phi) * Math.sin(theta);
    const z = r * Math.cos(phi);

    pos[i * 3] = x; pos[i * 3 + 1] = y; pos[i * 3 + 2] = z;

    const b = 0.55 + 0.45 * Math.random();
    col[i * 3] = b; col[i * 3 + 1] = b; col[i * 3 + 2] = b;
  }

  geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  geo.setAttribute("color", new THREE.BufferAttribute(col, 3));
  const mat = new THREE.PointsMaterial({ size: 1.6, vertexColors: true, transparent: true, opacity: 0.9, depthWrite: false });
  const pts = new THREE.Points(geo, mat);
  return pts;
}

// --- Particle system (Points + custom shader, budget & fade)
function createParticleSystem(maxParticles) {
  const geo = new THREE.BufferGeometry();
  const pos = new Float32Array(maxParticles * 3);
  const vel = new Float32Array(maxParticles * 3);
  const col = new Float32Array(maxParticles * 3);
  const life = new Float32Array(maxParticles);
  const size = new Float32Array(maxParticles);
  const kind = new Float32Array(maxParticles);

  geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  geo.setAttribute("aColor", new THREE.BufferAttribute(col, 3));
  geo.setAttribute("aLife", new THREE.BufferAttribute(life, 1));
  geo.setAttribute("aSize", new THREE.BufferAttribute(size, 1));
  geo.setAttribute("aKind", new THREE.BufferAttribute(kind, 1));

  const mat = new THREE.ShaderMaterial({
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
    uniforms: { uTime: { value: 0 } },
    vertexShader: `
      attribute vec3 aColor;
      attribute float aLife;
      attribute float aSize;
      attribute float aKind;
      varying vec3 vColor;
      varying float vLife;
      varying float vKind;
      void main() {
        vColor = aColor;
        vLife = aLife;
        vKind = aKind;
        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
        gl_PointSize = aSize * (300.0 / -mvPosition.z);
        gl_Position = projectionMatrix * mvPosition;
      }
    `,
    fragmentShader: `
      varying vec3 vColor;
      varying float vLife;
      varying float vKind;
      void main() {
        // soft sprite
        vec2 uv = gl_PointCoord * 2.0 - 1.0;
        float d = dot(uv, uv);
        float soft = exp(-2.6 * d);
        float alpha = clamp(vLife, 0.0, 1.0) * soft;

        // smoke: darker, less additive feel (still additive but dim)
        vec3 c = vColor;
        if (vKind > 1.5) { // smoke
          c *= 0.45;
          alpha *= 0.75;
        }
        gl_FragColor = vec4(c, alpha);
      }
    `
  });

  const points = new THREE.Points(geo, mat);

  let aliveCount = 0;
  const free = [];
  for (let i = 0; i < maxParticles; i++) free.push(i);

  function spawnBurst(req, budgetSoft) {
    const { kind: k, origin, axis, count, energy } = req;
    let want = count | 0;
    if (!Number.isFinite(want) || want <= 0) return;

    // budget 기반으로 clamp + 페이드(너무 많으면 적게 뿌림)
    want = Math.min(want, Math.max(0, budgetSoft - aliveCount));
    if (want <= 0) return;

    const ax = axis.clone().normalize();
    const tmp = new THREE.Vector3();
    for (let s = 0; s < want; s++) {
      if (free.length === 0) break;
      const id = free.pop();
      aliveCount++;

      // 방향: 축 주변으로 분포
      const cone = (k === "flash") ? 0.35 : (k === "fire" ? 0.55 : 0.8);
      const u = Math.random(), v = Math.random();
      const cosA = Math.cos(cone);
      const z = cosA + (1 - cosA) * u;
      const phi = 2 * Math.PI * v;
      const r = Math.sqrt(Math.max(0, 1 - z * z));
      tmp.set(r * Math.cos(phi), r * Math.sin(phi), z);

      // rotate tmp to axis
      const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), ax);
      tmp.applyQuaternion(q).normalize();

      const spdBase = (k === "flash") ? 26 : (k === "fire" ? 16 : 8);
      const spd = spdBase * (0.5 + 1.2 * Math.random()) * Math.min(2.0, 0.5 + Math.sqrt(energy) / 70);

      // position
      pos[id * 3] = origin.x + (Math.random() - 0.5) * 0.8;
      pos[id * 3 + 1] = origin.y + (Math.random() - 0.5) * 0.8;
      pos[id * 3 + 2] = origin.z + (Math.random() - 0.5) * 0.8;

      // velocity
      vel[id * 3] = tmp.x * spd;
      vel[id * 3 + 1] = tmp.y * spd;
      vel[id * 3 + 2] = tmp.z * spd;

      // color
      if (k === "flash") {
        col[id * 3] = 1.0; col[id * 3 + 1] = 0.95; col[id * 3 + 2] = 0.8;
        kind[id] = 0;
        size[id] = 22 + 18 * Math.random();
        life[id] = 0.65 + 0.35 * Math.random();
      } else if (k === "fire") {
        col[id * 3] = 1.0; col[id * 3 + 1] = 0.45 + 0.35 * Math.random(); col[id * 3 + 2] = 0.1;
        kind[id] = 1;
        size[id] = 14 + 14 * Math.random();
        life[id] = 0.8 + 0.5 * Math.random();
      } else { // smoke
        col[id * 3] = 0.35; col[id * 3 + 1] = 0.35; col[id * 3 + 2] = 0.35;
        kind[id] = 2;
        size[id] = 18 + 22 * Math.random();
        life[id] = 1.2 + 1.5 * Math.random();
      }
    }

    geo.attributes.position.needsUpdate = true;
    geo.attributes.aColor.needsUpdate = true;
    geo.attributes.aLife.needsUpdate = true;
    geo.attributes.aSize.needsUpdate = true;
    geo.attributes.aKind.needsUpdate = true;
  }

  function update(dt, gravityY = 0, wind = new THREE.Vector3()) {
    // update attributes on CPU (budget 제어로 수량 제한)
    for (let i = 0; i < maxParticles; i++) {
      if (life[i] <= 0) continue;
      life[i] -= dt;

      // motion
      vel[i * 3] += wind.x * dt;
      vel[i * 3 + 1] += (wind.y - gravityY) * dt;
      vel[i * 3 + 2] += wind.z * dt;

      pos[i * 3] += vel[i * 3] * dt;
      pos[i * 3 + 1] += vel[i * 3 + 1] * dt;
      pos[i * 3 + 2] += vel[i * 3 + 2] * dt;

      // drag
      const drag = (kind[i] > 1.5) ? 1.6 : 0.8;
      vel[i * 3] *= Math.exp(-drag * dt);
      vel[i * 3 + 1] *= Math.exp(-drag * dt);
      vel[i * 3 + 2] *= Math.exp(-drag * dt);

      // smoke rises a bit
      if (kind[i] > 1.5) vel[i * 3 + 1] += 1.6 * dt;

      // fade-out handled in shader via aLife (life used as alpha factor)
      if (life[i] <= 0) {
        life[i] = 0;
        free.push(i);
        aliveCount--;
      }
    }

    geo.attributes.position.needsUpdate = true;
    geo.attributes.aLife.needsUpdate = true;
  }

  return { points, spawnBurst, update, get aliveCount() { return aliveCount; } };
}

// --- Trails
function createTrailSystem(maxTrails, trailLen) {
  const trails = new Map(); // bodyIndex -> {geo,line,positions,colors,head,count}
  const material = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.95 });

  function ensure(scene, bodyIndex) {
    if (trails.has(bodyIndex)) return trails.get(bodyIndex);
    if (trails.size >= maxTrails) return null;

    const positions = new Float32Array(trailLen * 3);
    const colors = new Float32Array(trailLen * 3);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    geo.setDrawRange(0, 0);

    const line = new THREE.Line(geo, material);
    scene.add(line);

    const t = { geo, line, positions, colors, head: 0, count: 0 };
    trails.set(bodyIndex, t);
    return t;
  }

  function remove(scene, bodyIndex) {
    const t = trails.get(bodyIndex);
    if (!t) return;
    scene.remove(t.line);
    t.geo.dispose();
    trails.delete(bodyIndex);
  }

  function setTrailLen(newLen) {
    // 간단 구현: 기존 제거 후 재생성(성능 디그레이드 때만 호출)
    for (const [idx, t] of trails.entries()) {
      t.line.visible = false;
    }
  }

  /**
   * updateTrails()
   * - "끊기지 않고" 부드럽게: ring buffer
   * - 속도에 따라 밝기 변화(빠를수록 더 밝게)
   */
  function updateTrails(scene, bodies, params) {
    const { pos, vel, active, emissive, materialId } = bodies;

    // 어떤 바디를 trail로 둘지: 가장 큰 질량/가까운 것 등을 main에서 선택 가능
    for (const [idx, t] of trails.entries()) {
      if (!active[idx]) {
        // 자연스럽게 페이드: opacity 감소 대신 line 숨김
        t.line.visible = false;
        continue;
      }
      t.line.visible = true;

      const k = idx * 3;
      const x = pos[k], y = pos[k + 1], z = pos[k + 2];
      const vx = vel[k], vy = vel[k + 1], vz = vel[k + 2];
      const spd = Math.sqrt(vx * vx + vy * vy + vz * vz);

      // ring push
      const head = t.head;
      t.positions[head * 3] = x;
      t.positions[head * 3 + 1] = y;
      t.positions[head * 3 + 2] = z;

      // speed -> brightness
      const e = emissive[idx];
      const bright = Math.min(1, 0.12 + 0.018 * spd + 0.75 * e);

      // material tint (RGB)
      const tint = params.materials[materialId[idx]].trailTint;
      t.colors[head * 3] = tint.x * bright;
      t.colors[head * 3 + 1] = tint.y * bright;
      t.colors[head * 3 + 2] = tint.z * bright;

      t.head = (head + 1) % params.trailLength;
      t.count = Math.min(params.trailLength, t.count + 1);

      // write to geometry in correct order (to avoid breaks)
      // cost: O(L) per trail, but trail count is limited by budget
      const L = params.trailLength;
      const outPos = t.geo.attributes.position.array;
      const outCol = t.geo.attributes.color.array;
      const start = (t.head - t.count + L) % L;

      for (let i = 0; i < t.count; i++) {
        const src = (start + i) % L;
        outPos[i * 3] = t.positions[src * 3];
        outPos[i * 3 + 1] = t.positions[src * 3 + 1];
        outPos[i * 3 + 2] = t.positions[src * 3 + 2];

        // tail fade
        const fade = i / Math.max(1, t.count - 1);
        outCol[i * 3] = t.colors[src * 3] * fade;
        outCol[i * 3 + 1] = t.colors[src * 3 + 1] * fade;
        outCol[i * 3 + 2] = t.colors[src * 3 + 2] * fade;
      }

      t.geo.setDrawRange(0, t.count);
      t.geo.attributes.position.needsUpdate = true;
      t.geo.attributes.color.needsUpdate = true;
    }
  }

  return { ensure, remove, updateTrails, setTrailLen, trails };
}

export function createRenderSystem(container, params, bodies) {
  const scene = new THREE.Scene();

  const renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: "high-performance" });
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio));
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = params.exposure;
  container.appendChild(renderer.domElement);

  const camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.01, 80000);
  camera.position.set(0, 110, 260);

  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.06;
  controls.maxDistance = 12000;

  // Lights
  const amb = new THREE.AmbientLight(0xffffff, 0.15);
  scene.add(amb);
  const dir = new THREE.DirectionalLight(0xffffff, 0.45);
  dir.position.set(1, 1, 0.6);
  scene.add(dir);

  // Starfield background
  scene.add(makeStarfield(params.starCount, params.starRadius));

  // Bodies: Instanced spheres
  const sphereGeo = new THREE.IcosahedronGeometry(1, 2);
  const bodyMat = new THREE.MeshStandardMaterial({
    metalness: 0.1,
    roughness: 0.65,
    vertexColors: true,
    emissive: new THREE.Color(0xffffff),
    emissiveIntensity: 1.0
  });

  const inst = new THREE.InstancedMesh(sphereGeo, bodyMat, bodies.capacity);
  inst.instanceMatrix.setUsage(THREE.DynamicDrawUsage);

  const instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(bodies.capacity * 3), 3);
  inst.instanceColor = instanceColor;
  scene.add(inst);

  // Spawn preview arrow
  const arrow = new THREE.ArrowHelper(new THREE.Vector3(1,0,0), new THREE.Vector3(0,0,0), 10, 0x88ccff);
  arrow.visible = false;
  scene.add(arrow);

  // Trails
  const trailSys = createTrailSystem(params.maxTrails, params.trailLength);

  // Particles
  const particleSys = createParticleSystem(params.particleBudgetHard);
  scene.add(particleSys.points);

  // Postprocessing
  const composer = new EffectComposer(renderer);
  composer.addPass(new RenderPass(scene, camera));

  const bloom = new UnrealBloomPass(
    new THREE.Vector2(window.innerWidth, window.innerHeight),
    params.bloomStrength,
    params.bloomRadius,
    params.bloomThreshold
  );
  composer.addPass(bloom);

  const fxaa = new ShaderPass(FXAAShader);
  fxaa.material.uniforms["resolution"].value.set(1 / window.innerWidth, 1 / window.innerHeight);
  composer.addPass(fxaa);

  function updateBloomFromParams() {
    bloom.strength = params.bloomStrength;
    bloom.radius = params.bloomRadius;
    bloom.threshold = params.bloomThreshold;
  }

  function updateExposure() {
    renderer.toneMappingExposure = params.exposure;
  }

  function resize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
    composer.setSize(window.innerWidth, window.innerHeight);
    fxaa.material.uniforms["resolution"].value.set(1 / window.innerWidth, 1 / window.innerHeight);
  }
  window.addEventListener("resize", resize);

  /**
   * updateBodyInstances()
   * - LOD-1: 찌그러짐/손상도 기반 scale 변형
   * - emissive는 thermal emissiveLevel에 연동
   * - material별 baseColor 차등
   */
  const tmpMat = new THREE.Matrix4();
  const tmpQ = new THREE.Quaternion();
  const tmpS = new THREE.Vector3();
  const tmpP = new THREE.Vector3();
  const tmpE = new THREE.Euler();

  function updateBodyInstances() {
    const { pos, radius, active, materialId, emissive, damage, omega, meltFrac, vaporFrac } = bodies;
    let drawCount = 0;

    for (let i = 0; i < bodies.capacity; i++) {
      if (!active[i]) continue;

      const k = i * 3;
      tmpP.set(pos[k], pos[k + 1], pos[k + 2]);

      // Spin visualization: 간이로 Euler 증가분 반영(실제 rigid orientation은 별도 추적이 이상적이지만 여기선 근사)
      const ox = omega[k], oy = omega[k + 1], oz = omega[k + 2];
      tmpE.set(ox * 0.02, oy * 0.02, oz * 0.02);
      tmpQ.setFromEuler(tmpE);

      // LOD-1 deformation: damage/melt에 따라 타원체 스케일
      const d = Math.min(1, damage[i]);
      const melt = meltFrac[i];
      const vap = vaporFrac[i];

      const baseR = radius[i];
      const squash = 1 - 0.22 * d - 0.18 * melt;
      const stretch = 1 + 0.28 * d + 0.12 * (Math.random() - 0.5);

      // 고온 가스는 더 "퍼짐"
      const puff = 1 + 0.9 * vap;

      tmpS.set(baseR * stretch * puff, baseR * squash * puff, baseR * (0.92 + 0.18 * d) * puff);

      tmpMat.compose(tmpP, tmpQ, tmpS);
      inst.setMatrixAt(drawCount, tmpMat);

      // color: material base + emissive level
      const mat = params.materials[materialId[i]];
      const e = emissive[i];

      // "HDR-like" feel: emissive가 높을수록 더 밝게 (Bloom이 먹게)
      const c = mat.baseColor.clone().lerp(mat.glowColor, Math.min(1, 0.65 * e + 0.25 * vap));
      instanceColor.setXYZ(drawCount, c.x, c.y, c.z);

      drawCount++;
    }

    inst.count = drawCount;
    inst.instanceMatrix.needsUpdate = true;
    inst.instanceColor.needsUpdate = true;
  }

  function render(dt) {
    controls.update();
    particleSys.update(dt, 0.0, params.wind);

    updateBodyInstances();
    trailSys.updateTrails(scene, bodies, params);
    updateBloomFromParams();
    updateExposure();

    composer.render();
  }

  return {
    scene, camera, renderer, composer, controls,
    arrow,
    trailSys,
    particleSys,
    render,
    resize,
    updateBloomFromParams,
    updateExposure
  };
}
