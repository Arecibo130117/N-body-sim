// render.js
// - Bloom 필수 + FXAA
// - InstancedMesh: per-instance (color / emissive / roughness / metalness / seed)
// - Procedural rock surface (fbm noise + crater-ish ridges) via onBeforeCompile
// - Trails: teleport/jump 감지해서 reset → "스파이크" 방지, merge/reset 지원

import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { EffectComposer } from "three/addons/postprocessing/EffectComposer.js";
import { RenderPass } from "three/addons/postprocessing/RenderPass.js";
import { UnrealBloomPass } from "three/addons/postprocessing/UnrealBloomPass.js";
import { ShaderPass } from "three/addons/postprocessing/ShaderPass.js";
import { FXAAShader } from "three/addons/shaders/FXAAShader.js";

function makeStarfield(count = 6000, radius = 2200) {
  const geo = new THREE.BufferGeometry();
  const pos = new Float32Array(count * 3);
  const col = new Float32Array(count * 3);

  for (let i = 0; i < count; i++) {
    const u = Math.random(), v = Math.random();
    const theta = 2 * Math.PI * u;
    const phi = Math.acos(2 * v - 1);
    const r = radius * (0.75 + 0.25 * Math.pow(Math.random(), 0.22));

    pos[i * 3]     = r * Math.sin(phi) * Math.cos(theta);
    pos[i * 3 + 1] = r * Math.sin(phi) * Math.sin(theta);
    pos[i * 3 + 2] = r * Math.cos(phi);

    const b = 0.55 + 0.45 * Math.random();
    col[i * 3] = b; col[i * 3 + 1] = b; col[i * 3 + 2] = b;
  }

  geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  geo.setAttribute("color", new THREE.BufferAttribute(col, 3));
  const mat = new THREE.PointsMaterial({
    size: 1.6, vertexColors: true, transparent: true, opacity: 0.9, depthWrite: false
  });
  return new THREE.Points(geo, mat);
}

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
        vec2 uv = gl_PointCoord * 2.0 - 1.0;
        float d = dot(uv, uv);
        float soft = exp(-2.6 * d);
        float alpha = clamp(vLife, 0.0, 1.0) * soft;

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

    want = Math.min(want, Math.max(0, budgetSoft - aliveCount));
    if (want <= 0) return;

    const ax = axis.clone().normalize();
    const tmp = new THREE.Vector3();
    for (let s = 0; s < want; s++) {
      if (free.length === 0) break;
      const id = free.pop();
      aliveCount++;

      const cone = (k === "flash") ? 0.35 : (k === "fire" ? 0.55 : 0.8);
      const u = Math.random(), v = Math.random();
      const cosA = Math.cos(cone);
      const z = cosA + (1 - cosA) * u;
      const phi = 2 * Math.PI * v;
      const rr = Math.sqrt(Math.max(0, 1 - z * z));
      tmp.set(rr * Math.cos(phi), rr * Math.sin(phi), z);

      const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), ax);
      tmp.applyQuaternion(q).normalize();

      const spdBase = (k === "flash") ? 26 : (k === "fire" ? 16 : 8);
      const spd = spdBase * (0.5 + 1.2 * Math.random()) * Math.min(2.0, 0.5 + Math.sqrt(energy) / 70);

      pos[id * 3]     = origin.x + (Math.random() - 0.5) * 0.8;
      pos[id * 3 + 1] = origin.y + (Math.random() - 0.5) * 0.8;
      pos[id * 3 + 2] = origin.z + (Math.random() - 0.5) * 0.8;

      vel[id * 3]     = tmp.x * spd;
      vel[id * 3 + 1] = tmp.y * spd;
      vel[id * 3 + 2] = tmp.z * spd;

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
      } else {
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
    for (let i = 0; i < maxParticles; i++) {
      if (life[i] <= 0) continue;

      life[i] -= dt;

      vel[i * 3]     += wind.x * dt;
      vel[i * 3 + 1] += (wind.y - gravityY) * dt;
      vel[i * 3 + 2] += wind.z * dt;

      pos[i * 3]     += vel[i * 3] * dt;
      pos[i * 3 + 1] += vel[i * 3 + 1] * dt;
      pos[i * 3 + 2] += vel[i * 3 + 2] * dt;

      const drag = (kind[i] > 1.5) ? 1.6 : 0.8;
      const damp = Math.exp(-drag * dt);
      vel[i * 3] *= damp;
      vel[i * 3 + 1] *= damp;
      vel[i * 3 + 2] *= damp;

      if (kind[i] > 1.5) vel[i * 3 + 1] += 1.6 * dt;

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

/**
 * Trails: 스파이크 방지 핵심
 * - body가 순간이동(merge/큰 보정/재할당)하면 자동 reset
 * - resetTrail(idx, x,y,z) 외부에서 호출 가능
 */
function createTrailSystem(params) {
  const trails = new Map();

  function ensure(scene, bodyIndex) {
    let t = trails.get(bodyIndex);
    if (t) return t;
    if (trails.size >= params.maxTrails) return null;

    const L = params.trailLength;
    const positions = new Float32Array(L * 3);
    const colors = new Float32Array(L * 3);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    geo.setDrawRange(0, 0);

    // ✅ trail마다 material을 따로 가져서 "페이드" 가능
    const mat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.95 });
    const line = new THREE.Line(geo, mat);
    scene.add(line);

    t = {
      geo, line, mat,
      positions, colors,
      head: 0, count: 0, L,
      lastX: 0, lastY: 0, lastZ: 0,
      fade: 1.0
    };
    trails.set(bodyIndex, t);
    return t;
  }

  function remove(scene, bodyIndex) {
    const t = trails.get(bodyIndex);
    if (!t) return;
    scene.remove(t.line);
    t.geo.dispose();
    t.mat.dispose();
    trails.delete(bodyIndex);
  }

  function resetTrail(scene, bodyIndex, x, y, z) {
    const t = ensure(scene, bodyIndex);
    if (!t) return;
    // 버퍼를 현재 좌표로 채워 “먼 선분” 자체를 없앰
    for (let i = 0; i < t.L; i++) {
      t.positions[i * 3] = x;
      t.positions[i * 3 + 1] = y;
      t.positions[i * 3 + 2] = z;

      t.colors[i * 3] = 0;
      t.colors[i * 3 + 1] = 0;
      t.colors[i * 3 + 2] = 0;
    }
    t.head = 0;
    t.count = 0;
    t.lastX = x; t.lastY = y; t.lastZ = z;
    t.fade = 1.0;
    t.geo.setDrawRange(0, 0);
    t.geo.attributes.position.needsUpdate = true;
    t.geo.attributes.color.needsUpdate = true;
  }

  function updateTrails(scene, bodies) {
    const { pos, vel, active, emissive, materialId, radius } = bodies;

    for (const [idx, t] of trails.entries()) {
      // trail length 변경 시 재생성
      if (t.L !== params.trailLength) {
        remove(scene, idx);
        continue;
      }

      if (!active[idx]) {
        // ✅ 갑자기 툭 사라지지 않게 페이드
        t.fade *= Math.exp(-2.2 * params._dtRealForFade);
        t.mat.opacity = 0.95 * t.fade;
        if (t.fade < 0.03) remove(scene, idx);
        continue;
      }

      t.fade = 1.0;
      t.mat.opacity = 0.95;

      const k = idx * 3;
      const x = pos[k], y = pos[k + 1], z = pos[k + 2];

      // ✅ “텔레포트/점프” 감지 → 스파이크 방지 핵심
      const dx = x - t.lastX, dy = y - t.lastY, dz = z - t.lastZ;
      const jump = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const maxJump = params.trailTeleportFactor * Math.max(1.0, radius[idx]);
      if (!Number.isFinite(jump) || jump > maxJump) {
        resetTrail(scene, idx, x, y, z);
      }

      t.lastX = x; t.lastY = y; t.lastZ = z;

      const vx = vel[k], vy = vel[k + 1], vz = vel[k + 2];
      const spd = Math.sqrt(vx * vx + vy * vy + vz * vz);

      const head = t.head;
      t.positions[head * 3] = x;
      t.positions[head * 3 + 1] = y;
      t.positions[head * 3 + 2] = z;

      const e = emissive[idx];
      const bright = Math.min(1, 0.08 + 0.016 * spd + 0.85 * e);

      const tint = params.materials[materialId[idx]].trailTint;
      t.colors[head * 3] = tint.x * bright;
      t.colors[head * 3 + 1] = tint.y * bright;
      t.colors[head * 3 + 2] = tint.z * bright;

      t.head = (head + 1) % params.trailLength;
      t.count = Math.min(params.trailLength, t.count + 1);

      const L = params.trailLength;
      const outPos = t.geo.attributes.position.array;
      const outCol = t.geo.attributes.color.array;
      const start = (t.head - t.count + L) % L;

      for (let i = 0; i < t.count; i++) {
        const src = (start + i) % L;
        outPos[i * 3] = t.positions[src * 3];
        outPos[i * 3 + 1] = t.positions[src * 3 + 1];
        outPos[i * 3 + 2] = t.positions[src * 3 + 2];

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

  return { trails, ensure, remove, resetTrail, updateTrails };
}

/**
 * InstancedMesh material에 "암석 행성" 디테일을 넣기 위한 onBeforeCompile
 * - fbm noise로 albedo/roughness 변화 + vertex displacement(미세 크레이터 느낌)
 * - per-instance emissive/rough/metal/seed 지원
 */
function injectPlanetSurfaceShader(material, params) {
  material.onBeforeCompile = (shader) => {
    shader.uniforms.uDisp = { value: params.surfaceDisplacement };
    shader.uniforms.uNoiseScale = { value: params.surfaceNoiseScale };
    shader.uniforms.uEmissiveBoost = { value: params.emissiveBoost };

    shader.vertexShader = shader.vertexShader
      .replace(
        "#include <common>",
        `#include <common>
         attribute float aSeed;
         attribute float aEmissive;
         attribute float aRough;
         attribute float aMetal;
         varying float vSeed;
         varying float vEmissive;
         varying float vRough;
         varying float vMetal;
         varying float vN;

         // value noise (짧고 충분히 “행성 표면” 설득력)
         float hash11(float p){ p=fract(p*0.1031); p*=p+33.33; p*=p+p; return fract(p); }
         float hash31(vec3 p){
           p = fract(p*0.1031);
           p += dot(p, p.yzx + 33.33);
           return fract((p.x+p.y)*p.z);
         }
         float noise3(vec3 p){
           vec3 i=floor(p), f=fract(p);
           f=f*f*(3.0-2.0*f);
           float n000=hash31(i+vec3(0,0,0));
           float n100=hash31(i+vec3(1,0,0));
           float n010=hash31(i+vec3(0,1,0));
           float n110=hash31(i+vec3(1,1,0));
           float n001=hash31(i+vec3(0,0,1));
           float n101=hash31(i+vec3(1,0,1));
           float n011=hash31(i+vec3(0,1,1));
           float n111=hash31(i+vec3(1,1,1));
           float nx00=mix(n000,n100,f.x);
           float nx10=mix(n010,n110,f.x);
           float nx01=mix(n001,n101,f.x);
           float nx11=mix(n011,n111,f.x);
           float nxy0=mix(nx00,nx10,f.y);
           float nxy1=mix(nx01,nx11,f.y);
           return mix(nxy0,nxy1,f.z);
         }
         float fbm(vec3 p){
           float a=0.5; float s=0.0;
           for(int i=0;i<5;i++){
             s += a*noise3(p);
             p *= 2.02;
             a *= 0.5;
           }
           return s;
         }`
      )
      .replace(
        "#include <begin_vertex>",
        `#include <begin_vertex>
         vSeed = aSeed;
         vEmissive = aEmissive;
         vRough = aRough;
         vMetal = aMetal;

         // 암석 표면: 기본 fbm + ridge(크레이터 테두리 느낌)
         vec3 p = position * uNoiseScale + vec3(aSeed*17.3, aSeed*7.1, aSeed*11.7);
         float n = fbm(p);
         float ridge = 1.0 - abs(2.0*n - 1.0);
         float crater = pow(ridge, 3.0);

         vN = clamp(0.55*n + 0.45*crater, 0.0, 1.0);

         // 미세 변위(행성 디테일)
         float disp = (vN - 0.5) * 2.0 * uDisp;
         transformed += normal * disp;`
      );

    shader.fragmentShader = shader.fragmentShader
      .replace(
        "#include <common>",
        `#include <common>
         varying float vSeed;
         varying float vEmissive;
         varying float vRough;
         varying float vMetal;
         varying float vN;
         uniform float uEmissiveBoost;`
      )
      // baseColor 변조: 더 이상 “별 점”처럼 흰색으로 뭉개지지 않게, 표면 톤 변화
      .replace(
        "#include <color_fragment>",
        `#include <color_fragment>
         // vN 기반 표면 톤 변화(암석/빙질 모두 디테일)
         diffuseColor.rgb *= (0.82 + 0.30 * vN);`
      )
      // roughness/metalness per-instance 적용
      .replace(
        "float roughnessFactor = roughness;",
        "float roughnessFactor = clamp(roughness * vRough * (0.90 + 0.25 * vN), 0.02, 1.0);"
      )
      .replace(
        "float metalnessFactor = metalness;",
        "float metalnessFactor = clamp(metalness * vMetal, 0.0, 1.0);"
      )
      // emissive per-instance 추가 → “고온일수록 발광 + bloom”
      .replace(
        "#include <emissivemap_fragment>",
        `#include <emissivemap_fragment>
         totalEmissiveRadiance += diffuseColor.rgb * vEmissive * uEmissiveBoost;`
      );

    material.userData._shader = shader;
  };
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

  scene.add(new THREE.AmbientLight(0xffffff, 0.18));
  const dir = new THREE.DirectionalLight(0xffffff, 0.75);
  dir.position.set(1, 1, 0.6);
  scene.add(dir);

  scene.add(makeStarfield(params.starCount, params.starRadius));

  // ✅ 행성 디테일이 보이려면 기본 지오메트리도 좀 촘촘해야 함
  const sphereGeo = new THREE.IcosahedronGeometry(1, 5);

  const bodyMat = new THREE.MeshStandardMaterial({
    metalness: 0.08,
    roughness: 0.75,
    vertexColors: true,
    emissive: new THREE.Color(0x000000),   // ✅ 기본 emissive 제거 (별처럼 뭉개지는 원인)
    emissiveIntensity: 1.0
  });

  // ✅ “암석 행성” 디테일 셰이더 주입
  injectPlanetSurfaceShader(bodyMat, params);

  const inst = new THREE.InstancedMesh(sphereGeo, bodyMat, bodies.capacity);
  inst.instanceMatrix.setUsage(THREE.DynamicDrawUsage);

  // per-instance attributes
  const instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(bodies.capacity * 3), 3);
  inst.instanceColor = instanceColor;

  const aSeed = new THREE.InstancedBufferAttribute(new Float32Array(bodies.capacity), 1);
  const aEmissive = new THREE.InstancedBufferAttribute(new Float32Array(bodies.capacity), 1);
  const aRough = new THREE.InstancedBufferAttribute(new Float32Array(bodies.capacity), 1);
  const aMetal = new THREE.InstancedBufferAttribute(new Float32Array(bodies.capacity), 1);

  inst.geometry.setAttribute("aSeed", aSeed);
  inst.geometry.setAttribute("aEmissive", aEmissive);
  inst.geometry.setAttribute("aRough", aRough);
  inst.geometry.setAttribute("aMetal", aMetal);

  scene.add(inst);

  const arrow = new THREE.ArrowHelper(new THREE.Vector3(1,0,0), new THREE.Vector3(0,0,0), 10, 0x88ccff);
  arrow.visible = false;
  scene.add(arrow);

  const trailSys = createTrailSystem(params);
  const particleSys = createParticleSystem(params.particleBudgetHard);
  scene.add(particleSys.points);

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

  function resize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
    composer.setSize(window.innerWidth, window.innerHeight);
    fxaa.material.uniforms["resolution"].value.set(1 / window.innerWidth, 1 / window.innerHeight);
  }
  window.addEventListener("resize", resize);

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

      // 스핀(간이) → 표면 디테일이 돌면서 행성처럼 보임
      const ox = omega[k], oy = omega[k + 1], oz = omega[k + 2];
      tmpE.set(ox * 0.02, oy * 0.02, oz * 0.02);
      tmpQ.setFromEuler(tmpE);

      const d = Math.min(1, damage[i]);
      const melt = meltFrac[i];
      const vap = vaporFrac[i];

      const baseR = radius[i];

      // 찌그러짐/박리/용융 느낌을 더 “행성”에서 보이게
      const squash = 1 - 0.26 * d - 0.20 * melt;
      const stretch = 1 + 0.30 * d;
      const puff = 1 + 0.70 * vap;

      tmpS.set(baseR * stretch * puff, baseR * squash * puff, baseR * (0.92 + 0.22 * d) * puff);

      tmpMat.compose(tmpP, tmpQ, tmpS);
      inst.setMatrixAt(drawCount, tmpMat);

      const mat = params.materials[materialId[i]];
      const e = emissive[i];

      // base color는 “암석/금속/얼음”에 맞게, 고온은 glow 쪽으로
      const c = mat.baseColor.clone().lerp(mat.glowColor, Math.min(1, 0.85 * e + 0.30 * vap));
      instanceColor.setXYZ(drawCount, c.x, c.y, c.z);

      // ✅ per-instance emissive: 고온일수록 bloom
      aEmissive.setX(drawCount, Math.min(1, 0.08 + 1.25 * e + 0.6 * vap));

      // ✅ rough/metal: 재료별 룩 차이 + 용융되면 거칠기 감소
      const rough = Math.max(0.08, Math.min(1, mat.roughnessBase * (1 - 0.45 * melt)));
      const metal = Math.max(0.0, Math.min(1, mat.metalnessBase));

      aRough.setX(drawCount, rough);
      aMetal.setX(drawCount, metal);

      // ✅ seed: 같은 재료라도 표면 패턴이 전부 달라지게
      // index 기반 + 재료 기반 섞어서 안정적으로 고정
      const seed = ((i * 1103515245 + materialId[i] * 12345) >>> 0) / 4294967295;
      aSeed.setX(drawCount, seed);

      drawCount++;
    }

    inst.count = drawCount;
    inst.instanceMatrix.needsUpdate = true;
    inst.instanceColor.needsUpdate = true;
    aSeed.needsUpdate = true;
    aEmissive.needsUpdate = true;
    aRough.needsUpdate = true;
    aMetal.needsUpdate = true;

    // 셰이더 유니폼 업데이트(런타임 UI)
    const sh = bodyMat.userData._shader;
    if (sh) {
      sh.uniforms.uDisp.value = params.surfaceDisplacement;
      sh.uniforms.uNoiseScale.value = params.surfaceNoiseScale;
      sh.uniforms.uEmissiveBoost.value = params.emissiveBoost;
    }
  }

  function render(dtReal) {
    // ✅ Add Body Mode일 때 카메라 회전/이동 잠금
    controls.enabled = !params.addBodyMode;
    controls.update();

    // trails 페이드용 dt 공유
    params._dtRealForFade = dtReal;

    particleSys.update(dtReal, 0.0, params.wind);

    bloom.strength = params.bloomStrength;
    bloom.radius = params.bloomRadius;
    bloom.threshold = params.bloomThreshold;

    renderer.toneMappingExposure = params.exposure;

    updateBodyInstances();
    trailSys.updateTrails(scene, bodies);

    composer.render();
  }

  return { scene, camera, renderer, composer, controls, arrow, trailSys, particleSys, render, resize };
}
