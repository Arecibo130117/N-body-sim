import * as THREE from "three";
import { planetVertex, planetFragment } from "./shaders.js";
import { MaterialPresets, MaterialType } from "../sim/Body.js";

export class BodyRenderer {
  constructor(scene){
    this.scene = scene;
    this.meshes = new Map(); // body.id -> mesh
    this._geoCache = new Map();
  }

  _getGeo(detail){
    if (this._geoCache.has(detail)) return this._geoCache.get(detail);
    const geo = new THREE.IcosahedronGeometry(1, detail);
    this._geoCache.set(detail, geo);
    return geo;
  }

  ensure(body){
    let m = this.meshes.get(body.id);
    if (m) return m;

    const preset = MaterialPresets[body.material] || MaterialPresets.rock;
    const detail = body.radius > 220 ? 5 : (body.radius > 90 ? 4 : 3);

    const mat = new THREE.ShaderMaterial({
      vertexShader: planetVertex,
      fragmentShader: planetFragment,
      uniforms: {
        uTime: { value: 0 },
        uRadius: { value: body.radius },
        uDamage: { value: body.damage },
        uBaseColor: { value: preset.color.clone() },
        uTempK: { value: body.temperatureK },
        uIsGas: { value: (body.material === MaterialType.GAS) ? 1.0 : 0.0 },
        uCraterCount: { value: 0 },
        uCraterDir: { value: new Array(32).fill(0).map(()=> new THREE.Vector3(0,1,0)) },
        uCraterR: { value: new Float32Array(32) },
        uCraterDepth: { value: new Float32Array(32) },
      }
    });

    m = new THREE.Mesh(this._getGeo(detail), mat);
    m.frustumCulled = true;
    m.userData.bodyId = body.id;
    this.scene.add(m);
    this.meshes.set(body.id, m);
    return m;
  }

  remove(bodyId){
    const m = this.meshes.get(bodyId);
    if (!m) return;
    this.scene.remove(m);
    m.geometry?.dispose?.();
    m.material?.dispose?.();
    this.meshes.delete(bodyId);
  }

  sync(bodies, time){
    const alive = new Set(bodies.map(b=>b.id));
    for (const [id] of this.meshes){
      if (!alive.has(id)) this.remove(id);
    }

    for (const b of bodies){
      const mesh = this.ensure(b);
      mesh.position.copy(b.position);
      mesh.scale.setScalar(b.radius);
      mesh.rotation.setFromVector3(b.spinAxis.clone().multiplyScalar(b.spinAngle));

      const u = mesh.material.uniforms;
      u.uTime.value = time;
      u.uRadius.value = b.radius;
      u.uDamage.value = b.damage;
      u.uTempK.value = b.temperatureK;

      // base color drift slightly with temperature for drama
      const preset = MaterialPresets[b.material] || MaterialPresets.rock;
      u.uBaseColor.value.copy(preset.color).lerp(new THREE.Color("#3a2c1b"), Math.min(0.4, Math.max(0, (b.temperatureK-1200)/2000)));

      // craters
      const count = Math.min(32, b.craters.length);
      u.uCraterCount.value = count;
      const dirs = u.uCraterDir.value;
      const rArr = u.uCraterR.value;
      const dArr = u.uCraterDepth.value;
      for (let i=0;i<count;i++){
        dirs[i].copy(b.craters[i].dir);
        rArr[i] = b.craters[i].r;
        dArr[i] = b.craters[i].depth * Math.exp(-b.craters[i].age*0.01); // very slow "healing"
      }
      for (let i=count;i<32;i++){
        dirs[i].set(0,1,0);
        rArr[i] = 0;
        dArr[i] = 0;
      }
    }
  }

  pick(raycaster, bodies){
    const meshes = [...this.meshes.values()];
    const hits = raycaster.intersectObjects(meshes, false);
    if (!hits.length) return null;
    const id = hits[0].object.userData.bodyId;
    return bodies.find(b=>b.id === id) || null;
  }
}
