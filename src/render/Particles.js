import * as THREE from "three";
import { clamp, randRange } from "../utils/math.js";

function makeCircleTex(){
  const size = 64;
  const canvas = document.createElement("canvas");
  canvas.width = canvas.height = size;
  const ctx = canvas.getContext("2d");
  const g = ctx.createRadialGradient(size/2, size/2, 0, size/2, size/2, size/2);
  g.addColorStop(0, "rgba(255,255,255,1)");
  g.addColorStop(0.2, "rgba(255,255,255,0.9)");
  g.addColorStop(0.5, "rgba(255,255,255,0.25)");
  g.addColorStop(1, "rgba(255,255,255,0)");
  ctx.fillStyle = g;
  ctx.fillRect(0,0,size,size);
  const tex = new THREE.CanvasTexture(canvas);
  tex.needsUpdate = true;
  return tex;
}

export class ParticleSystem {
  constructor(scene){
    this.scene = scene;
    this.max = 6000;
    this.alive = 0;

    this.pos = new Float32Array(this.max*3);
    this.vel = new Float32Array(this.max*3);
    this.life = new Float32Array(this.max);
    this.size = new Float32Array(this.max);
    this.kind = new Float32Array(this.max); // 0 dust, 1 spark, 2 flame

    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(this.pos, 3));
    geo.setAttribute("aSize", new THREE.BufferAttribute(this.size, 1));
    geo.setAttribute("aLife", new THREE.BufferAttribute(this.life, 1));
    geo.setAttribute("aKind", new THREE.BufferAttribute(this.kind, 1));

    this.tex = makeCircleTex();

    const mat = new THREE.ShaderMaterial({
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
      uniforms: {
        uTime: { value: 0 },
        uTex: { value: this.tex }
      },
      vertexShader: /* glsl */`
        attribute float aSize;
        attribute float aLife;
        attribute float aKind;
        varying float vLife;
        varying float vKind;
        void main(){
          vLife = aLife;
          vKind = aKind;
          vec4 mv = modelViewMatrix * vec4(position, 1.0);
          gl_PointSize = aSize * (300.0 / -mv.z);
          gl_Position = projectionMatrix * mv;
        }
      `,
      fragmentShader: /* glsl */`
        uniform sampler2D uTex;
        varying float vLife;
        varying float vKind;
        void main(){
          vec2 uv = gl_PointCoord;
          vec4 t = texture2D(uTex, uv);
          float a = t.a * smoothstep(0.0, 0.12, vLife) * smoothstep(1.0, 0.0, 1.0 - vLife);
          vec3 col = vec3(0.65,0.75,1.0);
          if (vKind > 1.5) col = vec3(1.0, 0.55, 0.08);
          else if (vKind > 0.5) col = vec3(1.0, 0.9, 0.6);
          else col = vec3(0.75,0.85,1.0);
          gl_FragColor = vec4(col, a);
        }
      `
    });

    this.points = new THREE.Points(geo, mat);
    this.points.frustumCulled = false;
    scene.add(this.points);
  }

  emit({ origin, count=100, spread=10, baseVelocity=new THREE.Vector3(), normal=null, energy=0.5, kind="dust" }){
    const k = (kind === "catastrophic") ? 1 : (kind === "mergeDust" ? 0 : 0);
    const kindVal = (kind === "catastrophic") ? 1 : 0;

    for (let i=0;i<count;i++){
      if (this.alive >= this.max) break;
      const idx = this.alive++;

      const dir = randomUnit();
      if (normal) dir.addScaledVector(normal, 0.7).normalize();
      const r = spread * (0.1 + Math.random());
      const px = origin.x + dir.x * r;
      const py = origin.y + dir.y * r;
      const pz = origin.z + dir.z * r;

      const spd = (3 + 45*Math.random()) * (0.2 + energy);
      const vx = baseVelocity.x + dir.x * spd;
      const vy = baseVelocity.y + dir.y * spd;
      const vz = baseVelocity.z + dir.z * spd;

      this.pos[idx*3+0]=px; this.pos[idx*3+1]=py; this.pos[idx*3+2]=pz;
      this.vel[idx*3+0]=vx; this.vel[idx*3+1]=vy; this.vel[idx*3+2]=vz;

      this.life[idx]=1.0;
      this.size[idx]= (kind==="mergeDust" ? randRange(3,7) : randRange(2,9));
      this.kind[idx]= (kind==="catastrophic" ? randRange(0.7,1.6) : randRange(0.0,0.6));
    }
    this.points.geometry.attributes.position.needsUpdate = true;
    this.points.geometry.attributes.aSize.needsUpdate = true;
    this.points.geometry.attributes.aLife.needsUpdate = true;
    this.points.geometry.attributes.aKind.needsUpdate = true;
  }

  emitFlamesAround(body){
    // Burn effect: spawn a small set of particles around the hot body each frame.
    const count = Math.floor(6 + 18*Math.min(1, (body.temperatureK-1500)/1200));
    const origin = body.position;
    for (let i=0;i<count;i++){
      if (this.alive >= this.max) break;
      const idx = this.alive++;
      const dir = randomUnit();
      const r = body.radius * (0.85 + 0.35*Math.random());
      const px = origin.x + dir.x * r;
      const py = origin.y + dir.y * r;
      const pz = origin.z + dir.z * r;

      // Blow flames backwards relative to velocity
      const back = body.velocity.clone().normalize().multiplyScalar(-1);
      const d2 = dir.clone().addScaledVector(back, 0.8).normalize();

      const spd = randRange(4, 18) + 0.03*body.temperatureK;
      const vx = body.velocity.x + d2.x * spd;
      const vy = body.velocity.y + d2.y * spd;
      const vz = body.velocity.z + d2.z * spd;

      this.pos[idx*3+0]=px; this.pos[idx*3+1]=py; this.pos[idx*3+2]=pz;
      this.vel[idx*3+0]=vx; this.vel[idx*3+1]=vy; this.vel[idx*3+2]=vz;

      this.life[idx]=randRange(0.35, 1.0);
      this.size[idx]=randRange(10, 26);
      this.kind[idx]=2.0;
    }
  }

  step(dt, time){
    const pos = this.pos, vel=this.vel, life=this.life, size=this.size;
    let w=0;
    for (let i=0;i<this.alive;i++){
      let L = life[i] - dt*randRange(0.35, 0.95);
      if (L <= 0) continue;

      const ix=i*3;
      // drift
      pos[ix+0] += vel[ix+0]*dt;
      pos[ix+1] += vel[ix+1]*dt;
      pos[ix+2] += vel[ix+2]*dt;

      // drag
      vel[ix+0] *= 0.992;
      vel[ix+1] *= 0.992;
      vel[ix+2] *= 0.992;

      // lift for flames
      if (this.kind[i] > 1.5){
        vel[ix+1] += 6.0*dt;
      }

      // compact alive particles
      if (w !== i){
        const wx=w*3;
        pos[wx+0]=pos[ix+0]; pos[wx+1]=pos[ix+1]; pos[wx+2]=pos[ix+2];
        vel[wx+0]=vel[ix+0]; vel[wx+1]=vel[ix+1]; vel[wx+2]=vel[ix+2];
        life[w]=L;
        size[w]=size[i];
        this.kind[w]=this.kind[i];
      } else {
        life[i]=L;
      }
      w++;
    }
    this.alive = w;
    this.points.material.uniforms.uTime.value = time;

    this.points.geometry.attributes.position.needsUpdate = true;
    this.points.geometry.attributes.aLife.needsUpdate = true;
    this.points.geometry.attributes.aKind.needsUpdate = true;
  }
}

function randomUnit(){
  const u = Math.random()*2-1;
  const t = Math.random()*Math.PI*2;
  const r = Math.sqrt(1-u*u);
  return new THREE.Vector3(r*Math.cos(t), u, r*Math.sin(t));
}
