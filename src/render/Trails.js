import * as THREE from "three";

export class Trails {
  constructor(scene){
    this.scene = scene;
    this.lines = new Map(); // id -> Line
    this.maxPoints = 256;
    this.enabled = true;
  }

  setEnabled(v){
    this.enabled = v;
    for (const line of this.lines.values()){
      line.visible = v;
    }
  }

  _ensure(id){
    let line = this.lines.get(id);
    if (line) return line;

    const geo = new THREE.BufferGeometry();
    const pos = new Float32Array(this.maxPoints * 3);
    geo.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    geo.setDrawRange(0, 0);

    const mat = new THREE.LineBasicMaterial({ transparent: true, opacity: 0.45 });
    line = new THREE.Line(geo, mat);
    line.frustumCulled = false;
    line.visible = this.enabled;
    this.scene.add(line);
    this.lines.set(id, line);
    return line;
  }

  sync(bodies){
    const alive = new Set(bodies.map(b=>b.id));
    for (const [id, line] of this.lines){
      if (!alive.has(id)){
        this.scene.remove(line);
        line.geometry.dispose();
        line.material.dispose();
        this.lines.delete(id);
      }
    }

    for (const b of bodies){
      const line = this._ensure(b.id);
      const attr = line.geometry.getAttribute("position");
      const pts = b.trail;
      const start = Math.max(0, pts.length - this.maxPoints);
      const count = pts.length - start;
      for (let i=0;i<count;i++){
        const p = pts[start+i];
        attr.setXYZ(i, p.x, p.y, p.z);
      }
      attr.needsUpdate = true;
      line.geometry.setDrawRange(0, count);
    }
  }
}
