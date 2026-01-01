import * as THREE from "three";

/**
 * Barnes-Hut octree for approximating gravitational forces in O(N log N).
 * This implementation is tuned for moderate N (~1k). For huge N you still need more optimizations.
 */
class Node {
  constructor(center, halfSize){
    this.center = center.clone();
    this.halfSize = halfSize;
    this.children = null; // length 8
    this.mass = 0;
    this.com = new THREE.Vector3(); // center of mass
    this.body = null; // leaf body reference
  }
  isLeaf(){ return this.children === null; }
}

function childIndex(center, p){
  let idx = 0;
  if (p.x >= center.x) idx |= 1;
  if (p.y >= center.y) idx |= 2;
  if (p.z >= center.z) idx |= 4;
  return idx;
}

function childCenter(center, halfSize, idx){
  const q = halfSize * 0.5;
  return new THREE.Vector3(
    center.x + (idx & 1 ? q : -q),
    center.y + (idx & 2 ? q : -q),
    center.z + (idx & 4 ? q : -q),
  );
}

export class BarnesHutOctree {
  constructor({ center=new THREE.Vector3(), halfSize=1e6, theta=0.6 } = {}){
    this.root = new Node(center, halfSize);
    this.theta = theta;
  }

  rebuild(bodies){
    // Reset root
    const c = this.root.center.clone();
    const h = this.root.halfSize;
    this.root = new Node(c, h);

    for (const b of bodies){
      if (!b.alive) continue;
      this._insert(this.root, b);
    }
    this._accumulate(this.root);
  }

  _subdivide(node){
    node.children = new Array(8);
    for (let i=0;i<8;i++){
      node.children[i] = new Node(childCenter(node.center, node.halfSize, i), node.halfSize*0.5);
    }
  }

  _insert(node, body){
    if (node.isLeaf()){
      if (node.body === null){
        node.body = body;
        return;
      }
      // already has a body -> subdivide and reinsert
      const old = node.body;
      node.body = null;
      this._subdivide(node);
      this._insert(node, old);
      this._insert(node, body);
      return;
    }
    const idx = childIndex(node.center, body.position);
    this._insert(node.children[idx], body);
  }

  _accumulate(node){
    if (node.isLeaf()){
      if (node.body){
        node.mass = node.body.mass;
        node.com.copy(node.body.position);
      } else {
        node.mass = 0;
        node.com.set(0,0,0);
      }
      return;
    }
    let mass = 0;
    const com = new THREE.Vector3();
    for (const ch of node.children){
      this._accumulate(ch);
      if (ch.mass > 0){
        com.addScaledVector(ch.com, ch.mass);
        mass += ch.mass;
      }
    }
    node.mass = mass;
    if (mass > 0) com.multiplyScalar(1/mass);
    node.com.copy(com);
  }

  /**
   * Adds gravitational acceleration contribution from the tree to `outAccel`.
   */
  accumulateAccel(body, G, softening, outAccel){
    const stack = [this.root];
    const p = body.position;
    while (stack.length){
      const node = stack.pop();
      if (!node || node.mass <= 0) continue;
      if (node.isLeaf()){
        const other = node.body;
        if (!other || other === body) continue;
        const dx = other.position.x - p.x;
        const dy = other.position.y - p.y;
        const dz = other.position.z - p.z;
        const r2 = dx*dx + dy*dy + dz*dz + softening*softening;
        const invR = 1/Math.sqrt(r2);
        const invR3 = invR*invR*invR;
        const a = G * other.mass * invR3;
        outAccel.x += dx * a;
        outAccel.y += dy * a;
        outAccel.z += dz * a;
        continue;
      }
      // Size / distance criterion
      const dx = node.com.x - p.x;
      const dy = node.com.y - p.y;
      const dz = node.com.z - p.z;
      const dist = Math.sqrt(dx*dx + dy*dy + dz*dz) + 1e-9;
      const s = node.halfSize * 2; // node size
      if ((s / dist) < this.theta){
        const r2 = dist*dist + softening*softening;
        const invR = 1/Math.sqrt(r2);
        const invR3 = invR*invR*invR;
        const a = G * node.mass * invR3;
        outAccel.x += dx * a;
        outAccel.y += dy * a;
        outAccel.z += dz * a;
      } else {
        for (const ch of node.children) stack.push(ch);
      }
    }
  }
}
