import { human } from "../utils/math.js";
import * as THREE from "three";

export class UI {
  constructor(sim, renderer){
    this.sim = sim;
    this.renderer = renderer;

    this.el = document.createElement("div");
    this.el.className = "hud";
    document.getElementById("app").appendChild(this.el);

    this.left = document.createElement("div");
    this.left.className = "panel";
    this.left.innerHTML = `
      <h1>N-body Simulator <span class="badge">Planet-scale • collisions • heat</span></h1>
      <div class="row">
        <button class="btn primary" data-act="pause">Pause</button>
        <button class="btn" data-act="reset">Reset</button>
        <button class="btn" data-act="random">Random</button>
      </div>
      <div class="row">
        <button class="btn" data-act="addRock">+ Rock</button>
        <button class="btn" data-act="addIce">+ Ice</button>
        <button class="btn" data-act="addGas">+ Gas</button>
      </div>
      <div class="sep"></div>
      <div class="small">
        • 좌클릭: 천체 선택 / 드래그: 카메라 회전 / 휠: 줌<br/>
        • SHIFT+클릭: 그 지점에 새 천체 투하(질량/속도는 우측 패널에서 조절)<br/>
        • 충돌: 탄성/파편/재응집(merge) + 발열/연소(VFX)<br/>
      </div>
      <div class="kv" id="kv"></div>
    `;
    this.el.appendChild(this.left);

    this.right = document.createElement("div");
    this.right.className = "rightPanel";
    this.right.innerHTML = `
      <h2>Selected Body <span class="badge" id="selName">None</span></h2>
      <div class="mini">
        <div class="title">Properties</div>
        <div class="grid" id="props"></div>
      </div>
      <div class="sep"></div>
      <div class="mini">
        <div class="title">Spawn Settings (SHIFT+Click)</div>
        <div class="grid">
          <div class="item">Mass <b id="spMass">—</b></div>
          <div class="item">Radius <b id="spRad">—</b></div>
          <div class="item">Speed <b id="spSpd">—</b></div>
          <div class="item">Dir <b id="spDir">—</b></div>
        </div>
        <div class="row">
          <button class="btn" data-act="spawnMassDown">Mass -</button>
          <button class="btn" data-act="spawnMassUp">Mass +</button>
        </div>
        <div class="row">
          <button class="btn" data-act="spawnSpeedDown">Speed -</button>
          <button class="btn" data-act="spawnSpeedUp">Speed +</button>
        </div>
      </div>
    `;
    this.el.appendChild(this.right);

    this.footer = document.createElement("div");
    this.footer.className = "footerTag";
    this.footer.textContent = "JavaScript N-body Implementation v0.1 (Barnes-Hut + fragmentation + thermal FX)";
    this.el.appendChild(this.footer);

    this.kv = this.left.querySelector("#kv");
    this.selNameEl = this.right.querySelector("#selName");
    this.propsEl = this.right.querySelector("#props");

    // spawn settings
    this.spawn = {
      mass: 2e12,
      speed: 18,
      dir: new THREE.Vector3(1,0,0)
    };
    this._renderSpawn();

    this.selected = null;
    this.paused = false;

    this.el.addEventListener("click", (e)=>{
      const btn = e.target.closest("button[data-act]");
      if (!btn) return;
      this._act(btn.dataset.act);
    });
  }

  _act(act){
    const sim = this.sim;
    switch(act){
      case "pause":
        this.paused = !this.paused;
        this.left.querySelector('[data-act="pause"]').textContent = this.paused ? "Resume" : "Pause";
        break;
      case "reset":
        sim.reset();
        this.renderer.seedDefaultSystem();
        this.selected = null;
        break;
      case "random":
        sim.reset();
        this.renderer.seedRandomSystem();
        this.selected = null;
        break;
      case "addRock": this.renderer.spawnSimple("rock"); break;
      case "addIce":  this.renderer.spawnSimple("ice"); break;
      case "addGas":  this.renderer.spawnSimple("gas"); break;

      case "spawnMassDown": this.spawn.mass *= 0.6; this._renderSpawn(); break;
      case "spawnMassUp": this.spawn.mass *= 1.6; this._renderSpawn(); break;
      case "spawnSpeedDown": this.spawn.speed *= 0.7; this._renderSpawn(); break;
      case "spawnSpeedUp": this.spawn.speed *= 1.3; this._renderSpawn(); break;
    }
  }

  setSelected(body){
    this.selected = body;
    this.selNameEl.textContent = body ? body.name : "None";
  }

  _renderSpawn(){
    const r = this.right;
    r.querySelector("#spMass").textContent = human(this.spawn.mass);
    r.querySelector("#spSpd").textContent = this.spawn.speed.toFixed(1);
    r.querySelector("#spDir").textContent = `${this.spawn.dir.x.toFixed(1)},${this.spawn.dir.y.toFixed(1)},${this.spawn.dir.z.toFixed(1)}`;
    // radius shown as if rock density
    const dens = 3500;
    const rad = Math.cbrt((3*this.spawn.mass)/(4*Math.PI*dens));
    r.querySelector("#spRad").textContent = rad.toFixed(1);
  }

  updateStats(fps, sim){
    const alive = sim.bodies.length;
    const now = sim.time;

    this.kv.innerHTML = `
      <div>FPS</div><b>${fps.toFixed(0)}</b>
      <div>Sim Time</div><b>${(now/3600).toFixed(2)} h</b>
      <div>Bodies</div><b>${alive}</b>
      <div>G</div><b>${sim.params.G.toExponential(2)}</b>
      <div>Time Scale</div><b>${sim.params.timeScale.toFixed(0)}x</b>
      <div>Integrator</div><b>${sim.params.integrator}</b>
    `;

    if (this.selected){
      const b = this.selected;
      const v = b.velocity.length();
      const dSun = b.position.length();
      const props = [
        ["Mass", human(b.mass)],
        ["Radius", b.radius.toFixed(2)],
        ["Speed", v.toFixed(2)],
        ["Temp (K)", b.temperatureK.toFixed(0)],
        ["Material", b.material],
        ["Damage", (b.damage*100).toFixed(0)+"%"],
        ["Craters", String(b.craters.length)],
        ["Dist", dSun.toFixed(1)]
      ];
      this.propsEl.innerHTML = props.map(([k,val]) => `<div class="item">${k} <b>${val}</b></div>`).join("");
    } else {
      this.propsEl.innerHTML = `<div class="item">Click a body to inspect.</div>`;
    }
  }
}
