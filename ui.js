// ui.js
import GUI from "https://cdn.jsdelivr.net/npm/lil-gui@0.19/+esm";

export function createUI(params, callbacks) {
  const gui = new GUI({ width: 360, title: "N-Body Impact Simulator" });

  // Physics
  const fPhys = gui.addFolder("Physics (Gravity / Integrator)");
  fPhys.add(params, "G", 0.001, 8.0, 0.001);
  fPhys.add(params, "eps", 0.01, 10.0, 0.01).name("softening ε");
  fPhys.add(params, "dt", 0.001, 0.06, 0.001);
  fPhys.add(params, "substeps", 1, 8, 1);
  fPhys.add(params, "collisionPositionCorrection", 0.0, 1.0, 0.01).name("pos corr");
  fPhys.add(params, "enableCollisions");
  fPhys.add(params, "enableFracture");
  fPhys.add(params, "enableThermal");
  fPhys.open();

  // Visual
  const fVis = gui.addFolder("Visual (Bloom / Trails / AA)");
  fVis.add(params, "exposure", 0.4, 2.4, 0.01);
  fVis.add(params, "bloomStrength", 0.0, 2.5, 0.01);
  fVis.add(params, "bloomThreshold", 0.0, 1.0, 0.01);
  fVis.add(params, "bloomRadius", 0.0, 1.0, 0.01);
  fVis.add(params, "trailLength", 32, 256, 1);
  fVis.add(params, "maxTrails", 16, 256, 1);
  fVis.add(params, "starCount", 500, 15000, 100);
  fVis.open();

  // Budget / Degrade
  const fBud = gui.addFolder("Performance / Budget");
  fBud.add(params, "particleBudgetSoft", 500, params.particleBudgetHard, 50).name("particle soft cap");
  fBud.add(params, "maxFragmentsPerEvent", 8, 160, 1);
  fBud.add(params, "maxSpallFragments", 8, 160, 1);
  fBud.add(params, "autoDegrade");
  fBud.add(params, "degradeTargetFPS", 20, 60, 1);
  fBud.open();

  // Spawn
  const fSpawn = gui.addFolder("Runtime Spawn");
  fSpawn.add(params, "addBodyMode").name("Add Body Mode");
  fSpawn.add(params, "spawnMass", 0.1, 800, 0.1);
  fSpawn.add(params, "spawnRadius", 0.2, 30, 0.1);
  fSpawn.add(params, "spawnTemp", 200, 8000, 10);
  fSpawn.add(params, "spawnSpeedScale", 0.1, 60, 0.1);
  fSpawn.add(params, "spawnMaterialPreset", ["rock", "ice", "metal", "gas"]);
  fSpawn.add({ spawn: () => callbacks.spawnAtOrigin?.() }, "spawn").name("Spawn at Origin");
  fSpawn.open();

  // Presets
  const fPre = gui.addFolder("Presets");
  fPre.add({ galaxy: () => callbacks.loadPreset?.("galaxy") }, "galaxy").name("Galaxy disk");
  fPre.add({ impact: () => callbacks.loadPreset?.("impact") }, "impact").name("High-speed impact test");
  fPre.add({ reaccrete: () => callbacks.loadPreset?.("reaccrete") }, "reaccrete").name("Re-accretion test");
  fPre.open();

  // Toggle actions
  gui.add({ reset: () => callbacks.reset?.() }, "reset").name("Reset (clear)");

  // Overlay
  const overlay = document.getElementById("overlay");
  function setOverlay(text) { overlay.textContent = text; }

  return { gui, setOverlay };
}
