export const planetVertex = /* glsl */`
precision highp float;

uniform float uTime;
uniform float uRadius;
uniform float uDamage;

uniform int uCraterCount;
uniform vec3 uCraterDir[32];
uniform float uCraterR[32];
uniform float uCraterDepth[32];

varying vec3 vN;
varying vec3 vPos;
varying float vCrater;

float hash(vec3 p){
  p = fract(p * 0.3183099 + vec3(0.1, 0.2, 0.3));
  p *= 17.0;
  return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
}
float noise(vec3 p){
  vec3 i = floor(p);
  vec3 f = fract(p);
  f = f*f*(3.0-2.0*f);
  float n000 = hash(i + vec3(0,0,0));
  float n100 = hash(i + vec3(1,0,0));
  float n010 = hash(i + vec3(0,1,0));
  float n110 = hash(i + vec3(1,1,0));
  float n001 = hash(i + vec3(0,0,1));
  float n101 = hash(i + vec3(1,0,1));
  float n011 = hash(i + vec3(0,1,1));
  float n111 = hash(i + vec3(1,1,1));
  float n00 = mix(n000, n100, f.x);
  float n10 = mix(n010, n110, f.x);
  float n01 = mix(n001, n101, f.x);
  float n11 = mix(n011, n111, f.x);
  float n0 = mix(n00, n10, f.y);
  float n1 = mix(n01, n11, f.y);
  return mix(n0, n1, f.z);
}
float fbm(vec3 p){
  float a = 0.55;
  float f = 1.0;
  float n = 0.0;
  for(int i=0;i<5;i++){
    n += a*noise(p*f);
    f *= 2.0;
    a *= 0.55;
  }
  return n;
}

float craterField(vec3 nrm){
  float c = 0.0;
  for (int i=0;i<32;i++){
    if (i >= uCraterCount) break;
    float d = 1.0 - dot(nrm, normalize(uCraterDir[i])); // 0 at center
    float ang = sqrt(max(0.0, d*2.0)); // approx angle distance
    float r = uCraterR[i] / uRadius;   // relative
    float x = ang / max(1e-4, r);
    float bowl = exp(-x*x*2.0);
    c += bowl * uCraterDepth[i] / uRadius;
  }
  return c;
}

void main(){
  vec3 pos = position;
  vec3 nrm = normalize(normal);

  // base roughness via fbm
  float n = fbm(nrm*4.0 + uTime*0.03);
  float ridge = 1.0 - abs(2.0*n - 1.0);
  float rough = mix(n, ridge, 0.35);

  float c = craterField(nrm);
  vCrater = c;

  // "deformation": craters dent inward, damage adds wobble
  float wobble = (fbm(nrm*2.0 + uTime*0.12) - 0.5) * uDamage * 0.35;
  float disp = (rough - 0.5) * 0.08 + wobble - c;

  pos += nrm * (disp * uRadius);

  vN = normalize(normalMatrix * (nrm));
  vec4 wp = modelMatrix * vec4(pos, 1.0);
  vPos = wp.xyz;

  gl_Position = projectionMatrix * viewMatrix * wp;
}
`;

export const planetFragment = /* glsl */`
precision highp float;

uniform vec3 uBaseColor;
uniform float uTempK;
uniform float uDamage;
uniform float uIsGas;
uniform float uTime;

varying vec3 vN;
varying vec3 vPos;
varying float vCrater;

vec3 heatColor(float T){
  // Very rough "blackbody-ish" mapping (game look).
  // 500K: dull, 1200K: orange, 2000K: yellow-white, 3500K: white-blue
  float t = clamp((T - 500.0) / 3000.0, 0.0, 1.0);
  vec3 c1 = vec3(0.12,0.02,0.00);
  vec3 c2 = vec3(1.00,0.28,0.03);
  vec3 c3 = vec3(1.00,0.95,0.55);
  vec3 c4 = vec3(0.65,0.80,1.00);
  vec3 a = mix(c1, c2, smoothstep(0.0,0.45,t));
  vec3 b = mix(c3, c4, smoothstep(0.65,1.0,t));
  return mix(a,b, smoothstep(0.45,1.0,t));
}

float fresnel(vec3 n, vec3 v){
  return pow(1.0 - max(0.0, dot(n, v)), 4.0);
}

float hash(vec3 p){
  p = fract(p * 0.3183099 + vec3(0.1, 0.2, 0.3));
  p *= 17.0;
  return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
}
float noise(vec3 p){
  vec3 i = floor(p);
  vec3 f = fract(p);
  f = f*f*(3.0-2.0*f);
  float n000 = hash(i + vec3(0,0,0));
  float n100 = hash(i + vec3(1,0,0));
  float n010 = hash(i + vec3(0,1,0));
  float n110 = hash(i + vec3(1,1,0));
  float n001 = hash(i + vec3(0,0,1));
  float n101 = hash(i + vec3(1,0,1));
  float n011 = hash(i + vec3(0,1,1));
  float n111 = hash(i + vec3(1,1,1));
  float n00 = mix(n000, n100, f.x);
  float n10 = mix(n010, n110, f.x);
  float n01 = mix(n001, n101, f.x);
  float n11 = mix(n011, n111, f.x);
  float n0 = mix(n00, n10, f.y);
  float n1 = mix(n01, n11, f.y);
  return mix(n0, n1, f.z);
}
float fbm(vec3 p){
  float a = 0.55;
  float f = 1.0;
  float n = 0.0;
  for(int i=0;i<5;i++){
    n += a*noise(p*f);
    f *= 2.0;
    a *= 0.55;
  }
  return n;
}

void main(){
  vec3 n = normalize(vN);
  vec3 v = normalize(cameraPosition - vPos);
  vec3 l = normalize(vec3(0.6, 0.8, 0.2));
  float ndl = max(0.0, dot(n, l));

  vec3 base = uBaseColor;
  // surface darkening in craters + damage
  base *= (1.0 - 0.55*clamp(vCrater*8.0, 0.0, 1.0));
  base *= (1.0 - 0.30*uDamage);

  // pseudo spec
  float spec = pow(max(0.0, dot(reflect(-l,n), v)), 32.0) * (0.15 + 0.35*(1.0-uIsGas));

  // molten cracks when hot + damaged
  float hot = smoothstep(900.0, 2200.0, uTempK);
  float cracks = fbm(n*8.0 + uTime*0.1);
  cracks = smoothstep(0.58, 0.78, cracks) * hot * (0.25 + 0.75*uDamage);

  vec3 emiss = heatColor(uTempK) * (0.10 + 0.95*hot) * (0.15 + cracks*2.0);

  // gas giants: thicker rim + softer lighting
  float rim = fresnel(n, v) * (0.65*uIsGas);
  vec3 gasTint = vec3(0.25,0.45,1.0) * rim;

  vec3 col = base*(0.18 + 1.25*ndl) + spec + emiss + gasTint;

  // subtle haze
  col += 0.03 * rim;

  gl_FragColor = vec4(col, 1.0);
}
`;
