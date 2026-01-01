export function clamp(x, a, b){ return Math.min(b, Math.max(a, x)); }
export function lerp(a,b,t){ return a + (b-a)*t; }

export function randRange(a,b){ return a + Math.random()*(b-a); }

export function seededRand(seed){
  // xorshift32
  let x = (seed|0) || 123456789;
  return () => {
    x ^= x << 13; x |= 0;
    x ^= x >>> 17;
    x ^= x << 5; x |= 0;
    return ((x >>> 0) / 4294967296);
  };
}

export function human(n){
  const abs = Math.abs(n);
  if (abs >= 1e12) return (n/1e12).toFixed(2) + " T";
  if (abs >= 1e9) return (n/1e9).toFixed(2) + " B";
  if (abs >= 1e6) return (n/1e6).toFixed(2) + " M";
  if (abs >= 1e3) return (n/1e3).toFixed(2) + " K";
  return n.toFixed(2);
}
