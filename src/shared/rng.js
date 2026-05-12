// Seeded PRNG (mulberry32) + Gaussian (Box–Muller).
// Deterministic given the seed; reseed via .seed(n).

export function makeRng(seed = 1) {
  let s = seed >>> 0;
  let cachedGauss = null;

  function next() {
    s = (s + 0x6D2B79F5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  function gauss() {
    if (cachedGauss !== null) {
      const g = cachedGauss;
      cachedGauss = null;
      return g;
    }
    let u1, u2;
    do { u1 = next(); } while (u1 < 1e-12);
    u2 = next();
    const r = Math.sqrt(-2 * Math.log(u1));
    const th = 2 * Math.PI * u2;
    cachedGauss = r * Math.sin(th);
    return r * Math.cos(th);
  }

  return {
    uniform: next,
    gauss,
    seed(n) { s = n >>> 0; cachedGauss = null; },
  };
}
