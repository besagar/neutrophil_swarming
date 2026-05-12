// Fixed-capacity ring buffer over Float64Array. O(1) push, O(n) snapshot.
// Multiple parallel channels share the same write index.

export function makeRing(capacity, channels = 1) {
  const cap = capacity;
  const bufs = Array.from({ length: channels }, () => new Float64Array(cap));
  let n = 0;          // total writes
  return {
    capacity: cap,
    channels,
    get count() { return Math.min(n, cap); },
    get total() { return n; },
    push(...vals) {
      const i = n % cap;
      for (let c = 0; c < channels; c++) bufs[c][i] = vals[c];
      n++;
    },
    clear() { n = 0; for (const b of bufs) b.fill(0); },
    // Snapshot oldest→newest into provided typed arrays (reused across frames).
    snapshot(out) {
      const c = Math.min(n, cap);
      const start = n >= cap ? n % cap : 0;
      for (let ch = 0; ch < channels; ch++) {
        const src = bufs[ch];
        const dst = out[ch];
        if (start + c <= cap) {
          dst.set(src.subarray(start, start + c));
        } else {
          const head = cap - start;
          dst.set(src.subarray(start, cap), 0);
          dst.set(src.subarray(0, c - head), head);
        }
      }
      return c;
    },
    // Allocate snapshot buffers (one per channel).
    allocSnapshot() {
      return Array.from({ length: channels }, () => new Float64Array(cap));
    },
  };
}
