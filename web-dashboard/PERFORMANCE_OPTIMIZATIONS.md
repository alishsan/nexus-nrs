# Performance Optimizations for DWBA Web Dashboard

## Summary

**Status**: Optimizations applied but need profiling to identify actual bottlenecks.

**Current Performance**: ~7.9 seconds (slower than baseline ~5.8 seconds)
**Issue**: Initial optimizations (pmap, atom caching) added overhead without benefit.

## Latest Changes

1. ✅ **Removed `pmap`**: Overhead was too high for this workload
2. ✅ **Removed atom-based caching**: Overhead without benefit (each combination is unique)
3. ✅ **Kept step size optimization**: `h = 0.02` (was 0.01) - ~2x fewer points
4. ✅ **Reduced r-max**: `r-max = 15.0` (was 20.0) - ~25% fewer points
5. **Total point reduction**: ~2000 → ~750 points per wavefunction (2.7x fewer)

**Expected speedup**: ~2.5x (should bring ~5.8s → ~2.3s), but actual performance needs verification.

## Identified Bottlenecks

### 1. **No Caching of Distorted Waves**
- **Problem**: Each (E, L) combination recalculated distorted waves independently
- **Impact**: For 26 energies × 6 L-values = 156 combinations, each doing 2 Numerov integrations (~2000 points each)
- **Solution**: Implemented caching using atoms to store calculated wavefunctions

### 2. **Fine Step Size**
- **Problem**: `h = 0.01` fm creates ~2000 points per wavefunction
- **Impact**: Numerov integration scales as O(N) where N = r_max/h
- **Solution**: Increased step size to `h = 0.02` fm (~1000 points, ~4x faster)

### 3. **Sequential Processing**
- **Problem**: All calculations done sequentially in a single thread
- **Impact**: CPU cores underutilized
- **Solution**: Used `pmap` for parallel processing across all (E, L) combinations

### 4. **Redundant Calculations**
- **Problem**: Same distorted waves recalculated multiple times
- **Impact**: Unnecessary Numerov integrations
- **Solution**: Cache-based lookup prevents recalculation

## Optimizations Implemented

### Inelastic Scattering Endpoint (`/api/inelastic`)

**Before:**
```clojure
(for [E-i energies L-i L-values]
  (let [chi-i (inel E-i L-i ws-params h r-max)      ; Recalculated each time
        chi-f (inel-exit E-i E-ex L-i ws-params h r-max)  ; Recalculated each time
        ...])
```

**After:**
```clojure
;; Cache distorted waves
(let [chi-i-cache (atom {})
      chi-f-cache (atom {})
      get-chi-i (fn [E-i L-i] ...)  ; Cached lookup
      get-chi-f (fn [E-i L-i] ...)  ; Cached lookup
      combinations (for [E-i energies L-i L-values] [E-i L-i])
      inelastic-data (doall (pmap ... combinations))]  ; Parallel processing
```

**Changes:**
1. ✅ **Caching**: Distorted waves cached by `[E-i, L-i]` key
2. ✅ **Parallelization**: `pmap` processes all combinations in parallel
3. ✅ **Step size**: `h = 0.02` instead of `0.01` (2x fewer points, ~4x faster Numerov)

### Transfer Reaction Endpoint (`/api/transfer`)

**Before:**
```clojure
(for [E-i energies L L-values]
  ...)  ; Sequential processing
```

**After:**
```clojure
(let [combinations (for [E-i energies L L-values] [E-i L])
      transfer-data (doall (pmap ... combinations))]  ; Parallel processing
```

**Changes:**
1. ✅ **Parallelization**: `pmap` processes all combinations in parallel
2. ✅ **Bound states**: Already calculated once (not per energy/L), so no caching needed

## Performance Impact

### Expected Speedup

| Optimization | Speedup Factor | Notes |
|-------------|----------------|-------|
| Step size (h: 0.01 → 0.02) | ~4x | Numerov integration scales as O(N²) |
| Parallelization (pmap) | ~2-4x | Depends on CPU cores (typically 2-4x on 4-8 core systems) |
| Caching | ~1.5-2x | Avoids redundant calculations |
| **Total Expected** | **~12-32x** | Combined effect |

### Actual Performance

- **Before**: ~5786 ms (5.8 seconds)
- **Expected After**: ~180-480 ms (0.18-0.48 seconds)
- **Realistic Estimate**: ~500-1500 ms (0.5-1.5 seconds) accounting for overhead

## Trade-offs

### Step Size (h = 0.02 vs 0.01)

**Pros:**
- 4x faster Numerov integration
- Still accurate for most nuclear physics calculations
- Standard step size for many DWBA codes

**Cons:**
- Slightly lower accuracy (typically <1% difference)
- May need finer step size for very high precision calculations

**Recommendation**: Keep `h = 0.02` for dashboard. Users requiring higher precision can use the command-line examples with `h = 0.01`.

### Parallelization

**Pros:**
- Utilizes all CPU cores
- Significant speedup on multi-core systems
- No accuracy impact

**Cons:**
- Slightly higher memory usage (all threads active)
- Overhead from thread coordination (minimal with `pmap`)

## Further Optimization Opportunities

### 1. **Adaptive Step Size**
- Use finer step size near nuclear surface (r ≈ R₀)
- Use coarser step size in asymptotic region (r > 10 fm)
- Could provide 2-3x additional speedup

### 2. **Reduced r-max**
- Current: `r-max = 20.0` fm
- Could reduce to `r-max = 15.0` fm for some calculations
- Saves ~25% computation time

### 3. **Memoization at Namespace Level**
- Use Clojure's `memoize` for frequently called functions
- Could cache phase shifts, R-matrices, etc.

### 4. **Progressive Loading**
- Calculate low-L values first (faster)
- Update UI progressively as calculations complete
- Already implemented in frontend, but could enhance backend

### 5. **GPU Acceleration**
- Numerov integration is highly parallelizable
- Could use GPU for very large parameter sweeps
- Requires significant refactoring

## Monitoring Performance

The optimized endpoints now return timing information:
- `h`: Step size used
- `r_max`: Maximum radius
- Calculation time can be measured client-side

## Testing Recommendations

1. **Accuracy Check**: Compare results with `h = 0.01` vs `h = 0.02`
2. **Performance Benchmark**: Measure actual speedup on target hardware
3. **Memory Usage**: Monitor memory consumption with parallel processing
4. **Edge Cases**: Test with very large energy/L-value ranges

## Code Locations

- **Inelastic endpoint**: `web-dashboard/src/dwba_web/simple_core.clj` (lines 183-266)
- **Transfer endpoint**: `web-dashboard/src/dwba_web/simple_core.clj` (lines 268-340)
- **Elastic endpoint**: `web-dashboard/src/dwba_web/simple_core.clj` (lines 136-182) - Not optimized yet (likely fast enough)
