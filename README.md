# Nexus-NRS (Nuclear Reaction Suite)

A Clojure library for nuclear reaction calculations: scattering, transfer (DWBA), inelastic reactions, and halo nuclei. Uses Riccati–Numerov integration and optical potentials for distorted waves.

## Features

- **Scattering**: S-matrix, K-matrix, phase shifts, R-matrix
- **Transfer reactions**: DWBA (POST), bound-state overlap, global optical potentials
- **Inelastic scattering**: Form factors, coupled channels
- **Halo nuclei**: Bound states, ANC extraction, Numerov hybrid vs finite start comparison
- **Numerics**: Riccati–Hankel/Bessel starts, Wronskian conservation
- **Web dashboard**: Interactive parameters and plots (elastic, inelastic, transfer) — [nexus-nrs.uz](https://www.nexus-nrs.uz)

### Design stance

**Correct DWBA (and related) physics in code** comes first. External codes or tabulated listings (e.g. DWUCK4) are useful **regression references** and for spotting bugs or missing mechanisms — not something to mimic with **ad hoc** angular factors or Coulomb-on-σ shortcuts. When shapes or scales disagree, the path forward is **more faithful theory in the amplitudes** (distorted waves, finite range, partial waves, coupling), not hacks tuned to a listing.

## Installation

- **Java** 8+
- **Leiningen** 2.0+

```bash
git clone https://github.com/alishsan/nexus-nrs.git
cd nexus-nrs
lein deps
```

### REPL error: `fastmath/special$Si (wrong name: fastmath/special$si)`

This comes from a **bad or mixed `fastmath` JAR** (often **`3.0.0-alpha4-SNAPSHOT`**) on **case-insensitive** disks (default macOS APFS): two inner classes differ only by case (`Si` vs `si`).

1. This repo pins **`[generateme/fastmath "3.0.0-alpha4"]`** (non-SNAPSHOT). After pulling, run:
   ```bash
   lein clean
   rm -rf ~/.m2/repository/generateme/fastmath/3.0.0-alpha4-SNAPSHOT
   lein deps
   lein repl
   ```
2. If it persists, remove the whole cache folder `~/.m2/repository/generateme/fastmath/` and `lein deps` again.

## Usage

- **REPL**: `lein repl` (starts in `dwba.core`).
- **Tests**: `lein test`.
- **DWUCK-style benchmark**: Ca40(d,p) kinematics vs DW4TST.DAT / listing — namespace `dwba.benchmark.ca40-dwuck`, tests `lein test :only dwba.ca40-dwuck-benchmark-test`. In the REPL use `(require '[dwba.benchmark.ca40-dwuck :as c])` then `(c/ca40-dp-kinematics)` — **not** `(dwba/benchmark/...)` (slashes are division).
- **Examples**: See `examples/` (halo nuclei, **16O(p,d)** `example_16Opd.clj`, **16O(d,p)** `example_16Odp.clj` — uses **`dwba.benchmark.o16-dp-handbook`** (handbook ZR + Austern (5.6), like **`ca40-pd-handbook`**), **11Li(p,d)10Li** kinematics `example_11Li_pd_10Li.clj` (¹⁰Li unbound; see file header + *Phys. Lett. B* 755 (2016) 481), 11Li elastic/inelastic, etc.).
  - **Ca40(d,p) transfer vs DWUCK4**: `examples/plot_ca40_dp_dwuck.clj` embeds listing **Inelsig** as **fm²/sr** and multiplies by **`dwuck-inelsig-fm2-sr->mb-sr`** (10) for **mb/sr** curves vs Nexus. **`ca40-dp-dsigma-mb-sr`** is **DWBA only**: Coulomb in **`distorted-wave-optical`** χ, zero-range POST **`transfer-amplitude-post`**, and by default **`transfer-differential-cross-section-angular-coherent`** (θ-dependent reduced multipole angular factor). Pass **`:angular-mode :m-sum`** for orbital **m-sum**, which is **θ-flat** when **l_i=0** (constant dσ vs θ). **`ca40-dp-flux-scale-to-embedded-dwuck`** matches scale at **30° CM** (~10¹–10²; max-norm χ vs DWUCK). Optional **`:cm-asymmetry-kappa`** is exploratory **1+κ cos θ**. See `dwba.benchmark.ca40-dwuck`.
- **Paper data**: `(load-file "generate_paper_data.clj")`, `(load-file "generate_halo_numerov_comparison.clj")`.
- **Web dashboard**: **Live:** https://www.nexus-nrs.uz — **local:** `./start-dashboard.sh` or `cd web-dashboard && lein run` → http://localhost:3000

## License

EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0 (see [LICENSE](LICENSE)).
