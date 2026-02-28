# Nexus-NRS (Nuclear Reaction Suite)

A Clojure library for nuclear reaction calculations: scattering, transfer (DWBA), inelastic reactions, and halo nuclei. Uses Riccati–Numerov integration and optical potentials for distorted waves.

## Features

- **Scattering**: S-matrix, K-matrix, phase shifts, R-matrix
- **Transfer reactions**: DWBA (POST), bound-state overlap, global optical potentials
- **Inelastic scattering**: Form factors, coupled channels
- **Halo nuclei**: Bound states, ANC extraction, Numerov hybrid vs finite start comparison
- **Numerics**: Riccati–Hankel/Bessel starts, Wronskian conservation
- **Web dashboard**: Interactive parameters and plots (elastic, inelastic, transfer)

## Installation

- **Java** 8+
- **Leiningen** 2.0+

```bash
git clone https://github.com/alishsan/nexus-nrs.git
cd nexus-nrs
lein deps
```

## Usage

- **REPL**: `lein repl` (starts in `dwba.core`).
- **Tests**: `lein test`.
- **Examples**: See `examples/` (halo nuclei, 16O(p,d), 11Li elastic/inelastic, etc.).
- **Paper data**: `(load-file "generate_paper_data.clj")`, `(load-file "generate_halo_numerov_comparison.clj")`.
- **Web dashboard**: `./start-dashboard.sh` or `cd web-dashboard && lein run` → http://localhost:3000

## License

EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0 (see [LICENSE](LICENSE)).
