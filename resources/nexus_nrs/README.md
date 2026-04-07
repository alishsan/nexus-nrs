# `resources/nexus_nrs/`

## `X4sGetSubent.html`

Saved HTML export from the IAEA Nuclear Data Section **EXFOR** retrieval tool (**X4sGetSubent**, `www-nds.iaea.org`), generated 2026-03-30.

| | |
|---|---|
| **EXFOR entry** | **O1198** |
| **Data subentry** | **O1198006** |
| **Reaction** | **¹⁶O(p,elastic)¹⁶O** — differential cross section **dσ/dΩ** vs **θ_CM** |
| **E_lab (p)** | **35.2 MeV** |
| **Quantity** | `DATA-CM` in **mb/sr** (CM) |

**Reference:** E. Fabrici *et al.*, *Phys. Rev. C* **21**, **844** (1980); companion paper **21**, **830** (1980). DOI [10.1103/PhysRevC.21.844](https://doi.org/10.1103/PhysRevC.21.844).

The numeric table is also embedded in **`examples/example_16Opp_elastic.clj`** as **`exfor-o1198006-16opp-dsigma`** for comparison with the code’s optical model.

**Note:** Linked scripts (`/exfor/x4js/…`) and CSS point at the IAEA site; open locally for static text only, or use the example for numbers.

## `standard_optical_potentials.json`

Schematic / tutorial Woods–Saxon vectors — see `src/nexus_nrs/optical_potentials.clj` (`load-standard-optical-potentials`).
