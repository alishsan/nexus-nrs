# Benchmark Reactions for Paper

This document lists all benchmark reactions used for validation and comparison in the paper.

## Overview

The benchmark reactions cover:
- **Transfer reactions** (DWBA POST formulation)
- **Elastic scattering** (S-matrix formalism)
- **Inelastic scattering** (DWBA with transition form factors)
- **Halo nuclei** (bound states, ANC extraction, transfer reactions)

---

## 1. Transfer Reactions

### 1.1 16O(p,d)15O

**Reaction Type**: Neutron pickup (p,d)

**Parameters**:
- Incident energy: E_lab = 20.0 MeV
- Initial bound state: neutron in 16O (l=1, E=-15.67 MeV)
- Final bound state: neutron in deuteron (l=0, E=-2.214 MeV)
- Q-value: negative (endothermic pickup reaction)

**Nuclear Parameters**:
- 16O bound state: V₀ = 62 MeV, R₀ = 2.7 fm, a₀ = 0.6 fm, l=1
- Deuteron bound state: V₀ = 50 MeV, R₀ = 1.5 fm, a₀ = 0.6 fm, l=0

**Optical Potentials**:
- Entrance channel: CH89 (proton on 16O)
- Exit channel: Daehnick80 (deuteron on 15O)

**Key Results**:
- Normalized overlap (form factor)
- Transfer amplitudes T_L for L = 0 to 7
- Differential cross section dσ/dΩ(θ)
- Expected: dσ/dΩ(70°) ≈ 0.5 mb/sr at 20 MeV

**File**: `examples/example_16Opd.clj`

---

### 1.2 ¹⁰Be(d,p)¹¹Be

**Reaction Type**: Neutron transfer (d,p)

**Parameters**:
- Lab energy: E_lab = 10.0 MeV
- Q-value: Q = 4.49 MeV
- Final state: ¹¹Be halo nucleus (E_b = 504 keV, l=0, S_{1/2})

**Nuclear Parameters**:
- Target: ¹⁰Be (A=10, Z=4)
- Residual: ¹¹Be (A=11, Z=4)
- Entrance: deuteron + ¹⁰Be
- Exit: proton + ¹¹Be

**Optical Potentials**:
- Entrance channel: Daehnick80 (deuteron)
- Exit channel: CH89 (proton)

**Key Results**:
- Bound state wavefunction for ¹¹Be (halo nucleus)
- Transfer amplitude (POST formulation)
- Differential cross section
- Comparison with experimental data

**File**: `examples/halo_nuclei_examples.clj` (Example 7)

---

## 2. Elastic Scattering

### 2.1 11Li(d,d) Elastic Scattering

**Reaction Type**: Elastic scattering

**Parameters**:
- Beam energy: E/A = 7.1 MeV → E_d = 14.2 MeV (lab)
- Target: 11Li (halo nucleus)

**Nuclear Parameters**:
- Woods-Saxon: V₀ = 50.0 MeV, R₀ = 2.5 fm, a₀ = 0.6 fm
- Reduced mass: μ = 1585.5 MeV/c² (d+11Li)
- Coulomb: Z₁Z₂ = 1×3 = 3

**Key Results**:
- S-matrix for L = 0 to 10
- Phase shifts δ_L
- Differential cross section dσ/dΩ(θ)
- Total cross section
- Partial wave contributions
- Angular distribution (forward peak, oscillations)

**File**: `examples/example_11Li_dd_elastic.clj`

---

## 3. Inelastic Scattering

### 3.1 11Li(d,d') Monopole Excitation

**Reaction Type**: Inelastic scattering (L=0 monopole "breathing" mode)

**Parameters**:
- Beam energy: E/A = 7.1 MeV → E_d = 14.2 MeV (lab)
- Excitation energy: E_ex = 2.09 MeV
- Multipolarity: λ = 0 (monopole)
- Ground state: 11Li is 3/2⁻

**Nuclear Parameters**:
- Woods-Saxon: V₀ = 50.0 MeV, R₀ = 2.5 fm, a₀ = 0.6 fm
- Deformation parameter: β₀ = 0.10 (monopole)
- Reduced mass: μ = 1585.5 MeV/c²

**Key Results**:
- Transition form factor F₀(r) = β₀ · R₀ · dV/dr
- Inelastic scattering amplitude T_inel
- Differential cross section dσ/dΩ
- Angular distribution (should be isotropic for L=0)
- Legendre polynomial expansion (only P₀ term)

**Reference**: "A measurement of scattering of 11Li on a deuteron target at beam energy of E/A = 7.1 MeV (at the center of the target) was performed at the IRIS facility at TRIUMF. The differential cross section of the 2.09 MeV state when interpreted by distorted wave Born approximation calculation (DWBA) shows the first evidence of soft monopole (L = 0) excitation."

**File**: `examples/example_11Li_dd_inelastic.clj`

---

### 3.2 ¹²C(α,α')¹²C* (2⁺ State)

**Reaction Type**: Inelastic scattering (quadrupole excitation)

**Parameters**:
- Incident energy: E = 10.0 MeV
- Excitation energy: E_ex = 4.44 MeV (first 2⁺ state)
- Multipolarity: λ = 2 (quadrupole)
- Deformation parameter: β₂ = 0.25

**Nuclear Parameters**:
- Woods-Saxon: V₀ = 50.0 MeV, R₀ = 2.0 fm, a₀ = 0.6 fm
- Target: ¹²C

**Key Results**:
- Transition form factor F₂(r)
- Reduced matrix element
- B(E2) transition strength
- Differential cross section

**File**: `test/dwba/inelastic_validation_test.clj`

---

## 4. Halo Nuclei

### 4.1 ¹¹Be (Neutron Halo)

**Properties**:
- Binding energy: E_b = 504 keV
- Angular momentum: l = 0
- State: S_{1/2}
- Core: ¹⁰Be
- Halo: 1 neutron

**Key Calculations**:
- Bound state wavefunction
- Adaptive matching radius
- Decay length κ⁻¹
- ANC extraction
- Transfer reaction (see 1.2)

**File**: `examples/halo_nuclei_examples.clj` (Example 1)

---

### 4.2 ⁸B (Proton Halo)

**Properties**:
- Binding energy: E_b = 137 keV
- Angular momentum: l = 2
- Core: ⁷Be
- Halo: 1 proton

**Key Calculations**:
- Bound state wavefunction (p-wave)
- Adaptive matching radius
- Decay length κ⁻¹
- Coulomb effects (charged halo)

**File**: `examples/halo_nuclei_examples.clj` (Example 2)

---

## 5. Additional Benchmark Cases

### 5.1 Low-Energy Scattering

**Test Case**: 50 keV scattering
- Energy: E = 0.05 MeV
- Tests method at very low energies near threshold
- Large radius required: r_max = 50.0 fm

**File**: `examples/halo_nuclei_examples.clj` (Example 5)

---

### 5.2 Coulomb Scattering

**Test Case**: α+p at 10 MeV
- Demonstrates Riccati-Bessel initialization with Coulomb effects
- Sommerfeld parameter η calculation
- Nuclear + Coulomb scattering

**File**: `examples/halo_nuclei_examples.clj` (Example 4)

---

## 6. Comparison Cases

### 6.1 LEA vs Zero-Range

**Comparison**: Local Energy Approximation (LEA) vs Zero-Range for (d,p) reactions
- Tests different approximations for transfer reactions
- Comparison of transfer amplitudes and cross sections

**File**: `examples/example_lea_vs_zero_range.clj`

---

### 6.2 Post vs Prior Formulations

**Comparison**: POST vs PRIOR formulations for transfer reactions
- Tests different DWBA formulations
- Comparison of transfer amplitudes

**File**: `examples/example_post_prior.clj`

---

### 6.3 Numerov Start Comparison

**Comparison**: Hybrid (Bessel) vs Finite start near zero
- Tests different Numerov initialization methods
- Important for halo nuclei with extended wavefunctions

**File**: `examples/halo_nuclei_examples.clj` (Example 7, mentioned)

---

## Summary Table

| Reaction | Type | Energy | Key Feature |
|----------|------|--------|-------------|
| 16O(p,d)15O | Transfer | 20 MeV | Neutron pickup, well-studied |
| ¹⁰Be(d,p)¹¹Be | Transfer | 10 MeV | Halo nucleus final state |
| 11Li(d,d) | Elastic | 14.2 MeV | Halo nucleus target |
| 11Li(d,d') | Inelastic | 14.2 MeV | Monopole L=0 excitation |
| ¹²C(α,α')¹²C* | Inelastic | 10 MeV | Quadrupole excitation |
| ¹¹Be | Halo | - | Neutron halo, l=0 |
| ⁸B | Halo | - | Proton halo, l=2 |

---

## Notes

1. All reactions use appropriate optical potentials (CH89 for protons, Daehnick80 for deuterons)
2. Numerical parameters: typically h = 0.01 fm, r_max = 20-30 fm
3. Mass factors and Coulomb parameters are reaction-specific
4. Results should be compared with experimental data where available
5. Convergence tests should be performed for all benchmark cases
