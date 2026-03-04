# Experimental Data Comparison Table

This document provides a comparison between calculated and experimental values for benchmark reactions.

## Legend

- **Calc**: Calculated value from Nexus-NRS
- **Exp**: Experimental value from literature
- **Ratio**: Calc/Exp ratio
- **Agreement**: Qualitative assessment (✓ = good, ~ = fair, ✗ = poor, — = not available)

---

## 1. Transfer Reactions

### 1.1 16O(p,d)15O

| Parameter | Angle/Energy | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|--------------|------------|--------------|-------|-----------|-----------|
| dσ/dΩ (mb/sr) | θ = 0° | — | — | — | — | — |
| dσ/dΩ (mb/sr) | θ = 70° | — | ~0.5 | — | Code comment | — |
| dσ/dΩ (mb/sr) | θ = 90° | — | — | — | — | — |
| E_lab (MeV) | — | 20.0 | 20.0 | — | — | ✓ |
| Q-value (MeV) | — | — | -3.97 | — | — | — |

**Notes**:
- Well-studied benchmark reaction
- Expected dσ/dΩ(70°) ≈ 0.5 mb/sr at 20 MeV (from code comment)
- Need to run calculation to fill in calculated values

---

### 1.2 ¹⁰Be(d,p)¹¹Be

| Parameter | Angle/Energy | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|--------------|------------|--------------|-------|-----------|-----------|
| dσ/dΩ (mb/sr) | θ = 0° (L=0) | — | — | — | — | — |
| dσ/dΩ (mb/sr) | θ = 90° (L=0) | — | ~0.07 | — | Code comment | — |
| E_lab (MeV) | — | 10.0 | 10.0-21.4 | — | PRC 88, 064612 (2013) | ✓ |
| Q-value (MeV) | — | 4.49 | 4.49 | — | — | ✓ |
| Spectroscopic Factor S | Ground state | — | 0.71(5) | — | PRC 88, 064612 (2013) | — |
| Spectroscopic Factor S | 1st excited | — | 0.62(4) | — | PRC 88, 064612 (2013) | — |

**Notes**:
- Transfer to halo nucleus final state (¹¹Be)
- Experimental spectroscopic factors from PRC 88, 064612 (2013)
- Typical order of magnitude: ~0.07 mb/sr (from code comment, needs verification)
- Reference: Phys. Rev. C 88, 064612 (2013) - "Reactions of a ¹⁰Be beam on proton and deuteron targets"

---

## 2. Elastic Scattering

### 2.1 11Li(d,d) Elastic Scattering

| Parameter | Angle/Energy | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|--------------|------------|--------------|-------|-----------|-----------|
| dσ/dΩ (mb/sr) | θ = 0° | — | — | — | — | — |
| dσ/dΩ (mb/sr) | θ = 90° | — | — | — | — | — |
| dσ/dΩ (mb/sr) | θ = 180° | — | — | — | — | — |
| σ_total (mb) | — | — | — | — | — | — |
| E_lab (MeV) | — | 14.2 | 14.2 (E/A=7.1) | — | — | ✓ |
| Phase shift δ₀ (rad) | L=0 | — | — | — | — | — |

**Notes**:
- Elastic scattering on halo nucleus target
- E/A = 7.1 MeV → E_d = 14.2 MeV
- Need to run calculation to fill in values
- Expected forward peak and oscillatory angular distribution

---

## 3. Inelastic Scattering

### 3.1 11Li(d,d') Monopole Excitation (L=0)

| Parameter | Angle/Energy | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|--------------|------------|--------------|-------|-----------|-----------|
| dσ/dΩ (mb/sr) | θ = 0° | — | — | — | — | — |
| dσ/dΩ (mb/sr) | θ = 90° | — | — | — | — | — |
| dσ/dΩ (mb/sr) | Isotropic (L=0) | — | — | — | — | — |
| E_ex (MeV) | — | 2.09 | 2.09 | — | TRIUMF IRIS | ✓ |
| E_lab (MeV) | — | 14.2 | 14.2 (E/A=7.1) | — | TRIUMF IRIS | ✓ |
| β₀ (deformation) | — | 0.10 | — | — | — | — |

**Notes**:
- First evidence of soft monopole (L=0) excitation
- Reference: "A measurement of scattering of 11Li on a deuteron target at beam energy of E/A = 7.1 MeV (at the center of the target) was performed at the IRIS facility at TRIUMF"
- Angular distribution should be isotropic for L=0
- Deformation parameter β₀ = 0.10 (may need adjustment to match data)

---

### 3.2 ¹²C(α,α')¹²C* (2⁺ State)

| Parameter | Angle/Energy | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|--------------|------------|--------------|-------|-----------|-----------|
| dσ/dΩ (mb/sr) | θ = 0° | — | — | — | — | — |
| dσ/dΩ (mb/sr) | θ = 90° | — | — | — | — | — |
| E_ex (MeV) | — | 4.44 | 4.44 | — | — | ✓ |
| B(E2) (e²·fm⁴) | — | — | 3-5 | — | Typical exp. | — |
| β₂ (deformation) | — | 0.25 | 0.25 | — | — | ✓ |
| E_lab (MeV) | — | 10.0 | 10.0 | — | — | ✓ |

**Notes**:
- Quadrupole excitation (λ=2)
- Typical experimental B(E2) ≈ 3-5 e²·fm⁴ for ¹²C
- Deformation parameter β₂ = 0.25

---

## 4. Halo Nuclei Properties

### 4.1 ¹¹Be (Neutron Halo)

| Parameter | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|------------|--------------|-------|-----------|-----------|
| E_b (MeV) | 0.504 | 0.504 | 1.00 | — | ✓ |
| l | 0 | 0 | — | — | ✓ |
| Decay length κ⁻¹ (fm) | — | — | — | — | — |
| ANC (fm⁻¹/²) | — | — | — | — | — |
| Matching radius (fm) | — | — | — | — | — |

**Notes**:
- Neutron halo nucleus
- Binding energy matches experimental value
- Need to calculate decay length and ANC

---

### 4.2 ⁸B (Proton Halo)

| Parameter | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|------------|--------------|-------|-----------|-----------|
| E_b (MeV) | 0.137 | 0.137 | 1.00 | — | ✓ |
| l | 2 | 2 | — | — | ✓ |
| Decay length κ⁻¹ (fm) | — | — | — | — | — |
| ANC (fm⁻¹/²) | — | — | — | — | — |
| Matching radius (fm) | — | — | — | — | — |

**Notes**:
- Proton halo nucleus (charged)
- Binding energy matches experimental value
- p-wave (l=2) halo

---

## 5. Summary Statistics

| Reaction Type | Number of Cases | Cases with Data | Cases with Good Agreement | Cases Needing Calculation |
|---------------|-----------------|-----------------|---------------------------|---------------------------|
| Transfer | 2 | 0 | 0 | 2 |
| Elastic | 1 | 0 | 0 | 1 |
| Inelastic | 2 | 0 | 0 | 2 |
| Halo Properties | 2 | 2 | 2 | 0 |
| **Total** | **7** | **2** | **2** | **5** |

---

## 6. References

1. **16O(p,d)15O**: Well-studied benchmark reaction (needs specific reference)
2. **¹⁰Be(d,p)¹¹Be**: Phys. Rev. C 88, 064612 (2013) - "Reactions of a ¹⁰Be beam on proton and deuteron targets"
3. **11Li(d,d')**: TRIUMF IRIS facility measurement (E/A = 7.1 MeV)
4. **11Li(d,d)**: Same experiment as above
5. **¹²C(α,α')**: Standard benchmark for inelastic scattering
6. **Halo nuclei properties**: Standard nuclear data tables

---

## 7. Action Items

### Calculations Needed:
- [ ] Run 16O(p,d)15O calculation and extract dσ/dΩ at various angles
- [ ] Run ¹⁰Be(d,p)¹¹Be calculation and extract dσ/dΩ
- [ ] Run 11Li(d,d) elastic scattering calculation
- [ ] Run 11Li(d,d') inelastic scattering calculation
- [ ] Run ¹²C(α,α') inelastic scattering calculation
- [ ] Calculate decay lengths and ANCs for halo nuclei

### Experimental Data Needed:
- [ ] Find specific experimental dσ/dΩ values for 16O(p,d)15O at 20 MeV
- [ ] Find experimental dσ/dΩ values for ¹⁰Be(d,p)¹¹Be at 10 MeV
- [ ] Find experimental elastic scattering data for 11Li(d,d)
- [ ] Find experimental inelastic scattering data for 11Li(d,d') at 2.09 MeV
- [ ] Find experimental B(E2) value for ¹²C(α,α')¹²C* reaction
- [ ] Find experimental ANC values for ¹¹Be and ⁸B

### Validation:
- [ ] Compare calculated vs experimental cross sections
- [ ] Calculate χ² values where multiple data points available
- [ ] Assess systematic uncertainties
- [ ] Document any discrepancies and possible sources

---

## 8. Notes on Comparison

1. **Transfer Reactions**: 
   - Spectroscopic factors are key for comparison
   - Angular distributions sensitive to optical potentials
   - Zero-range vs finite-range effects

2. **Elastic Scattering**:
   - Forward angle data most sensitive to optical potential
   - Backward angles sensitive to nuclear structure
   - Total cross section provides integral check

3. **Inelastic Scattering**:
   - Deformation parameters may need adjustment
   - Angular distributions test transition form factors
   - B(Eλ) values provide normalization

4. **Halo Nuclei**:
   - Binding energies are well-known (good check)
   - ANCs are key observables for transfer reactions
   - Decay lengths test asymptotic behavior

---

## 9. Format for Adding New Data

When adding new experimental data, use this format:

```markdown
| Parameter | Angle/Energy | Calculated | Experimental | Ratio | Reference | Agreement |
|-----------|--------------|------------|--------------|-------|-----------|-----------|
| dσ/dΩ (mb/sr) | θ = X° | Y.YY | Z.ZZ | Y.YY/Z.ZZ | Ref. | ✓/~/✗ |
```

Where:
- **Parameter**: Physical quantity (dσ/dΩ, E_b, S, etc.)
- **Angle/Energy**: Specific angle or energy point
- **Calculated**: Value from Nexus-NRS calculation
- **Experimental**: Value from literature/experiment
- **Ratio**: Calculated/Experimental
- **Reference**: Citation or source
- **Agreement**: ✓ (good, <10% difference), ~ (fair, 10-30%), ✗ (poor, >30%), — (not available)
