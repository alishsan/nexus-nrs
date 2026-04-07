# Paper Materials

This folder contains materials for the Nexus-NRS paper, including benchmark reactions and experimental comparisons.

## LaTeX manuscript (MDPI)

- **`MDPI_template_ACS.zip`**: official MDPI ACS-style template; extract with `unzip -o MDPI_template_ACS.zip` so `Definitions/mdpi.cls` and **`template.tex`** exist beside **`mdpi_proceedings_manuscript.tex`**.
- **`mdpi_proceedings_manuscript.tex`**: full article using `\documentclass[particles,proceedingpaper,submit,pdftex,moreauthors]{Definitions/mdpi}`. Build from this folder: `pdflatex mdpi_proceedings_manuscript.tex` (run twice for references/LastPage). Place **`nexus_nrs_dashboard.png`** here for the figure.
- Before submission, switch **`submit`** → **`accept`** in `\documentclass[...]` per MDPI instructions; fill `\datereceived`, `\dateaccepted`, etc.
- **Graphical abstract (website / TOC):** `graphical_abstract.png` — single self-explanatory image; upload where the submission form asks for it (do **not** place it inside the manuscript-only ZIP unless the instructions say so). Typical specs: raster **PNG** or **TIFF**, **≥600 px** wide, readable at thumbnail size; edit the PNG if you want journal-specific colors or branding.
- **Zip for MDPI “Manuscript (Word/ZIP)” upload (main text only):** `Nexus-NRS_HALO40_manuscript_MDPI.zip` — contains `mdpi_proceedings_manuscript.tex`, main-text figure `nexus_nrs_dashboard.png`, and `Definitions/` (class files required to compile). **No PDF, no supplementary files** inside the zip (upload supplementary materials separately if the journal asks). Regenerate:  
  `cd paper && zip -r Nexus-NRS_HALO40_manuscript_MDPI.zip mdpi_proceedings_manuscript.tex nexus_nrs_dashboard.png Definitions -x "*.DS_Store"`

## Files

1. **`benchmark_reactions.md`** - Complete list of benchmark reactions with parameters, nuclear properties, and key results
2. **`experimental_comparison_table.md`** - Comparison table between calculated and experimental values
3. **`run_benchmarks.clj`** - Script to run all benchmark calculations and generate plots
4. **`run_benchmarks.sh`** - Shell script wrapper for easy execution
5. **`README.md`** - This file

## Quick Reference

### Benchmark Reactions

1. **Transfer Reactions**:
   - 16O(p,d)15O (neutron pickup)
   - ¹⁰Be(d,p)¹¹Be (transfer to halo nucleus)

2. **Elastic Scattering**:
   - 11Li(d,d) (elastic on halo nucleus)

3. **Inelastic Scattering**:
   - 11Li(d,d') (monopole L=0 excitation)
   - ¹²C(α,α')¹²C* (quadrupole excitation)

4. **Halo Nuclei**:
   - ¹¹Be (neutron halo, E_b = 504 keV)
   - ⁸B (proton halo, E_b = 137 keV)

### Status

- ✅ Benchmark reactions documented
- ✅ Experimental comparison table created
- ⏳ Calculations need to be run to fill in values
- ⏳ Experimental data needs to be collected from literature

## Next Steps

1. Run calculations for all benchmark reactions
2. Collect experimental data from literature
3. Fill in comparison table
4. Calculate agreement metrics (χ², ratios, etc.)
5. Document any discrepancies

## Usage

### Running Benchmark Calculations

To run all benchmark calculations and generate plots:

**Option 1: Using the shell script**
```bash
./paper/run_benchmarks.sh
```

**Option 2: Using Clojure directly**
```bash
clj -M -e "(load-file \"paper/run_benchmarks.clj\")" -e "(paper.run-benchmarks/run-all-benchmarks)"
```

**Option 3: From REPL (Recommended)**
```clojure
;; First, ensure you're in the project root directory
;; Then load the helper script:
(load-file "paper/load_benchmarks.clj")

;; Or load directly with absolute path:
(load-file (str (System/getProperty "user.dir") "/paper/run_benchmarks.clj"))

;; Then run:
(paper.run-benchmarks/run-all-benchmarks)
```

**Option 4: From REPL (Alternative - if path issues)**
```clojure
;; Change to project root first
(System/setProperty "user.dir" "/path/to/nexus-nrs")
(load-file "paper/run_benchmarks.clj")
(paper.run-benchmarks/run-all-benchmarks)
```

### Output

The script generates:
- Individual plots for each reaction in `paper/plots/`
- Summary plot with all reactions in `paper/plots/summary.png`
- Results data in `paper/benchmark_results.edn` (EDN format)

### Adding New Benchmark Reactions

1. Add to `benchmark_reactions.md`
2. Add corresponding row(s) to `experimental_comparison_table.md`
3. Add calculation function to `run_benchmarks.clj`
4. Run calculation and fill in values
5. Find experimental data and add to table
