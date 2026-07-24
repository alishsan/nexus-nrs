(ns dwba.resonance-transfer-test
  "Sanity tests for the 'quasi-bound resonance' DWBA transfer pipeline
   (`dwba.resonance-transfer`) — ¹¹Li(p,d)¹⁰Li* (¹⁰Li 1p₁/₂ resonance, E_r ≈ 0.62 MeV).

   These are SANITY checks (finiteness, sign, plausible order of magnitude), not
   regression checks against a fitted reference curve -- the physics here is a genuine
   first-principles estimate (see namespace docstring in
   `src/dwba/resonance_transfer.clj` for the caveats)."
  (:require [clojure.test :refer :all]
            [dwba.resonance-transfer :as rt]))

(deftest resonance-search-converges-test
  (let [search (rt/find-resonance-v0)]
    (is (:converged? search) (pr-str search))
    (is (Double/isFinite (double (:v0 search))))
    (is (< 5.0 (:v0 search) 100.0) "well depth should be a plausible nuclear WS depth")))

(deftest resonance-width-plausible-test
  "Wigner Gamma should be finite, positive, and within an order of magnitude of the
   measured Gamma = 0.33 +/- 0.07 MeV (Sanetullaev et al. 2016) -- NOT required to match
   closely, since a single real-potential well depth cannot independently fit both E_r
   and Gamma."
  (let [search (rt/find-resonance-v0)
        gamma (rt/resonance-width-wigner :v0 (:v0 search))]
    (is (Double/isFinite gamma))
    (is (pos? gamma))
    (is (< 0.05 gamma 3.0)
        (format "Gamma_calc=%.4f MeV should be within a broad plausible range around Gamma_lit=%.2f MeV"
                gamma rt/Gamma-lit-10Li))))

(deftest quasi-bound-wavefunction-normalized-test
  (let [search (rt/find-resonance-v0)
        wf (rt/quasi-bound-resonance-wavefunction :v0 (:v0 search))
        u (:u wf)
        h (:h wf)
        norm (* h (reduce + (map #(* (double %) (double %)) u)))]
    (is (pos? (:r-cut wf)))
    (is (> (count u) 10))
    (is (< 0.5 norm 2.0) (format "Simpson-normalized integral should be ~1 (got %.4f)" norm))))

(deftest kinematics-sane-test
  (let [kin (rt/li11-pd-kinematics)]
    (is (pos? (:E-cm-i kin)))
    (is (pos? (:E-cm-f kin)))
    (is (pos? (:k-i kin)))
    (is (pos? (:k-f kin)))
    (is (pos? (:E-lab-p-equiv kin)))
    (is (pos? (:E-lab-d-equiv-per-nucleon kin)))))

(deftest dsigma-finite-positive-test
  "End-to-end smoke test: dsigma/dOmega must be finite and positive at a representative
   angle, on a coarser grid than the example script (for test speed)."
  (let [search (rt/find-resonance-v0 :h 0.05 :r-max 30.0)
        wf (rt/quasi-bound-resonance-wavefunction :v0 (:v0 search) :h 0.05 :r-max 30.0)
        d (rt/li11-pd-dsigma-mb-sr 39.6 :resonance wf :h 0.05 :r-max 40.0 :L-max 2)]
    (is (Double/isFinite d))
    (is (pos? d))))

(deftest angular-comparison-all-finite-positive-test
  "All 7 EXFOR angles: dsigma finite/positive, ratio finite/positive, and the calc/exp
   ratio stays within a broad sanity band (no fitted scale factor is applied -- see
   namespace docstring; this only guards against gross unit/normalization errors,
   e.g. off by 1e6 or negative)."
  (let [search (rt/find-resonance-v0 :h 0.05 :r-max 30.0)
        wf (rt/quasi-bound-resonance-wavefunction :v0 (:v0 search) :h 0.05 :r-max 30.0)
        cmp (rt/li11-pd-angular-comparison :resonance wf :h 0.05 :r-max 40.0 :L-max 2)]
    (is (= 7 (count cmp)))
    (doseq [{:keys [dsigma-mb-sr exp-mb-sr ratio]} cmp]
      (is (Double/isFinite dsigma-mb-sr))
      (is (pos? dsigma-mb-sr))
      (is (pos? exp-mb-sr))
      (is (Double/isFinite ratio))
      (is (< 1.0e-4 ratio 1.0e2)
          (format "ratio=%.4g should be within a broad sanity band (no fit applied)" ratio)))))
