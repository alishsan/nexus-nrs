(ns dwba.austern-eq56-and-rutherford-test
  "**Elastic Rutherford** vs **Austern (5.6)** — different observables.

  - **Rutherford** is **elastic** point-Coulomb **dσ/dΩ** (non-relativistic CM):
    **((Z₁Z₂e²)/(4E_cm))² / sin⁴(θ/2)** (**mb/sr** here with **×10** via **`transfer-rutherford-dsigma-mb-sr`**).
  - **Austern (5.6)** (**`austern-reduced-amplitude-beta-sum-eq-5-6`**) builds the **transfer** DWBA
    reduced amplitude **β_{sj}^{ℓm}** from **I_{L_β L_α}**, **e^{i(σ_{αL_α}+σ_{βL_β})}**, Clebsch–Gordan
    factors, and **Y_{L_β}^m(Θ,0)**. Turning **nuclear OM** off in χ still leaves a **rearrangement**
    amplitude (often vanishing if overlaps / coupling go to zero) — it does **not** collapse **(5.6)** to
    **elastic** Rutherford.

  **Elastic** check belongs in **`functions/differential-cross-section`** (partial-wave **S** with Coulomb
  only), not in the **(5.5)+(5.6)** pickup pipeline."
  (:require [clojure.test :refer :all]
            [complex :as c]
            [dwba.transfer :as t]))

(defn- approx-eq? [^double a ^double b ^double tol]
  (< (Math/abs (- a b)) tol))

(deftest transfer-rutherford-analytic-90deg-matches-helper
  (let [th (/ Math/PI 2.0)
        ;; Z₁Z₂e² = 28.8 MeV·fm, E_cm = 18 MeV → (28.8/72)² / sin⁴(45°) × 10 = 6.4 mb/sr
        expected 6.4
        code (t/transfer-rutherford-dsigma-mb-sr 1 20 18.0 th)]
    (is (approx-eq? (double code) expected 1e-6))))

(deftest austern-eq-5-6-beta-is-not-elastic-rutherford-cross-section
  (let [rows [{:L-alpha 0 :L-beta 0 :sigma-alpha 0.0 :sigma-beta 0.0 :I 1.0}]
        th (/ Math/PI 2.0)
        beta (t/austern-reduced-amplitude-beta-sum-eq-5-6 0 0 th rows)
        bmag (double (c/mag beta))
        y00 (/ 1.0 (Math/sqrt (* 4.0 Math/PI)))
        ruth (t/transfer-rutherford-dsigma-mb-sr 1 20 18.0 th)]
    (is (approx-eq? bmag y00 1e-9)
        "single L_α=L_β=ℓ=m=0 row, σ=0, I=1 ⇒ β ∝ Y_00")
    ;; |β|² is ~0.08; Rutherford at 90° is 6.4 mb/sr — not the same object or units pathway
    (is (not (approx-eq? (* bmag bmag) ruth 0.5)))))
