(ns dwba.ca40-pd-austern-test
  "Smoke tests for **⁴⁰Ca(p,d)** **(5.6)** angular DWBA (`dwba.benchmark.ca40-pd-austern`)."
  (:require [clojure.test :refer :all]
            [dwba.benchmark.ca40-pd-austern :as pd]))

(deftest ca40-pd-kinematics-q-value-open-exit
  (let [{:keys [e-cm-f Q-mev k-i k-f]} (pd/ca40-pd-kinematics 18.0)]
    (is (< (Math/abs (+ Q-mev 5.847)) 0.02) "crude Q ≈ −5.85 MeV")
    (is (> e-cm-f 0.5))
    (is (pos? k-i))
    (is (pos? k-f))))

(deftest ca40-pd-dsigma-austern-eq-56-finite-test
  (let [s (pd/ca40-pd-dsigma-austern-eq-56-mb-sr 35.0 :h 0.1 :r-max 16.0 :L-max 4)]
    (is (Double/isFinite (double s)))
    (is (pos? (double s))
        "non-zero: T_m from D₀√(2ℓ+1)β_m and spin weight like ca40-dp (no broken 6j prefactor)")))

(deftest ca40-pd-angular-curve-nonempty-test
  (let [curve (pd/ca40-pd-angular-curve-austern-eq-56-mb-sr [10.0 90.0] :h 0.1 :r-max 16.0 :L-max 4)]
    (is (= 2 (count curve)))
    (is (every? #(Double/isFinite (:differential_cross_section_mb_sr %)) curve))))

(deftest ca40-pd-L-alpha-only-slice-test
  "**`:L-alpha-only`** keeps one entrance partial wave; σ is **≤** full sum when that slice is non-empty."
  (let [opts [:h 0.1 :r-max 16.0 :L-max 6]
        s35-full (apply pd/ca40-pd-dsigma-austern-eq-56-mb-sr 35.0 opts)
        s35-L3 (apply pd/ca40-pd-dsigma-austern-eq-56-mb-sr 35.0 :L-alpha-only 3 opts)
        curve (pd/ca40-pd-angular-curve-austern-eq-56-mb-sr [20.0] :L-alpha-only 5 :h 0.1 :r-max 16.0 :L-max 6)]
    (is (pos? (double s35-full)))
    (is (Double/isFinite (double s35-L3)))
    (is (<= (double s35-L3) (* 1.0001 (double s35-full))))
    (is (every? #(Double/isFinite (:differential_cross_section_mb_sr %)) curve))))
