(ns dwba.form-factors-test
  "Tests for transfer form factor calculations."
  (:require [clojure.test :refer :all]
            [dwba.form-factors :as ff]
            [dwba.transfer :as t]
            [functions :refer :all]
            [fastmath.core :as m]))

(def ws-params [50.0 2.0 0.6])  ; V0=50 MeV, R0=2.0 fm, a0=0.6 fm
(def r-max 20.0)
(def h 0.01)
(def ^:private no-so {:no-spin-orbit true})

(defn- R-from-reduced-u
  " **`ff/overlap-integral`** expects radial **R(r)**; **`solve-bound-state`** returns reduced **u = rR**."
  [u h]
  (mapv (fn [^long i]
          (let [r (* (double i) h)
                uu (double (get u i))]
            (if (< r 1e-14)
              0.0
              (/ uu r))))
        (range (count u))))

(deftest form-factor-at-r-basic-test
  (testing "form-factor-at-r calculates form factor at specific r"
    (let [phi-i [0.0 0.1 0.2 0.3]
          phi-f [0.0 0.2 0.4 0.6]
          r 0.02
          result (ff/form-factor-r r phi-i phi-f h)]
      (is (number? result) "Should return a number")
      (is (>= result 0.0) "Form factor should be non-negative for real wavefunctions"))))

(deftest overlap-integral-same-state-test
  (testing "Overlap integral of state with itself gives norm squared"
    (let [result-1s (t/solve-bound-state ws-params 1 0 nil r-max h no-so)
          phi-1s (R-from-reduced-u (:normalized-wavefunction result-1s) h)]
      (when (seq phi-1s)
        (let [norm-squared (ff/overlap-integral phi-1s phi-1s r-max h)]
          (is (< (Math/abs (- norm-squared 1.0)) 0.1)
              (format "Norm squared should be close to 1.0 for normalized state: got %.6f" norm-squared)))))))

(deftest overlap-integral-orthogonal-states-test
  (testing "Overlap integral between orthogonal states should be small"
    ;; Use 1p (not 2s): excited s search can return the same radial root as 1s.
    (let [result-1s (t/solve-bound-state ws-params 1 0 nil r-max h no-so)
          result-1p (t/solve-bound-state ws-params 1 1 nil r-max h no-so)]
      (when (and (seq (:normalized-wavefunction result-1s))
                (seq (:normalized-wavefunction result-1p)))
        (let [phi-1s (R-from-reduced-u (:normalized-wavefunction result-1s) h)
              phi-1p (R-from-reduced-u (:normalized-wavefunction result-1p) h)
              overlap (ff/overlap-integral phi-1s phi-1p r-max h)]
          (is (< (Math/abs overlap) 0.5)
              (format "Overlap integral between orthogonal states should be small: got %.6f" overlap)))))))

(deftest form-factor-function-test
  (testing "form-factor-function returns vector of form factor values"
    (let [result-1s (t/solve-bound-state ws-params 1 0 nil r-max h no-so)
          phi-1s (R-from-reduced-u (:normalized-wavefunction result-1s) h)]
      (when (seq phi-1s)
        (let [distribution (ff/form-factor-function phi-1s phi-1s h)]
          (is (seq distribution) "Should return non-empty vector")
          (is (= (count distribution) (count phi-1s)) "Should have same length as wavefunction")
          (is (every? number? distribution) "All values should be numbers"))))))

(deftest normalized-overlap-test
  (testing "normalized-overlap gives normalized overlap coefficient"
    (let [result-1s (t/solve-bound-state ws-params 1 0 nil r-max h no-so)
          phi-1s (R-from-reduced-u (:normalized-wavefunction result-1s) h)]
      (when (seq phi-1s)
        (let [normalized-overlap-val (ff/normalized-overlap phi-1s phi-1s r-max h)]
          (is (number? normalized-overlap-val) "Should return a number")
          (is (<= (Math/abs normalized-overlap-val) (+ 1.0 1e-6))
              "Normalized overlap should be ≤ 1.0")
          ;; State with itself should give close to 1.0
          (is (> normalized-overlap-val 0.9)
              (format "Normalized overlap of state with itself should be close to 1.0: got %.6f" normalized-overlap-val))))))

(deftest momentum-space-overlap-test
  (testing "momentum-space-overlap calculates momentum-space overlap integral"
    (let [result-1s (t/solve-bound-state ws-params 1 0 nil r-max h no-so)
          phi-1s (R-from-reduced-u (:normalized-wavefunction result-1s) h)
          q 0.5]  ; Momentum transfer in fm⁻¹
      (when (seq phi-1s)
        (let [momentum-overlap (ff/momentum-space-overlap phi-1s phi-1s r-max h q)]
          (is (number? momentum-overlap) "Should return a number")
          ;; At q=0, should reduce to regular overlap integral
          (let [q0-overlap (ff/momentum-space-overlap phi-1s phi-1s r-max h 0.0)
                regular-overlap (ff/overlap-integral phi-1s phi-1s r-max h)]
            (is (< (Math/abs (- q0-overlap regular-overlap)) 0.1)
                (format "At q=0, momentum-space overlap should equal regular overlap integral"))))))))

(deftest overlap-integral-symmetry-test
  (testing "Overlap integral should satisfy O(φ_i, φ_f) = O*(φ_f, φ_i) for complex wavefunctions"
    (let [result-1s (t/solve-bound-state ws-params 1 0 nil r-max h no-so)
          result-1p (t/solve-bound-state ws-params 1 1 nil r-max h no-so)]
      (when (and (seq (:normalized-wavefunction result-1s))
                (seq (:normalized-wavefunction result-1p)))
        (let [phi-1s (R-from-reduced-u (:normalized-wavefunction result-1s) h)
              phi-1p (R-from-reduced-u (:normalized-wavefunction result-1p) h)
              overlap-12 (ff/overlap-integral phi-1s phi-1p r-max h)
              overlap-21 (ff/overlap-integral phi-1p phi-1s r-max h)]
          ;; For real wavefunctions, should be equal
          (is (< (Math/abs (- overlap-12 overlap-21)) 0.01)
              (format "Overlap integral should be symmetric: O(1s,1p)=%.6f, O(1p,1s)=%.6f" overlap-12 overlap-21)))))))

(deftest overlap-integral-convergence-test
  (testing "Overlap integral should converge with smaller step size"
    (let [result-1s (t/solve-bound-state ws-params 1 0 nil r-max h no-so)
          phi-1s (R-from-reduced-u (:normalized-wavefunction result-1s) h)]
      (when (seq phi-1s)
        (let [overlap-h1 (ff/overlap-integral phi-1s phi-1s r-max h)
              h2 (* h 2.0)
              ;; Resample wavefunction for h2 (simple downsampling)
              phi-1s-h2 (R-from-reduced-u (take-nth 2 (:normalized-wavefunction result-1s)) h2)
              overlap-h2 (ff/overlap-integral phi-1s-h2 phi-1s-h2 r-max h2)]
          ;; Results should be reasonably close
          (is (< (Math/abs (- overlap-h1 overlap-h2)) 0.1)
              (format "Overlap integral should be reasonably stable with step size: h=%.3f gives %.6f, h=%.3f gives %.6f"
                     h overlap-h1 h2 overlap-h2))))))))

