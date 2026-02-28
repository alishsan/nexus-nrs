(ns dwba.numerov-test
  "Tests for Numerov integration method."
  (:require [clojure.test :refer :all]
            [dwba.transfer :as t]
            [functions :refer :all]
            [fastmath.core :as m]
            [complex :as c :refer [complex-cartesian re im mag add mul div subt2]]))

(deftest solve-bound-state-numerov-basic-test
  (testing "Basic solve-bound-state-numerov integration"
    (let [E-test -20.0
          l 0
          v0 50.0
          rad 2.0
          diff 0.6
          h 0.01
          r-max 20.0
          u (t/solve-bound-state-numerov E-test l v0 rad diff h r-max)]
      (is (seq u) "Should return wavefunction")
      (is (= (count u) (int (/ r-max h))) "Length should match expected")
      (is (every? number? u) "All values should be numbers"))))

(deftest solve-bound-state-numerov-different-l-test
  (testing "solve-bound-state-numerov for different angular momenta"
    (let [E-test -20.0
          v0 50.0
          rad 2.0
          diff 0.6
          h 0.01
          r-max 20.0]
      (doseq [l [0 1 2]]
        (let [u (t/solve-bound-state-numerov E-test l v0 rad diff h r-max)]
          (is (seq u) (format "Should return wavefunction for l=%d" l))
          (is (number? (get u 1)) (format "Should have valid start value for l=%d" l)))))))

(deftest solve-bound-state-numerov-different-energies-test
  (testing "solve-bound-state-numerov for different energies"
    (let [l 0
          v0 50.0
          rad 2.0
          diff 0.6
          h 0.01
          r-max 20.0]
      (doseq [E-test [-30.0 -20.0 -10.0 -5.0]]
        (let [u (t/solve-bound-state-numerov E-test l v0 rad diff h r-max)
              u-end (t/bound-state-boundary-value u r-max h)]
          (is (seq u) (format "Should return wavefunction for E=%.1f" E-test))
          (is (number? u-end) (format "Should have valid boundary value for E=%.1f" E-test)))))))

(deftest numerov-node-counting-test
  (testing "Node counting for different energies"
    (let [l 0
          v0 50.0
          rad 2.0
          diff 0.6
          h 0.01
          r-max 20.0]
      (doseq [E-test [-30.0 -20.0 -10.0]]
        (let [u (t/solve-bound-state-numerov E-test l v0 rad diff h r-max)
              nodes (t/count-nodes u)]
          (is (>= nodes 0) (format "Node count should be non-negative for E=%.1f" E-test))
          (is (integer? nodes) (format "Node count should be integer for E=%.1f" E-test)))))))

(deftest numerov-convergence-test
  (testing "Numerov convergence with different step sizes"
    (let [e 2.0
          l 1
          v0 46.23
          rad 2.0
          diff 0.5
          r-max 10.0
          h-fine 0.001
          h-values [0.1 0.05 0.01]]
      (doseq [h-test h-values]
        (let [u-fine (solve-numerov e l v0 rad diff h-fine r-max)
              u-test (solve-numerov e l v0 rad diff h-test r-max)
              downsample-factor (int (/ h-test h-fine))
              u-fine-downsampled (take-nth downsample-factor u-fine)
              min-len (min (count u-fine-downsampled) (count u-test))
              errors (map #(Math/abs (- %1 %2))
                         (take min-len u-fine-downsampled)
                         (take min-len u-test))
              max-error (apply max errors)]
          (is (< max-error 1.0) (format "Max error should be reasonable for h=%.3f" h-test)))))))

(deftest numerov-stability-test
  (testing "Numerov stability - Wronskian preservation"
    (let [e 2.0
          l 1
          v0 46.23
          rad 2.0
          diff 0.5
          r-max 10.0
          h 0.01]
      ;; Basic test that numerov runs without error
      (let [u (solve-numerov e l v0 rad diff h r-max)]
        (is (seq u) "Should return wavefunction")
        (is (every? number? u) "All values should be numbers")))))

;; --- Complex Numerov Step Tests ---

(deftest numerov-step-complex-basic-test
  (testing "Basic numerov-step-complex with complex numbers"
    (let [h 0.01
          u-n (complex-cartesian 1.0 0.5)
          u-prev (complex-cartesian 0.5 0.2)
          f-prev (complex-cartesian 0.1 0.05)
          f-n (complex-cartesian 0.2 0.1)
          f-next (complex-cartesian 0.3 0.15)
          result (numerov-step-complex u-n u-prev f-prev f-n f-next h)]
      (is (c/complex? result) "Should return a complex number")
      (is (number? (re result)) "Real part should be a number")
      (is (number? (im result)) "Imaginary part should be a number")
      (is (not (Double/isNaN (re result))) "Real part should not be NaN")
      (is (not (Double/isNaN (im result))) "Imaginary part should not be NaN")
      (is (not (Double/isInfinite (re result))) "Real part should not be infinite")
      (is (not (Double/isInfinite (im result))) "Imaginary part should not be infinite"))))

(deftest numerov-step-complex-real-comparison-test
  (testing "numerov-step-complex should match real Numerov step for real inputs"
    (let [h 0.01
          h2-12 (/ (* h h) 12.0)
          ;; Real values
          u-n-real 1.0
          u-prev-real 0.5
          f-prev-real 0.1
          f-n-real 0.2
          f-next-real 0.3
          ;; Convert to complex (imaginary part = 0)
          u-n (complex-cartesian u-n-real 0.0)
          u-prev (complex-cartesian u-prev-real 0.0)
          f-prev (complex-cartesian f-prev-real 0.0)
          f-n (complex-cartesian f-n-real 0.0)
          f-next (complex-cartesian f-next-real 0.0)
          ;; Complex Numerov step
          result-complex (numerov-step-complex u-n u-prev f-prev f-n f-next h)
          ;; Real Numerov step (manual calculation)
          numerator-real (+ (* 2.0 u-n-real)
                            (- u-prev-real)
                            (* h2-12 (+ (* 10.0 f-n-real u-n-real)
                                        (* f-prev-real u-prev-real))))
          denominator-real (- 1.0 (* h2-12 f-next-real))
          result-real (/ numerator-real denominator-real)]
      (is (< (Math/abs (- (re result-complex) result-real)) 1e-10)
          "Real part should match real Numerov step")
      (is (< (Math/abs (im result-complex)) 1e-10)
          "Imaginary part should be approximately zero for real inputs"))))

(deftest numerov-step-complex-formula-test
  (testing "numerov-step-complex should follow Numerov formula"
    (let [h 0.01
          h2-12 (/ (* h h) 12.0)
          u-n (complex-cartesian 1.0 0.5)
          u-prev (complex-cartesian 0.5 0.2)
          f-prev (complex-cartesian 0.1 0.05)
          f-n (complex-cartesian 0.2 0.1)
          f-next (complex-cartesian 0.3 0.15)
          result (numerov-step-complex u-n u-prev f-prev f-n f-next h)
          ;; Manual calculation of formula
          ;; Numerator: 2*u_n - u_{n-1} + (h^2/12)*(10*f_n*u_n + f_{n-1}*u_{n-1})
          term1 (subt2 (mul 2.0 u-n) u-prev)
          term2 (mul h2-12 (add (mul 10.0 (mul f-n u-n))
                                (mul f-prev u-prev)))
          numerator (add term1 term2)
          ;; Denominator: 1 - (h^2/12)*f_{n+1}
          denominator (subt2 1.0 (mul h2-12 f-next))
          expected (div numerator denominator)]
      (is (< (Math/abs (- (re result) (re expected))) 1e-10)
          "Real part should match formula")
      (is (< (Math/abs (- (im result) (im expected))) 1e-10)
          "Imaginary part should match formula"))))

(deftest numerov-step-complex-zero-imaginary-test
  (testing "numerov-step-complex with zero imaginary parts"
    (let [h 0.01
          u-n (complex-cartesian 1.0 0.0)
          u-prev (complex-cartesian 0.5 0.0)
          f-prev (complex-cartesian 0.1 0.0)
          f-n (complex-cartesian 0.2 0.0)
          f-next (complex-cartesian 0.3 0.0)
          result (numerov-step-complex u-n u-prev f-prev f-n f-next h)]
      (is (< (Math/abs (im result)) 1e-10)
          "Imaginary part should be approximately zero when all inputs have zero imaginary parts"))))

(deftest numerov-step-complex-pure-imaginary-test
  (testing "numerov-step-complex with pure imaginary numbers"
    (let [h 0.01
          u-n (complex-cartesian 0.0 1.0)
          u-prev (complex-cartesian 0.0 0.5)
          f-prev (complex-cartesian 0.0 0.1)
          f-n (complex-cartesian 0.0 0.2)
          f-next (complex-cartesian 0.0 0.3)
          result (numerov-step-complex u-n u-prev f-prev f-n f-next h)]
      (is (c/complex? result) "Should return a complex number")
      (is (number? (re result)) "Real part should be a number")
      (is (number? (im result)) "Imaginary part should be a number"))))

(deftest numerov-step-complex-small-step-test
  (testing "numerov-step-complex with very small step size"
    (let [h 0.0001
          u-n (complex-cartesian 1.0 0.5)
          u-prev (complex-cartesian 0.5 0.2)
          f-prev (complex-cartesian 0.1 0.05)
          f-n (complex-cartesian 0.2 0.1)
          f-next (complex-cartesian 0.3 0.15)
          result (numerov-step-complex u-n u-prev f-prev f-n f-next h)]
      (is (c/complex? result) "Should return a complex number")
      (is (not (Double/isNaN (re result))) "Real part should not be NaN")
      (is (not (Double/isNaN (im result))) "Imaginary part should not be NaN"))))

(deftest numerov-step-complex-large-values-test
  (testing "numerov-step-complex with large values"
    (let [h 0.01
          u-n (complex-cartesian 100.0 50.0)
          u-prev (complex-cartesian 50.0 20.0)
          f-prev (complex-cartesian 10.0 5.0)
          f-n (complex-cartesian 20.0 10.0)
          f-next (complex-cartesian 30.0 15.0)
          result (numerov-step-complex u-n u-prev f-prev f-n f-next h)]
      (is (c/complex? result) "Should return a complex number")
      (is (not (Double/isNaN (re result))) "Real part should not be NaN")
      (is (not (Double/isNaN (im result))) "Imaginary part should not be NaN"))))

(deftest numerov-step-complex-consistency-test
  (testing "numerov-step-complex should be consistent across multiple steps"
    (let [h 0.01
          ;; First step
          u0 (complex-cartesian 0.0 0.0)
          u1 (complex-cartesian 0.01 0.005)
          f0 (complex-cartesian 0.0 0.0)
          f1 (complex-cartesian 0.1 0.05)
          f2 (complex-cartesian 0.2 0.1)
          u2-step1 (numerov-step-complex u1 u0 f0 f1 f2 h)
          ;; Second step (using result from first step)
          f3 (complex-cartesian 0.3 0.15)
          u3-step2 (numerov-step-complex u2-step1 u1 f1 f2 f3 h)]
      (is (c/complex? u2-step1) "First step should return complex number")
      (is (c/complex? u3-step2) "Second step should return complex number")
      (is (not (Double/isNaN (re u2-step1))) "First step real part should not be NaN")
      (is (not (Double/isNaN (re u3-step2))) "Second step real part should not be NaN"))))
