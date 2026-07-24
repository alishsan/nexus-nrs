(ns dwba.rmatrix-lagrange-test
  (:require [clojure.test :refer :all]
            [functions :as f]
            [functions.rmatrix-lagrange :as lm]))

(deftest gauss-legendre-01-test
  (testing "Gauss–Legendre nodes integrate polynomials on (0,1)"
    (let [{:keys [nodes weights]} (lm/gauss-legendre-01 8)
          quad (fn [g] (reduce + 0.0 (map (fn [x w] (* w (g x))) nodes weights)))]
      (is (= (count nodes) 8))
      (is (< (Math/abs (- (quad (constantly 1.0)) 1.0)) 1e-10))
      (is (< (Math/abs (- (quad #(* % %)) (/ 1.0 3.0))) 1e-8)))))

(deftest lagrange-kinetic-matrix-symmetric-test
  (testing "Kinetic+Bloch matrix is symmetric"
    (let [T (lm/lagrange-kinetic-bloch-matrix 12 8.0)]
      (doseq [i (range 12) j (range 12)]
        (is (< (Math/abs (- (get-in T [i j]) (get-in T [j i]))) 1e-10))))))

(deftest lagrange-s-matrix-with-numerov-Ra-test
  (testing "Lagrange S-matrix with Numerov R a matches phase-shift0"
    (let [E 2.0 L 1 V [46.23 2.0 0.5] a 10.0
          delta-lm (lm/phase-shift-lagrange-numerov-Ra E V a L)
          delta0 (f/phase-shift0 E V a L)]
      (is (< (Math/abs (- delta-lm delta0)) 1e-6)))))

(deftest lagrange-vs-numerov-numerov-Ra-test
  (testing "Numerov phase shift matches r-matrix-a S-matrix pipeline"
    (let [{:keys [numerov numerov-r-matrix-a]}
          (lm/lagrange-vs-numerov-phase-shift 2.0 1 [46.23 2.0 0.5] 10.0 15)]
      (is (< (Math/abs (- numerov numerov-r-matrix-a)) 0.05)))))

(deftest phase-shift-lagrange-star-numerov-method-test
  (testing "phase-shift-lagrange* :numerov-r-matrix-a"
    (let [d1 (lm/phase-shift-lagrange* 5.0 [40.0 2.89 0.65] 8.0 0 15 :numerov-r-matrix-a)
          d2 (f/phase-shift0 5.0 [40.0 2.89 0.65] 8.0 0)]
      (is (< (Math/abs (- d1 d2)) 1e-6)))))
