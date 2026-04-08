(ns dwba.alpha-n-elastic-neutral-test
  "**α + n**, **Z₁Z₂ = 0**: Numerov **u** + R-matrix **S** must match **`functions/s-matrix`** per partial wave."
  (:require [clojure.test :refer [deftest is testing]]
            [complex :as c]
            [dwba.benchmark.alpha-n-elastic-neutral :as an]))

(defn- cabs2 [z]
  (+ (Math/pow (double (c/re z)) 2.0)
     (Math/pow (double (c/im z)) 2.0)))

(defn- S-distance [s-ref s-num]
  (Math/sqrt (cabs2 (c/subt s-ref s-num))))

(deftest alpha-n-neutral-S-numerov-matches-s-matrix-test
  (let [e-cm 10.0
        ;; Illustrative shallow real WS (same style as other tests); **no** imaginary part.
        ws [22.0 3.4 0.62]
        r-max 50.0
        ;; Finer **h** tightens **R** at **a** vs coarse **`r-matrix`** Euler (**dr = 0.001**).
        h 0.01
        L-max 10]
    (testing "per-L S: Numerov R-match vs functions/s-matrix (η=0)"
      (doseq [L (range 0 (inc L-max))
              :let [S-ref (an/alpha-n-sigma-L-ref e-cm ws (long L))
                    S-num (an/alpha-n-S-from-numerov-for-L e-cm ws (long L) r-max h)
                    d (S-distance S-ref S-num)]]
        (is (< d 2.5e-2)
            (format "L=%d |S_ref - S_num|=%.4e ref=%s num=%s"
                    L d (pr-str S-ref) (pr-str S-num)))))))
