(ns dwba.ca40-pd-handbook-test
  "Smoke tests for **⁴⁰Ca(p,d)** handbook ZR pipeline (`dwba.benchmark.ca40-pd-handbook`)."
  (:require [clojure.test :refer :all]
            [dwba.benchmark.ca40-pd-handbook :as ph]))

(deftest ca40-pd-handbook-kinematics-test
  (let [{:keys [e-cm-f Q-mev k-i k-f]} (ph/ca40-pd-kinematics 18.0)]
    (is (< (Math/abs (+ Q-mev 5.847)) 0.02))
    (is (> e-cm-f 0.5))
    (is (pos? k-i))
    (is (pos? k-f))))

(deftest ca40-pd-handbook-dsigma-finite-test
  (let [s (ph/ca40-pd-dsigma-handbook-mb-sr 35.0 :h 0.1 :r-max 16.0 :L-max 4)]
    (is (Double/isFinite (double s)))
    (is (pos? (double s)))))

(deftest ca40-pd-handbook-curve-nonempty-test
  (let [curve (ph/ca40-pd-angular-curve-handbook-mb-sr [15.0 75.0] :h 0.1 :r-max 16.0 :L-max 4)]
    (is (= 2 (count curve)))
    (is (every? #(Double/isFinite (:differential_cross_section_mb_sr %)) curve))))
