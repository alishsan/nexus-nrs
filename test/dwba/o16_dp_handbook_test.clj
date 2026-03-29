(ns dwba.o16-dp-handbook-test
  "Smoke tests for **¹⁶O(d,p)¹⁷O** handbook ZR pipeline (`dwba.benchmark.o16-dp-handbook`)."
  (:require [clojure.test :refer :all]
            [dwba.benchmark.o16-dp-handbook :as oh]))

(deftest o16-dp-handbook-kinematics-test
  (let [{:keys [e-cm-f Q-mev k-i k-f]} (oh/o16-dp-kinematics)]
    (is (pos? e-cm-f))
    (is (> Q-mev 0.0) "d,p on light target: Q should be exothermic (≈ +2–3 MeV)")
    (is (pos? k-i))
    (is (pos? k-f))))

(deftest o16-dp-handbook-dsigma-finite-test
  (let [s (oh/o16-dp-dsigma-handbook-mb-sr 35.0 :h 0.1 :r-max 16.0 :L-max 4)]
    (is (Double/isFinite (double s)))
    (is (pos? (double s)))))

(deftest o16-dp-handbook-curve-nonempty-test
  (let [curve (oh/o16-dp-angular-curve-handbook-mb-sr [15.0 75.0] :h 0.1 :r-max 16.0 :L-max 4)]
    (is (= 2 (count curve)))
    (is (every? #(Double/isFinite (:differential_cross_section_mb_sr %)) curve))))
