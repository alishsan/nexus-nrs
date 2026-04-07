(ns dwba.li11-pp-paper1708-test
  (:require [clojure.test :refer :all]
            [dwba.benchmark.li11-pp-paper1708 :as li11]))

(deftest li11-paper-dsigma-positive
  (testing "paper optical DWBA returns finite positive mb/sr (coarse grid)"
    (let [s (li11/li11-pp-dsigma-paper-mb-sr 90.0 :optical-set :V :L-max 6 :r-max 20 :h 0.05)]
      (is (pos? s))
      (is (Double/isFinite (double s))))))

(deftest li11-paper-kinematics-exit-energy
  (testing "E_cm - E_x positive for default ex 0.8 MeV"
    (let [k (:e-cm-f (li11/li11-pp-kinematics 6.0 0.8))]
      (is (> k 4.0)))))
