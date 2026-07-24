(ns dwba.o16-dp-radial-i-benchmark-test
  (:require [clojure.string :as str]
            [clojure.test :refer :all]
            [dwba.benchmark.o16-dp-radial-i-benchmark :as rib]
            [complex :as c]))

(deftest parse-ref-tsv-test
  (let [m (rib/parse-o16-dp-radial-i-ref-tsv
           "# c\nL_alpha L_beta ReI ImI\n0\t2\t1.5\t-0.25\n")]
    (is (= m {[0 2] {:re 1.5 :im -0.25}}))))

(deftest nexus-radial-i-tsv-roundtrip-smoke-test
  (let [rows (rib/o16-dp-radial-i-rows-benchmark :h 0.1 :r-max 12.0 :L-max 3 :attach-coulomb-sigma? true)
        s    (rib/o16-dp-radial-i-rows->tsv-string rows)]
    (is (pos? (count rows)))
    (is (str/includes? s "L_alpha"))
    (is (str/includes? s "Re_I_eff"))))

(deftest report-diff-self-consistent-test
  "Ref built from Nexus **I_eff** → relative magnitude error should be ~0."
  (let [rows (rib/o16-dp-radial-i-rows-benchmark :h 0.1 :r-max 12.0 :L-max 2 :attach-coulomb-sigma? true)
        ref (into {}
                  (map (fn [r]
                         (let [Ie (rib/o16-dp-row-I-eff r)]
                           [[(long (:L-alpha r)) (long (:L-beta r))]
                            {:re (c/re Ie) :im (c/im Ie)}]))
                  rows))
        res (rib/report-o16-dp-radial-i-diff rows ref :nexus-field :I-eff :print? false)]
    (is (pos? (:compared-pairs res)) "at least one (La,Lb) pair compared")
    (is (< (:max-rel-mag-err res) 1e-9) (pr-str res))))
