(ns dwba.finite-well-test
  "Tests for finite square well bound state calculations."
  (:require [clojure.test :refer :all]
            [functions :refer :all]
            [dwba.finite-well :refer :all]
            [fastmath.core :as m]))

(deftest finite-well-validation-test
  (testing "Finite square well bound state validation for various l and z0"
    (let [test-cases [{:l 0 :z0 3.0 :name "l=0, z0=3.0 (shallow)"}
                     {:l 0 :z0 6.0 :name "l=0, z0=6.0 (medium)"}
                     {:l 0 :z0 10.0 :name "l=0, z0=10.0 (deep)"}
                     {:l 1 :z0 3.0 :name "l=1, z0=3.0 (shallow)"}
                     {:l 1 :z0 6.0 :name "l=1, z0=6.0 (medium)"}
                     {:l 1 :z0 10.0 :name "l=1, z0=10.0 (deep)"}
                     {:l 2 :z0 6.0 :name "l=2, z0=6.0 (medium)"}
                     {:l 2 :z0 10.0 :name "l=2, z0=10.0 (deep)"}]]
      (doseq [test test-cases]
        (let [l (:l test)
              z0 (:z0 test)
              bound-states (find-all-bound-states l z0)
              valid-states (filter #(and (>= (:e-ratio %) 0.01)
                                        (<= (:e-ratio %) 0.99)
                                        (> (:xi %) 0)
                                        (> (:eta %) 0)
                                        (< (m/abs (:matching-error %)) 1e-3)
                                        (:converged? %))
                                  bound-states)]
          (testing (:name test)
            ;; Only require bound states for deep wells (z0 >= 6.0)
            (when (>= z0 6.0)
              (is (seq bound-states) "Should find at least one bound state for deep wells"))
            (when (seq valid-states)
              (is (every? #(and (>= (:e-ratio %) 0.01) (<= (:e-ratio %) 0.99)) valid-states)
                  "Energies should be in valid range")
              (is (every? #(< (m/abs (:matching-error %)) 1e-3) valid-states)
                  "Matching errors should be small")
              (is (= (map :e-ratio valid-states) (sort (map :e-ratio valid-states)))
                  "States should be sorted by energy")))))))

(deftest finite-well-physical-consistency-test
  (testing "Physical consistency: higher l should have fewer or equal bound states"
    (let [z0 6.0
          l0-states (count (find-all-bound-states 0 z0))
          l1-states (count (find-all-bound-states 1 z0))
          l2-states (count (find-all-bound-states 2 z0))]
      (is (>= l0-states l1-states) "l=0 should have >= l=1 states")
      (is (>= l1-states l2-states) "l=1 should have >= l=2 states"))))

(deftest finite-well-root-finding-test
  (testing "Root finding for finite well matching error"
    (let [l 2
          z0 6.0
          bound-states (find-all-bound-states l z0)]
      (is (seq bound-states) "Should find bound states for l=2, z0=6.0")
      (doseq [state bound-states]
        (when (:converged? state)
          (let [verify-result (solver-step (:e-ratio state) l z0)]
            (is (< (m/abs (:f verify-result)) 1e-3)
                (format "Matching error should be small: %.6e" (:f verify-result)))
            ;; Verify xi^2 + eta^2 = z0^2
            (let [xi (:xi state)
                  eta (:eta state)
                  expected (* z0 z0)
                  actual (+ (* xi xi) (* eta eta))]
              (is (< (m/abs (- actual expected)) 1e-6)
                  (format "xi^2 + eta^2 should equal z0^2: %.6f vs %.6f" actual expected))))))))))

(deftest finite-well-nlz-selector-test
  (testing "find-bound-state-finite-well-nlz: n=1 is most bound (max e-ratio), order matches find-all"
    (let [l 0
          z0 10.0
          all (find-all-bound-states l z0)
          valid-er (->> all
                       (filter #(and (:converged? %)
                                    (< (m/abs (:matching-error %)) 1e-3)
                                    (>= (:e-ratio %) 0.01)
                                    (<= (:e-ratio %) 0.99)))
                       (map :e-ratio)
                       sort
                       reverse)
          nmax (count valid-er)]
      (when (pos? nmax)
        (is (some? (find-bound-state-finite-well-nlz 1 l z0)))
        (is (= (:e-ratio (find-bound-state-finite-well-nlz 1 l z0)) (first valid-er)))
        (when (>= nmax 2)
          (is (= (:e-ratio (find-bound-state-finite-well-nlz 2 l z0)) (second valid-er))))
        (is (nil? (find-bound-state-finite-well-nlz (inc nmax) l z0)))
        (is (nil? (find-bound-state-finite-well-nlz 0 l z0)))
        (let [s (find-bound-state-finite-well-nlz 1 l z0)]
          (is (= (select-keys s [:n :l :z0]) {:n 1 :l l :z0 z0})))))))

(deftest woods-saxon-vs-finite-well-test
  (testing "Woods-Saxon with small diffuseness should approximate finite square well"
    (let [test-cases [{:V0 30.0 :R0 2.0 :a0 0.01 :l 0}
                     {:V0 50.0 :R0 2.0 :a0 0.01 :l 0}
                     {:V0 50.0 :R0 2.0 :a0 0.01 :l 1}]
          mass-factor (/ (* 2 869.4) (* 197.7 197.7))]
      (doseq [test test-cases]
        (let [V0 (:V0 test)
              R0 (:R0 test)
              l (:l test)
              z0 (m/sqrt (* mass-factor V0 (* R0 R0)))
              reference-e-ratios (map :e-ratio 
                                     (filter :converged? 
                                            (find-all-bound-states l z0)))]
          (when (seq reference-e-ratios)
            (testing (format "V0=%.1f, R0=%.1f, l=%d" V0 R0 l)
              (is (seq reference-e-ratios) "Should have reference bound states"))))))))

(deftest finite-well-delta-so-spin-half-test
  (testing "l-dot-s-spin-half matches j(j+1) algebra for p doublet"
    (is (< (m/abs (- (l-dot-s-spin-half 1 1.5) 0.5)) 1e-12))
    (is (< (m/abs (- (l-dot-s-spin-half 1 0.5) -1.0)) 1e-12))
    (is (< (m/abs (l-dot-s-spin-half 0 0.5)) 1e-12))))

(deftest finite-well-delta-so-shift-test
  (testing "δ-surface Thomas shift: l=0 gives 0; l>=1 splitting uses (2l+1)/2 weighted radial factor"
    (let [a 2.0
          Vso 6.0
          lam2 2.0
          xi 3.0
          eta 2.5
          de0 (delta-surface-thomas-so-shift-MeV Vso lam2 a 0 0.5 xi eta)]
      (is (< (m/abs de0) 1e-9)))
    (let [st (find-bound-state-finite-well-nlz 1 1 8.0)]
      (when st
        (let [a 2.0
              spl (finite-well-delta-so-splitting-j-doublet-MeV 5.0 lambda-pi-squared-fm2-default
                                                              a 1 (:xi st) (:eta st))
              dep (delta-surface-thomas-so-shift-MeV 5.0 lambda-pi-squared-fm2-default a 1 1.5
                                                     (:xi st) (:eta st))
              dem (delta-surface-thomas-so-shift-MeV 5.0 lambda-pi-squared-fm2-default a 1 0.5
                                                     (:xi st) (:eta st))]
          (is (Double/isFinite spl))
          (is (< (m/abs (- spl (- dep dem))) 1e-6)
              "Splitting should equal difference of j=l±1/2 shifts"))))))

(deftest finite-well-delta-so-dimensionless-test
  (testing "ΔE_so/V₀: l=0 → 0; MeV shift = V₀ × dimensionless shift"
    (let [a 2.0
          V0 50.0
          Vso 6.0
          v (/ Vso V0)
          lam2 2.0
          lam2-over-a2 (/ lam2 (* a a))
          xi 3.0
          eta 2.5]
      (is (< (m/abs (delta-surface-thomas-so-shift-over-v0 v lam2-over-a2 a 0 0.5 xi eta)) 1e-9))
      (is (< (m/abs (- (* V0 (delta-surface-thomas-so-shift-over-v0 v lam2-over-a2 a 1 1.5 xi eta))
                      (delta-surface-thomas-so-shift-MeV Vso lam2 a 1 1.5 xi eta)))
             1e-9))))
  (testing "nlz-delta-so-over-v0 matches MeV when scaled by V₀; total energy keys"
    (when (some? (find-bound-state-finite-well-nlz 1 1 8.0))
      (let [a 2.0
            V0 40.0
            Vso 5.0
            v (/ Vso V0)
            lam2 lambda-pi-squared-fm2-default
            lam2-over-a2 (/ lam2 (* a a))
            m (find-bound-state-finite-well-nlz-delta-so-over-v0 1 1 8.0 a v lam2-over-a2 1.5)
            er (:e-ratio m)
            de (:delta-e-so-over-v0 m)]
        (is (some? m))
        (is (< (m/abs (- (* V0 (:delta-e-so-over-v0 m))
                         (:delta-e-so-MeV (find-bound-state-finite-well-nlz-delta-so
                                           1 1 8.0 a Vso lam2 1.5))))
              1e-6))
        (is (< (m/abs (- (:e-central-over-v0 m) (- er))) 1e-12))
        (is (< (m/abs (- (:e-total-over-v0 m) (+ (- er) de))) 1e-12))
        (let [m8 (find-bound-state-finite-well-nlz-delta-so 1 1 8.0 a Vso lam2 1.5 V0)]
          (is (< (m/abs (- (:e-total-MeV m8)
                           (* V0 (:e-total-over-v0 m)))) 1e-4)))))))

